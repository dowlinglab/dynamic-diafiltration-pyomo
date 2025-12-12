#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

authors: Keshav Kasturi Rangan, Alexander Dowling
    
xlsx_experiment_reader.py

Functional-only parser for diafiltration Excel (.xlsx/.xls/.xlsm/.xlsb) experimental files.

Main entry point
----------------
    variables = extract_excel_variables(
        path,
        measure_sheet=None,
        index_sheet=None,
    )

Returns
-------
variables : dict[str, Any]
    A dictionary collecting all parsed quantities, including:

    # bookkeeping
    - "path"              : str
    - "measure_sheet"     : str
    - "index_sheet"       : Optional[str]

    # raw measurement data
    - "measurement_df"    : pandas.DataFrame (header-normalized sheet)

    # full-run 1D signals (length = N time points)
    - "time"              : np.ndarray
    - "mass"              : np.ndarray
    - "retentate_cond"    : np.ndarray (NaN if not present)
    - "permeate_cond"     : np.ndarray (NaN if not present)
    - "vial_marker"       : Optional[np.ndarray]

    # vial segmentation
    - "windows"                   : List[Tuple[int, int]]  (index ranges [start, end])
    - "per_vial_time"             : List[np.ndarray]
    - "per_vial_mass"             : List[np.ndarray]       (re-baselined, spike-filtered)
    - "per_vial_retentate_cond"   : List[np.ndarray]
    - "per_vial_permeate_cond"    : List[np.ndarray]

    # metadata / configuration
    - "metadata"          : dict[str, Any]   (free-text header info, rpm, datapoints, etc.)
    - "config_pairs"      : dict[str, Any]   (key–value from config-like sheets)
    - "data_config"       : dict[str, Any]   (PyDataConfig-like fields, including "extras")
"""

from __future__ import annotations

import os
import re
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# Low-level Excel measurement sheet handling
# =============================================================================

def _open_excel(path: str) -> pd.ExcelFile:
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
    ext = os.path.splitext(path)[1].lower()
    if ext not in {".xlsx", ".xls", ".xlsm", ".xlsb"}:
        raise ValueError(f"Not an Excel file: {path}")
    return pd.ExcelFile(path)


def _detect_header_row_generic(raw: pd.DataFrame, max_rows: int = 25) -> Optional[int]:
    """
    Scan the top rows of a sheet to find the header row containing both 'time'
    and ('mass' or 'cond') anywhere in the row.
    """
    nrows = min(max_rows, len(raw))
    for i in range(nrows):
        row = raw.iloc[i].astype(str).str.lower().tolist()
        has_time = any("time" in s for s in row)
        has_mass_or_cond = any(("mass" in s) or ("cond" in s) for s in row)
        if has_time and has_mass_or_cond:
            return i
    return None


def _detect_measurement_sheet(xl: pd.ExcelFile) -> Optional[str]:
    """
    Heuristic: first sheet whose top ~10 rows contain a header row with
    Time + Mass/Cond.
    """
    for sh in xl.sheet_names:
        df = xl.parse(sh, header=None, nrows=10)
        if _detect_header_row_generic(df, max_rows=10) is not None:
            return sh
    return None


def _detect_index_sheet(xl: pd.ExcelFile) -> Optional[str]:
    """
    Heuristic: sheet named like 'Index' or 'Summary' (case-insensitive).
    """
    for sh in xl.sheet_names:
        low = sh.lower()
        if "index" in low or "summary" in low:
            return sh
    return None


def _first_col_like(df: pd.DataFrame, needle: str) -> Optional[str]:
    """
    Return the name of the first column whose string representation contains
    `needle` (case-insensitive), or None if not found.
    """
    n = needle.lower()
    for c in df.columns:
        if n in str(c).lower():
            return c
    return None


def _first_col_from_list(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    """
    Return the first column whose name matches any candidate (case-insensitive),
    first via exact case-insensitive match, then via substring match.
    """
    lowers = {str(c).lower(): c for c in df.columns}
    for cand in candidates:
        cl = cand.lower()
        if cl in lowers:
            return lowers[cl]
    for cand in candidates:
        cl = cand.lower()
        for c in df.columns:
            if cl in str(c).lower():
                return c
    return None


def _read_measurement_df(
    xl: pd.ExcelFile,
    measure_sheet: Optional[str],
) -> Tuple[str, pd.DataFrame]:
    """
    Resolve the measurement sheet and return a normalized DataFrame:

        - inferred header row
        - header row used as df.columns
        - data rows below header
        - index reset to 0..N-1
    """
    if measure_sheet is not None:
        if measure_sheet not in xl.sheet_names:
            raise ValueError(f"Measurement sheet '{measure_sheet}' not in workbook.")
        sheet = measure_sheet
    else:
        sheet = _detect_measurement_sheet(xl)
        if sheet is None:
            raise ValueError("No measurement-like sheet (Time + Mass/Cond) found.")
    raw = xl.parse(sheet, header=None)
    hdr = _detect_header_row_generic(raw, max_rows=25)
    if hdr is None:
        raise ValueError(f"Could not detect header row in sheet '{sheet}'.")
    df = raw.iloc[hdr + 1 :].copy()
    df.columns = raw.iloc[hdr].astype(str).str.strip().tolist()
    df = df.reset_index(drop=True)
    return sheet, df


def _extract_time_mass_cond_vial(df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """
    From the normalized measurement df, extract:

        time, mass, retentate_cond, permeate_cond, vial_marker

    retentate/permeate may be all-NaN arrays if missing.
    vial_marker may be None.
    """
    # time
    tcol = _first_col_like(df, "time")
    if tcol is None:
        raise ValueError("Time column not found in measurement sheet.")
    time = pd.to_numeric(df[tcol], errors="coerce").to_numpy()

    # mass
    mcol = _first_col_from_list(
        df,
        ["Mass", "mass", "M (g)", "Permeate Mass", "Collection mass", "Weight"],
    )
    if mcol is None:
        raise ValueError("Mass column not found in measurement sheet.")
    mass = pd.to_numeric(df[mcol], errors="coerce").to_numpy()

    # conductivity
    rcol = _first_col_from_list(
        df,
        ["Retentate Cond", "Retentate_Cond", "RetentateCond", "Ret_Cond"],
    )
    pcol = _first_col_from_list(
        df,
        ["Permeate Cond", "Permeate_Cond", "PermeateCond", "Perm_Cond"],
    )
    rcon = (
        pd.to_numeric(df[rcol], errors="coerce").to_numpy()
        if rcol is not None
        else np.full_like(time, np.nan, dtype=float)
    )
    pcon = (
        pd.to_numeric(df[pcol], errors="coerce").to_numpy()
        if pcol is not None
        else np.full_like(time, np.nan, dtype=float)
    )

    # vial marker (e.g. 'Vial', 'Vial swap', 'VialSwap', etc.)
    vcol = _first_col_like(df, "vial") or _first_col_like(df, "vial swap")
    vial = (
        pd.to_numeric(df[vcol], errors="coerce").to_numpy()
        if vcol is not None
        else None
    )

    return time, mass, rcon, pcon, vial


# =============================================================================
# Metadata / config extraction (mirrors v8)
# =============================================================================

def _extract_excel_metadata(path: str, measure_sheet: Optional[str]) -> Dict[str, Any]:
    """
    Extract human-readable metadata from the rows ABOVE the header row
    in an Excel measurement sheet.

    Returns keys:
        - excel_metadata_rows: List[str]
        - feed_description: str (optional)
        - diafiltrate_description: str (optional)
        - notes: List[str] (optional)
        - stirring_rpm: int (optional)
        - datapoints: int (optional)
    """
    xl = _open_excel(path)
    sheet = measure_sheet or xl.sheet_names[0]
    raw = xl.parse(sheet, header=None)

    hdr = _detect_header_row_generic(raw, max_rows=25)
    if hdr is None:
        meta_block = raw.iloc[:10]
    else:
        meta_block = raw.iloc[:hdr]

    metadata_rows: List[str] = []
    extras: Dict[str, Any] = {}

    for _, row in meta_block.iterrows():
        tokens: List[str] = []
        for v in row:
            if pd.isna(v):
                continue
            s = str(v).strip()
            if not s:
                continue
            if s.lower().startswith("unnamed:"):
                continue
            tokens.append(s)
        if tokens:
            text = " ".join(tokens)
            metadata_rows.append(text)

    extras["excel_metadata_rows"] = metadata_rows

    notes: List[str] = []
    for txt in metadata_rows:
        low = txt.lower()

        if "feed" in low:
            extras.setdefault("feed_description", txt)
        if "diafiltrate" in low or "dialysate" in low:
            extras.setdefault("diafiltrate_description", txt)

        if "note" in low:
            notes.append(txt)

        m = re.search(r"(\d+)\s*rpm", low)
        if m:
            try:
                extras["stirring_rpm"] = int(m.group(1))
            except ValueError:
                pass

        if "datapoints" in low:
            m2 = re.search(r"datapoints\s*[:\-]\s*(\d+)", low)
            if m2:
                try:
                    extras["datapoints"] = int(m2.group(1))
                except ValueError:
                    pass

    if notes:
        extras["notes"] = notes

    return extras


def _extract_excel_config_pairs(
    path: str,
    measure_sheet: Optional[str],
    index_sheet: Optional[str],
) -> Dict[str, Any]:
    """
    Look for a 'config-like' sheet (Config, Setup, Parameters, etc.) and
    parse it as key–value pairs. Pattern:

        Row: [ key, value, ... ]

    The first non-empty cell is the key; the second non-empty cell is the value.
    """
    xl = _open_excel(path)
    candidates: List[str] = []
    for sh in xl.sheet_names:
        if sh == measure_sheet or sh == index_sheet:
            continue
        low = sh.lower()
        if any(k in low for k in ["config", "setup", "parameter", "params"]):
            candidates.append(sh)

    if not candidates:
        return {}

    config_pairs: Dict[str, Any] = {}
    for sh in candidates:
        df = xl.parse(sh, header=None)
        for _, row in df.iterrows():
            tokens = [str(v).strip() for v in row if pd.notna(v) and str(v).strip()]
            if len(tokens) < 2:
                continue
            key = tokens[0]
            val_raw = tokens[1]
            try:
                val = float(val_raw)
            except ValueError:
                val = val_raw
            k_norm = key.strip()
            if k_norm in config_pairs:
                old = config_pairs[k_norm]
                if isinstance(old, list):
                    old.append(val)
                else:
                    config_pairs[k_norm] = [old, val]
            else:
                config_pairs[k_norm] = val

    return config_pairs


def _build_data_config_dict_from_excel(
    path: str,
    measure_sheet: Optional[str],
    index_sheet: Optional[str],
) -> Dict[str, Any]:
    """
    Construct a PyDataConfig-like dict from:
      - free-text metadata rows (above header),
      - config-sheet key–value pairs.

    Returns a dict containing:
        - core numeric / array fields when recognizable
        - "extras": all remaining fields + metadata
    """
    meta_extras = _extract_excel_metadata(path, measure_sheet)
    cfg_pairs = _extract_excel_config_pairs(path, measure_sheet, index_sheet)

    numeric_fields = {
        "M_F0", "M_O", "nc", "ni", "delP", "Temp", "Am", "rho",
        "Lp0", "sigma0", "n", "nr", "n_extra",
    }
    array_fields = {"C_F0", "C_D", "B0", "theta0"}

    cfg_kwargs: Dict[str, Any] = {}
    extras: Dict[str, Any] = {}

    for key, val in cfg_pairs.items():
        key_norm = key.strip()
        if key_norm in numeric_fields:
            try:
                cfg_kwargs[key_norm] = float(val)
            except Exception:
                cfg_kwargs[key_norm] = val
        elif key_norm in array_fields:
            if isinstance(val, str):
                parts = re.split(r"[,\s]+", val.strip())
                parts = [p for p in parts if p]
                try:
                    arr = np.array([float(p) for p in parts], dtype=float)
                    cfg_kwargs[key_norm] = arr
                except Exception:
                    cfg_kwargs[key_norm] = val
            else:
                cfg_kwargs[key_norm] = np.atleast_1d(val)
        else:
            extras[key_norm] = val

    for k, v in meta_extras.items():
        if k in extras:
            old = extras[k]
            if isinstance(old, list):
                if isinstance(v, list):
                    old.extend(v)
                else:
                    old.append(v)
            else:
                extras[k] = [old, v]
        else:
            extras[k] = v

    cfg_kwargs["extras"] = extras
    return cfg_kwargs


# =============================================================================
# Vial window inference (functional version of IndexSheetStrategy / VialMarker)
# =============================================================================

def _infer_windows_from_index_sheet(
    path: str,
    index_sheet: Optional[str],
    time: np.ndarray,
) -> List[Tuple[int, int]]:
    if index_sheet is None:
        return []
    xl = _open_excel(path)
    if index_sheet not in xl.sheet_names:
        return []

    try:
        idx_df = xl.parse(index_sheet)
    except Exception:
        return []

    def first_col_like(df: pd.DataFrame, needle: str) -> Optional[str]:
        n = needle.lower()
        for c in df.columns:
            if n in str(c).lower():
                return c
        return None

    start_col = first_col_like(idx_df, "start")
    end_col = first_col_like(idx_df, "end")
    if start_col is None or end_col is None:
        return []

    starts = pd.to_numeric(idx_df[start_col], errors="coerce").dropna().astype(int).tolist()
    ends = pd.to_numeric(idx_df[end_col], errors="coerce").dropna().astype(int).tolist()

    windows: List[Tuple[int, int]] = []
    n = len(time)
    for s, e in zip(starts, ends):
        if 0 <= s < n and 0 <= e < n and s <= e:
            windows.append((s, e))
    return windows


def _infer_windows_from_vial_marker(vial_marker: Optional[np.ndarray]) -> List[Tuple[int, int]]:
    if vial_marker is None:
        return []
    v_arr = np.asarray(vial_marker).squeeze()
    if v_arr.size == 0:
        return []

    try:
        vnum = pd.to_numeric(pd.Series(v_arr), errors="coerce").fillna(-1).to_numpy()
    except Exception:
        v_str = np.array(list(map(str, v_arr)))
        vnum = (v_str != np.roll(v_str, 1)).astype(int).cumsum()

    boundaries: List[int] = [0]
    for i in range(1, len(vnum)):
        if vnum[i] != vnum[i - 1]:
            boundaries.append(i)
    boundaries.append(len(vnum))

    windows: List[Tuple[int, int]] = []
    for s, e in zip(boundaries[:-1], boundaries[1:]):
        if e - s >= 3:
            windows.append((int(s), int(e - 1)))
    return windows


def _infer_windows_default(time: np.ndarray) -> List[Tuple[int, int]]:
    if len(time) == 0:
        return []
    return [(0, len(time) - 1)]


def _infer_vial_windows(
    path: str,
    index_sheet: Optional[str],
    time: np.ndarray,
    vial_marker: Optional[np.ndarray],
) -> List[Tuple[int, int]]:
    """
    Apply strategies in order:
      1) index sheet
      2) vial marker
      3) single-window fallback
    """
    w = _infer_windows_from_index_sheet(path, index_sheet, time)
    if w:
        return w
    w = _infer_windows_from_vial_marker(vial_marker)
    if w:
        return w
    return _infer_windows_default(time)


# =============================================================================
# Per-vial slicing + mass cleaning
# =============================================================================

def _slice_per_vial_signals(
    time: np.ndarray,
    mass: np.ndarray,
    retentate_cond: np.ndarray,
    permeate_cond: np.ndarray,
    windows: List[Tuple[int, int]],
) -> Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray], List[np.ndarray]]:
    """
    Slice full-run arrays into per-vial windows and apply:

      - per-vial mass re-baselining
      - fixed 1.2 g spike filter on Δmass
    """
    per_vial_t: List[np.ndarray] = []
    per_vial_m: List[np.ndarray] = []
    per_vial_r: List[np.ndarray] = []
    per_vial_p: List[np.ndarray] = []

    for s, e in windows:
        sl = slice(s, e + 1)
        t = time[sl]
        m_raw = mass[sl]

        if m_raw.size and not np.isnan(m_raw[0]):
            m = m_raw - m_raw[0]
        else:
            m = m_raw.copy()

        m = np.where(m > 1.2, np.nan, m)

        r = retentate_cond[sl]
        p = permeate_cond[sl]

        per_vial_t.append(t)
        per_vial_m.append(m)
        per_vial_r.append(r)
        per_vial_p.append(p)

    return per_vial_t, per_vial_m, per_vial_r, per_vial_p


# =============================================================================
# Public API
# =============================================================================

def extract_excel_variables(
    path: str,
    measure_sheet: Optional[str] = None,
    index_sheet: Optional[str] = None,
) -> Dict[str, Any]:
    """
    High-level functional entry point.

    Parameters
    ----------
    path : str
        Path to Excel file.
    measure_sheet : Optional[str]
        Measurement sheet name (each sheet corresponds to one experiment). If None, auto-detect.
    index_sheet : Optional[str]
        Index sheet name override (if None, auto-detect).

    Returns
    -------
    variables : dict[str, Any]
        Dictionary collecting all parsed variables (see module docstring).
    """
    xl = _open_excel(path)

    # Resolve measurement & index sheets
    if measure_sheet is None:
        measure_sheet, df = _read_measurement_df(xl, measure_sheet=None)
    else:
        measure_sheet, df = _read_measurement_df(xl, measure_sheet=measure_sheet)

    if index_sheet is None:
        index_sheet = _detect_index_sheet(xl)

    # Extract core signals
    time, mass, rcon, pcon, vial = _extract_time_mass_cond_vial(df)

    # Vial windows
    windows = _infer_vial_windows(path, index_sheet, time, vial)

    # Per-vial slices
    per_t, per_m, per_r, per_p = _slice_per_vial_signals(time, mass, rcon, pcon, windows)

    # Metadata / configuration
    metadata = _extract_excel_metadata(path, measure_sheet)
    config_pairs = _extract_excel_config_pairs(path, measure_sheet, index_sheet)
    data_config = _build_data_config_dict_from_excel(path, measure_sheet, index_sheet)

    variables: Dict[str, Any] = {
# bookkeeping
        "path": path,
        "measure_sheet": measure_sheet,
        "index_sheet": index_sheet,
        # experiment-level metadata (when available in Excel)
        "dataset": data_config.get("dataset", None),
        "filename": data_config.get("filename", os.path.basename(path)),
        "mode": data_config.get("mode", None),
        "continuous_cF": data_config.get("continuous_cF", None),
        "n_obj": data_config.get("n_obj", None),
        # measurement data
        "measurement_df": df,
        # configuration
        "data_config": data_config,
        # full-run signals
        "time": time,
        "mass": mass,
        "cF_exp": rcon,  # retentate signal (uS/cm) from Excel
        "cV_avg": pcon,  # permeate signal (uS/cm) from Excel,
"vial_marker": vial,
        # vial segmentation
        "vial_ranges": windows,
        "windows": windows,  # deprecated alias (use vial_ranges)
        "per_vial_time": per_t,
        "per_vial_mass": per_m,
        # vial objects (MATLAB-like data_raw equivalent)
        "vials": [
            {
                "number": i + 1,
                "time": time[s : e + 1],
                "mass": mass[s : e + 1],
                "cF_exp": rcon[s : e + 1],
                "cV_avg": pcon[s : e + 1],
                "nr": int(e - s + 1),
                "extras": {},
            }
            for i, (s, e) in enumerate(windows)
        ],
        "per_vial_cF_exp": per_r,
        "per_vial_cV_avg": per_p,
        # metadata / config
        "metadata": metadata,
        "config_pairs": config_pairs,
        "data_config": data_config,
    }

    return variables


# Optional: simple CLI for inspection

# =============================================================================
# Plotting helpers (match MATLAB-style test script)
# =============================================================================

def plot_mass_and_conductivity(vars_: Dict[str, Any]) -> None:
    """Plot mass (solid) and cF_exp/cV_avg (solid/dashed) with vial-change markers.

    Uses:
        - vars_["time"], vars_["mass"], vars_["cF_exp"], vars_["cV_avg"]
        - vars_["vial_ranges"]
    """
    import matplotlib.pyplot as plt
    import numpy as np

    time = np.asarray(vars_["time"]).squeeze()
    mass = np.asarray(vars_["mass"]).squeeze()
    cF = np.asarray(vars_["cF_exp"]).squeeze()
    cV = np.asarray(vars_["cV_avg"]).squeeze()
    vial_ranges = list(vars_.get("vial_ranges", []))

    def _add_vial_change_lines(ax):
        for (s, _) in vial_ranges[1:]:
            if 0 <= s < len(time):
                ax.axvline(time[s], linestyle="--", linewidth=1)

    # Mass
    fig, ax = plt.subplots()
    ax.plot(time, mass, linestyle="-", label="mass")
    _add_vial_change_lines(ax)
    ax.set_xlabel("time")
    ax.set_ylabel("mass")
    ax.set_title("Mass (raw) with vial changes")
    ax.legend()

    # Conductivity / concentration
    fig, ax = plt.subplots()
    ax.plot(time, cF, linestyle="-", label="cF_exp (retentate)")
    ax.plot(time, cV, linestyle="--", label="cV_avg (permeate)")
    _add_vial_change_lines(ax)
    ax.set_xlabel("time")
    ax.set_ylabel("signal")
    ax.set_title("cF_exp / cV_avg with vial changes")
    ax.legend()

    plt.show()

def _cli():
    import argparse

    parser = argparse.ArgumentParser(description="Functional Excel diafiltration data extractor.")
    parser.add_argument("path", help="Path to Excel file (.xlsx/.xls/.xlsm/.xlsb)")
    parser.add_argument("--measure-sheet", type=str, default=None, help="Measurement sheet override")
    parser.add_argument("--index-sheet", type=str, default=None, help="Index sheet override")
    args = parser.parse_args()

    vars_dict = extract_excel_variables(args.path, args.measure_sheet, args.index_sheet)

    print(f"\n[info] Parsed variables from {args.path}")
    print("Keys:")
    for k in sorted(vars_dict.keys()):
        v = vars_dict[k]
        if isinstance(v, np.ndarray):
            print(f"  {k:30s} : ndarray shape={v.shape}")
        elif isinstance(v, list) and v and isinstance(v[0], np.ndarray):
            shapes = [x.shape for x in v]
            print(f"  {k:30s} : list of ndarrays, shapes={shapes}")
        else:
            print(f"  {k:30s} : type={type(v).__name__}")

if __name__ == "__main__":
    _cli()
