#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

authors: Keshav Kasturi Rangan, Alexander Dowling

experiment_dataload_mat_v1.py

Functional-only parser for diafiltration MATLAB (.mat) experimental files.

Design goals
------------
- Keep a *clean functional pipeline* (small helpers, explicit steps).
- Preserve *MATLAB-native variable names* as the public API (no aliases).
- Preserve *all MATLAB fields losslessly* via nested "_matlab" snapshots.
- Represent vial segmentation explicitly using "vial_ranges" (not "windows").

Supported MAT formats
---------------------
1) data_stru format (preferred, produced by load_data.m):
    - data_stru.dataset, data_stru.filename, data_stru.mode, data_stru.continuous_cF, data_stru.n_obj
    - data_stru.data_config (struct with many fields)
    - data_stru.data_raw(i).time, mass, cF_exp, cV_avg, nr, ...

2) flat-array format (fallback):
    - Time, Mass
    - optional cF_exp / cV_avg (or RetentateCond / PermeateCond variants)
    - optional Vial or VialSwap (or variants) for vial segmentation

Main entry point
----------------
    variables = extract_mat_variables(path)

Returns
-------
variables : dict[str, Any]
    Dictionary collecting all parsed quantities in a Python-friendly structure
    while preserving MATLAB-native names:

    # bookkeeping
    - "path"              : str
    - "format"            : str  ("data_stru" or "flat_arrays")

    # data_stru core fields (None in flat_arrays)
    - "dataset"           : Optional[float]
    - "filename"          : Optional[str]
    - "mode"              : Optional[str]
    - "continuous_cF"     : Optional[bool]
    - "n_obj"             : Optional[int]

    # config (MATLAB-like mapping + lossless snapshot)
    - "data_config"       : dict[str, Any]
        - includes "extras" for non-canonical fields
        - includes "_matlab" with all raw MATLAB fields
        - includes "n" and "nr" (computed if missing)

    # vial objects (rich per-vial dicts, even for flat_arrays we provide a single vial)
    - "vials"             : List[dict[str, Any]]
        Each vial dict uses MATLAB field names:
          - "number"  : int
          - "time"    : np.ndarray
          - "mass"    : np.ndarray
          - "cV_avg"  : np.ndarray  (permeate concentration measurement)
          - "cF_exp"  : np.ndarray  (final retentate concentration measurement)
          - "nr"      : int | Any
          - "extras"  : dict[str, Any] (non-core fields)
          - "_matlab" : dict[str, Any] (lossless snapshot of all MATLAB fields in data_raw(i))

    # full-run arrays (always provided by concatenating vials)
    - "time"              : np.ndarray
    - "mass"              : np.ndarray
    - "cV_avg"            : np.ndarray
    - "cF_exp"            : np.ndarray
    - "vial_marker"       : Optional[np.ndarray]  (present for data_stru concatenation; optional for flat_arrays)

    # vial segmentation / per-vial slices (mirrors Excel script structure but uses "vial" terminology)
    - "vial_ranges"             : List[Tuple[int, int]]
    - "per_vial_time"           : List[np.ndarray]
    - "per_vial_mass"           : List[np.ndarray]   (re-baselined + spike-filtered)
    - "per_vial_cV_avg"         : List[np.ndarray]
    - "per_vial_cF_exp"         : List[np.ndarray]

Notes
-----
- Vial range inference order mirrors the unified loader idea:
    1) If data_stru: vial_ranges come directly from vial boundaries (robust, canonical).
    2) Else if vial_marker exists: infer vial_ranges from changes in vial marker.
    3) Else fallback: single vial_range.
"""

from __future__ import annotations

import os
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scipy.io as sio


# =============================================================================
# MAT loading utilities
# =============================================================================

def _load_mat_scipy(path: str) -> Dict[str, Any]:
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
    if os.path.splitext(path)[1].lower() != ".mat":
        raise ValueError(f"Not a .mat file: {path}")
    return sio.loadmat(path, squeeze_me=True, struct_as_record=False)


def _mat_struct_to_dict(struct_obj: Any) -> Dict[str, Any]:
    """Convert a scipy-loaded MATLAB struct to a dict of its fields (one level)."""
    out: Dict[str, Any] = {}
    if struct_obj is None:
        return out
    if hasattr(struct_obj, "_fieldnames"):
        for fn in struct_obj._fieldnames:
            out[fn] = getattr(struct_obj, fn)
    return out


def _as_1d_float_array(x: Any) -> np.ndarray:
    if x is None:
        return np.array([], dtype=float)
    arr = np.asarray(x).squeeze()
    # Strings -> attempt numeric; if not possible return empty
    if arr.dtype.kind in {"U", "S", "O"}:
        try:
            arr = pd.to_numeric(pd.Series(arr.ravel()), errors="coerce").to_numpy()
        except Exception:
            return np.array([], dtype=float)
    return np.asarray(arr, dtype=float).squeeze()


def _grab_case_insensitive(mat: Dict[str, Any], *names: str) -> Optional[Any]:
    """Look up keys in a scipy loadmat dict ignoring case."""
    key_map = {k.lower(): k for k in mat.keys()}
    for n in names:
        k = key_map.get(n.lower())
        if k is not None:
            return mat[k]
    return None


# =============================================================================
# data_stru parsing (mirrors load_data.m / PyDataStru.from_mat)
# =============================================================================

def _parse_data_config(ds: Any) -> Dict[str, Any]:
    """Return a dict for data_config with canonical fields + extras + _matlab snapshot."""
    dc = getattr(ds, "data_config", None)

    cfg: Dict[str, Any] = {}
    extras: Dict[str, Any] = {}
    raw_all: Dict[str, Any] = _mat_struct_to_dict(dc)

    if dc is not None and hasattr(dc, "_fieldnames"):
        for fn in dc._fieldnames:
            val = getattr(dc, fn)
            if isinstance(val, (list, np.ndarray)):
                val = np.asarray(val).squeeze()

            # Keep known fields at top-level; everything else in extras
            if fn in {
                "M_F0", "C_F0", "M_O", "C_D",
                "nc", "namec", "ni",
                "delP", "Temp", "Am", "rho",
                "Lp0", "B0", "sigma0", "theta0",
                "n", "nr",
                "n_extra",  # legacy-ish
            }:
                cfg[fn] = val
            else:
                extras[fn] = val

    cfg["extras"] = extras
    cfg["_matlab"] = raw_all
    return cfg


def _parse_data_raw_vials(ds: Any) -> List[Dict[str, Any]]:
    """Return list of per-vial dicts using MATLAB field names + lossless _matlab snapshot."""
    dr = getattr(ds, "data_raw", None)
    if dr is None:
        return []

    dr_iter = dr.ravel() if isinstance(dr, np.ndarray) else [dr]

    vials: List[Dict[str, Any]] = []
    for i, v in enumerate(dr_iter, start=1):
        # Core fields
        t = _as_1d_float_array(getattr(v, "time", []))
        m = _as_1d_float_array(getattr(v, "mass", []))
        cV = getattr(v, "cV_avg", None)
        cF = getattr(v, "cF_exp", None)
        cV_arr = _as_1d_float_array(cV) if cV is not None else np.array([], dtype=float)
        cF_arr = _as_1d_float_array(cF) if cF is not None else np.array([], dtype=float)

        nr = getattr(v, "nr", None)
        number = getattr(v, "number", i)

        # Lossless snapshot of all MATLAB fields in this vial struct
        raw_all = _mat_struct_to_dict(v)

        # Extras: fields beyond the canonical load_data.m set
        vial_extras: Dict[str, Any] = {}
        if hasattr(v, "_fieldnames"):
            for fn in v._fieldnames:
                if fn not in {"time", "mass", "cF_exp", "cV_avg", "nr", "number"}:
                    vial_extras[fn] = getattr(v, fn)

        vials.append({
            "number": int(number),
            "time": t,
            "mass": m,
            "cV_avg": cV_arr,
            "cF_exp": cF_arr,
            "nr": int(nr) if isinstance(nr, (int, np.integer)) else nr,
            "extras": vial_extras,
            "_matlab": raw_all,
        })

    return vials


def _concatenate_vials_to_full_run(
    vials: List[Dict[str, Any]]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, List[Tuple[int, int]]]:
    """Concatenate vials to full-run arrays and return vial_marker + vial_ranges."""
    times: List[np.ndarray] = []
    masses: List[np.ndarray] = []
    cV_list: List[np.ndarray] = []
    cF_list: List[np.ndarray] = []
    vmarks: List[np.ndarray] = []
    vial_ranges: List[Tuple[int, int]] = []

    cursor = 0
    for idx, vial in enumerate(vials, start=1):
        t = np.asarray(vial["time"], dtype=float).squeeze()
        m = np.asarray(vial["mass"], dtype=float).squeeze()

        def _align(arr: np.ndarray, n: int) -> np.ndarray:
            a = np.asarray(arr, dtype=float).squeeze()
            if a.size == 0:
                return np.full(n, np.nan, dtype=float)
            if a.size >= n:
                return a[:n]
            return np.pad(a, (0, n - a.size), constant_values=np.nan)

        cV = _align(vial.get("cV_avg", np.array([], dtype=float)), t.size)
        cF = _align(vial.get("cF_exp", np.array([], dtype=float)), t.size)

        times.append(t)
        masses.append(m)
        cV_list.append(cV)
        cF_list.append(cF)
        vmarks.append(np.full(t.shape, idx, dtype=float))

        if t.size > 0:
            start = cursor
            end = cursor + t.size - 1
            vial_ranges.append((start, end))
            cursor += t.size

    if not times:
        empty = np.array([], dtype=float)
        return empty, empty, empty, empty, empty, []

    time = np.concatenate(times)
    mass = np.concatenate(masses)
    cV_avg = np.concatenate(cV_list)
    cF_exp = np.concatenate(cF_list)
    vial_marker = np.concatenate(vmarks)

    return time, mass, cV_avg, cF_exp, vial_marker, vial_ranges


# =============================================================================
# Vial range inference for flat arrays
# =============================================================================

def _infer_vial_ranges_from_vial_marker(vial_marker: Optional[np.ndarray]) -> List[Tuple[int, int]]:
    if vial_marker is None:
        return []
    v_arr = np.asarray(vial_marker).squeeze()
    if v_arr.size == 0:
        return []

    # Try numeric; fallback to change-points in string representation
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

    vial_ranges: List[Tuple[int, int]] = []
    for s, e in zip(boundaries[:-1], boundaries[1:]):
        if e - s >= 3:  # keep only non-trivial segments
            vial_ranges.append((int(s), int(e - 1)))
    return vial_ranges


def _infer_vial_ranges_default(time: np.ndarray) -> List[Tuple[int, int]]:
    if len(time) == 0:
        return []
    return [(0, len(time) - 1)]


def _infer_vial_ranges(
    format_name: str,
    time: np.ndarray,
    vial_marker: Optional[np.ndarray],
    data_stru_vial_ranges: Optional[List[Tuple[int, int]]] = None,
) -> List[Tuple[int, int]]:
    """Priority: canonical vial boundaries -> marker changes -> single segment."""
    if format_name == "data_stru" and data_stru_vial_ranges:
        return data_stru_vial_ranges

    vr = _infer_vial_ranges_from_vial_marker(vial_marker)
    if vr:
        return vr
    return _infer_vial_ranges_default(time)


# =============================================================================
# Per-vial slicing + mass cleaning (same logic as Excel script)
# =============================================================================

def _slice_per_vial_signals(
    time: np.ndarray,
    mass: np.ndarray,
    cV_avg: np.ndarray,
    cF_exp: np.ndarray,
    vial_ranges: List[Tuple[int, int]],
) -> Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray], List[np.ndarray]]:
    """
    Slice full-run arrays into per-vial segments and apply:
      - per-vial mass re-baselining
      - fixed 1.2 g spike filter on Δmass
    """
    per_vial_t: List[np.ndarray] = []
    per_vial_m: List[np.ndarray] = []
    per_vial_cV: List[np.ndarray] = []
    per_vial_cF: List[np.ndarray] = []

    for s, e in vial_ranges:
        sl = slice(s, e + 1)
        t = time[sl]
        m_raw = mass[sl]

        if m_raw.size and not np.isnan(m_raw[0]):
            m = m_raw - m_raw[0]
        else:
            m = m_raw.copy()

        # spike filter
        m = np.where(m > 1.2, np.nan, m)

        per_vial_t.append(t)
        per_vial_m.append(m)
        per_vial_cV.append(cV_avg[sl])
        per_vial_cF.append(cF_exp[sl])

    return per_vial_t, per_vial_m, per_vial_cV, per_vial_cF


# =============================================================================
# Flat-array parsing
# =============================================================================

def _parse_flat_arrays(mat: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """
    Parse a generic MAT file with flat arrays.
    Requires Time and Mass.
    cV_avg/cF_exp are optional (we try common variants).
    vial_marker optional.
    """
    time = _grab_case_insensitive(mat, "Time", "time", "t")
    mass = _grab_case_insensitive(mat, "Mass", "mass", "m")
    if time is None or mass is None:
        raise ValueError("MAT file missing required arrays: Time and Mass.")

    time = _as_1d_float_array(time)
    mass = _as_1d_float_array(mass)

    # Prefer MATLAB-native names, but accept common legacy synonyms
    cF = _grab_case_insensitive(mat, "cF_exp", "CF_exp", "RetentateCond", "Retentate_Cond", "retentate_cond", "ret_cond")
    cV = _grab_case_insensitive(mat, "cV_avg", "CV_avg", "PermeateCond", "Permeate_Cond", "permeate_cond", "perm_cond")

    cF_arr = _as_1d_float_array(cF) if cF is not None else np.full_like(time, np.nan, dtype=float)
    cV_arr = _as_1d_float_array(cV) if cV is not None else np.full_like(time, np.nan, dtype=float)

    vial = _grab_case_insensitive(mat, "Vial", "VialSwap", "vial", "vial_swap", "vialswap")
    vial_arr = _as_1d_float_array(vial) if vial is not None else None

    # Fit lengths to time
    n = time.size

    def _fit(x):
        if x is None:
            return None
        x = np.asarray(x).squeeze()
        if x.size >= n:
            return x[:n]
        return np.pad(x, (0, n - x.size), constant_values=np.nan)

    cF_arr = _fit(cF_arr)
    cV_arr = _fit(cV_arr)
    vial_arr = _fit(vial_arr)

    return time, mass, cV_arr, cF_arr, vial_arr


# =============================================================================
# Public API
# =============================================================================

def extract_mat_variables(path: str) -> Dict[str, Any]:
    """High-level functional entry point for MAT files."""
    mat = _load_mat_scipy(path)

    # ----------------------------------------------------------
    # Preferred: data_stru format
    # ----------------------------------------------------------
    if "data_stru" in mat:
        ds = mat["data_stru"]

        dataset = getattr(ds, "dataset", None)
        try:
            dataset = float(dataset) if dataset is not None else None
        except Exception:
            dataset = None

        filename = getattr(ds, "filename", None)
        filename = str(filename) if filename is not None else None

        mode = getattr(ds, "mode", None)
        mode = str(mode) if mode is not None else None

        continuous_cF = getattr(ds, "continuous_cF", None)
        if continuous_cF is not None:
            try:
                continuous_cF = bool(continuous_cF)
            except Exception:
                continuous_cF = None

        n_obj = getattr(ds, "n_obj", None)
        try:
            n_obj = int(n_obj) if n_obj is not None else None
        except Exception:
            n_obj = None

        data_config = _parse_data_config(ds)
        vials = _parse_data_raw_vials(ds)

        # Concatenate to full-run arrays + canonical vial_ranges + vial_marker
        time, mass, cV_avg, cF_exp, vial_marker, ds_vial_ranges = _concatenate_vials_to_full_run(vials)

        format_name = "data_stru"
        vial_ranges = _infer_vial_ranges(format_name, time, vial_marker, data_stru_vial_ranges=ds_vial_ranges)

    # ----------------------------------------------------------
    # Fallback: flat arrays
    # ----------------------------------------------------------
    else:
        dataset = None
        filename = None
        mode = None
        continuous_cF = None
        n_obj = None
        data_config = {"extras": {}, "_matlab": {}}

        time, mass, cV_avg, cF_exp, vial_marker = _parse_flat_arrays(mat)

        # Build a single vial object for completeness
        vials = [{
            "number": 1,
            "time": time,
            "mass": mass,
            "cV_avg": cV_avg,
            "cF_exp": cF_exp,
            "nr": int(np.isfinite(mass).sum()),
            "extras": {},
            "_matlab": {},
        }]

        format_name = "flat_arrays"
        vial_ranges = _infer_vial_ranges(format_name, time, vial_marker, data_stru_vial_ranges=None)

    # Per-vial slices with cleaning
    per_t, per_m, per_cV, per_cF = _slice_per_vial_signals(time, mass, cV_avg, cF_exp, vial_ranges)

    # Ensure MATLAB summary fields in data_config
    # n = number of vials
    data_config.setdefault("n", len(vials))
    # nr = total residual count (best-effort)
    if "nr" not in data_config:
        total_nr = 0
        for v in vials:
            vnr = v.get("nr", None)
            if isinstance(vnr, (int, np.integer)):
                total_nr += int(vnr)
            else:
                # fallback: count finite mass points in this vial
                total_nr += int(np.isfinite(np.asarray(v.get("mass", []))).sum())
        data_config["nr"] = int(total_nr)

    variables: Dict[str, Any] = {
        # bookkeeping
        "path": path,
        "format": format_name,

        # data_stru core (MATLAB names)
        "dataset": dataset,
        "filename": filename,
        "mode": mode,
        "continuous_cF": continuous_cF,
        "n_obj": n_obj,

        # config + vials
        "data_config": data_config,
        "vials": vials,

        # full-run signals (MATLAB names)
        "time": time,
        "mass": mass,
        "cV_avg": cV_avg,
        "cF_exp": cF_exp,
        "vial_marker": vial_marker,

        # vial segmentation
        "vial_ranges": vial_ranges,
        "per_vial_time": per_t,
        "per_vial_mass": per_m,
        "per_vial_cV_avg": per_cV,
        "per_vial_cF_exp": per_cF,
    }

    return variables


# Optional CLI for inspection
def _cli():
    import argparse

    parser = argparse.ArgumentParser(description="Functional MAT diafiltration data extractor.")
    parser.add_argument("path", help="Path to MAT file (.mat)")
    args = parser.parse_args()

    vars_dict = extract_mat_variables(args.path)

    print(f"\n[info] Parsed variables from {args.path}")
    print("format:", vars_dict["format"])
    print("keys:")
    for k in sorted(vars_dict.keys()):
        v = vars_dict[k]
        if isinstance(v, np.ndarray):
            print(f"  {k:30s} : ndarray shape={v.shape}")
        elif isinstance(v, list) and v and isinstance(v[0], np.ndarray):
            shapes = [x.shape for x in v]
            print(f"  {k:30s} : list of ndarrays, shapes={shapes}")
        elif k == "vials":
            print(f"  {k:30s} : {len(v)} vial dict(s)")
        else:
            print(f"  {k:30s} : type={type(v).__name__}")


if __name__ == "__main__":
    _cli()
