#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

authors: Keshav Kasturi Rangan, Alexander Dowling 
    
mat_experiment_reader.py

Functional-only parser for diafiltration MATLAB (.mat) experimental files.

Supports two formats:

1) data_stru format (preferred):
    - data_stru.dataset, data_stru.filename, data_stru.mode, data_stru.conductivity_cF, data_stru.n_obj
    - data_stru.data_config (struct with many fields)
    - data_stru.data_raw(i).time, mass, cF_exp, cV_avg, nr, ...

2) flat-array format (fallback):
    - Time, Mass
    - optional RetentateCond, PermeateCond
    - optional Vial or VialSwap (or variants) for windowing

Main entry point
----------------
    variables = extract_mat_variables(path)

Returns
-------
variables : dict[str, Any]
    Dictionary collecting all parsed quantities:

    # bookkeeping
    - "path"              : str
    - "format"            : str  ("data_stru" or "flat_arrays")

    # data_stru core fields (None in flat_arrays)
    - "dataset"           : Optional[float]
    - "filename"          : Optional[str]
    - "mode"              : Optional[str]
    - "conductivity_cF"   : Optional[bool]
    - "n_obj"             : Optional[int]

    # config
    - "data_config"       : dict[str, Any]   (PyDataConfig-like mapping with "extras")

    # vial objects (rich per-vial dicts, even for flat_arrays we provide a single vial)
    - "vials"             : List[dict[str, Any]]

    # full-run arrays (always provided by concatenating vials)
    - "time"              : np.ndarray
    - "mass"              : np.ndarray
    - "retentate_cond"    : np.ndarray
    - "permeate_cond"     : np.ndarray
    - "vial_marker"       : Optional[np.ndarray]  (always present for data_stru concatenation)

    # windowing / per-vial slices (mirrors Excel script structure)
    - "windows"                   : List[Tuple[int, int]]
    - "per_vial_time"             : List[np.ndarray]
    - "per_vial_mass"             : List[np.ndarray]       (re-baselined + spike-filtered)
    - "per_vial_retentate_cond"   : List[np.ndarray]
    - "per_vial_permeate_cond"    : List[np.ndarray]

Notes
-----
- Window inference order mirrors the unified loader idea:
    1) If data_stru: windows come directly from vial boundaries (robust, canonical).
    2) Else if vial_marker exists: infer windows from changes in vial marker.
    3) Else fallback: single window.
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
    """
    Convert a scipy-loaded MATLAB struct (mat_struct-like) to a dict of fields.
    Does not recursively expand nested structs beyond one level (by design);
    we keep "extras" for unknown fields.
    """
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
    # Strings -> empty
    if arr.dtype.kind in {"U", "S", "O"}:
        try:
            # allow numeric strings
            arr = pd.to_numeric(pd.Series(arr.ravel()), errors="coerce").to_numpy()
        except Exception:
            return np.array([], dtype=float)
    return np.asarray(arr, dtype=float).squeeze()


def _grab_case_insensitive(mat: Dict[str, Any], *names: str) -> Optional[Any]:
    """
    Look up keys in a scipy loadmat dict ignoring case.
    """
    key_map = {k.lower(): k for k in mat.keys()}
    for n in names:
        k = key_map.get(n.lower())
        if k is not None:
            return mat[k]
    return None


# =============================================================================
# data_stru parsing (mirrors PyDataStru.from_mat)
# =============================================================================

def _parse_data_config(ds: Any) -> Dict[str, Any]:
    """
    Returns a PyDataConfig-like dict with an "extras" dict.
    """
    dc = getattr(ds, "data_config", None)
    cfg: Dict[str, Any] = {}
    extras: Dict[str, Any] = {}

    if dc is not None and hasattr(dc, "_fieldnames"):
        for fn in dc._fieldnames:
            val = getattr(dc, fn)
            # normalize arrays
            if isinstance(val, (list, np.ndarray)):
                val = np.asarray(val).squeeze()
            # Keep canonical known fields at top-level; everything else in extras
            if fn in {
                "M_F0", "C_F0", "M_O", "C_D", "nc", "namec", "ni",
                "delP", "Temp", "Am", "rho",
                "Lp0", "B0", "sigma0", "theta0",
                "n", "nr", "n_extra",
            }:
                cfg[fn] = val
            else:
                extras[fn] = val

    cfg["extras"] = extras
    return cfg


def _parse_data_raw_vials(ds: Any) -> List[Dict[str, Any]]:
    """
    Returns list of per-vial dicts (rich), one per data_raw(i).
    """
    dr = getattr(ds, "data_raw", None)
    if dr is None:
        return []

    if isinstance(dr, np.ndarray):
        dr_iter = dr.ravel()
    else:
        dr_iter = [dr]

    vials: List[Dict[str, Any]] = []
    for i, v in enumerate(dr_iter, start=1):
        t = _as_1d_float_array(getattr(v, "time", []))
        m = _as_1d_float_array(getattr(v, "mass", []))

        cF = getattr(v, "cF_exp", None)
        cV = getattr(v, "cV_avg", None)
        cF_arr = _as_1d_float_array(cF) if cF is not None else None
        cV_arr = _as_1d_float_array(cV) if cV is not None else None

        nr = getattr(v, "nr", None)
        number = getattr(v, "number", i)

        vial_extras: Dict[str, Any] = {}
        if hasattr(v, "_fieldnames"):
            for fn in v._fieldnames:
                if fn not in {"time", "mass", "cF_exp", "cV_avg", "nr", "number"}:
                    vial_extras[fn] = getattr(v, fn)

        vials.append({
            "number": int(number),
            "time": t,
            "mass": m,
            "retentate_cond": cF_arr,   # matches your ontology (cF_exp)
            "permeate_cond": cV_arr,    # matches your ontology (cV_avg)
            "nr": int(nr) if isinstance(nr, (int, np.integer)) else nr,
            "extras": vial_extras,
        })

    return vials


def _concatenate_vials_to_full_run(vials: List[Dict[str, Any]]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, List[Tuple[int, int]]]:
    """
    Create full-run arrays by concatenating vials in order and provide:
      - vial_marker array (1..N)
      - windows list of (start,end) indices per vial
    """
    times: List[np.ndarray] = []
    masses: List[np.ndarray] = []
    rcons: List[np.ndarray] = []
    pcons: List[np.ndarray] = []
    vmarks: List[np.ndarray] = []
    windows: List[Tuple[int, int]] = []

    cursor = 0
    for idx, vial in enumerate(vials, start=1):
        t = np.asarray(vial["time"], dtype=float).squeeze()
        m = np.asarray(vial["mass"], dtype=float).squeeze()

        # conductivity may be None -> fill NaN to vial length
        r = vial.get("retentate_cond", None)
        p = vial.get("permeate_cond", None)

        def _align(arr, n):
            if arr is None:
                return np.full(n, np.nan, dtype=float)
            a = np.asarray(arr, dtype=float).squeeze()
            if a.size >= n:
                return a[:n]
            return np.pad(a, (0, n - a.size), constant_values=np.nan)

        r = _align(r, t.size)
        p = _align(p, t.size)

        times.append(t)
        masses.append(m)
        rcons.append(r)
        pcons.append(p)
        vmarks.append(np.full(t.shape, idx, dtype=float))

        if t.size > 0:
            start = cursor
            end = cursor + t.size - 1
            windows.append((start, end))
            cursor += t.size

    if not times:
        return (
            np.array([], dtype=float),
            np.array([], dtype=float),
            np.array([], dtype=float),
            np.array([], dtype=float),
            np.array([], dtype=float),
            [],
        )

    time = np.concatenate(times)
    mass = np.concatenate(masses)
    rcon = np.concatenate(rcons)
    pcon = np.concatenate(pcons)
    vmark = np.concatenate(vmarks)

    return time, mass, rcon, pcon, vmark, windows


# =============================================================================
# Window inference for flat arrays (mirrors vial marker strategy)
# =============================================================================

def _infer_windows_from_vial_marker(vial_marker: Optional[np.ndarray]) -> List[Tuple[int, int]]:
    if vial_marker is None:
        return []
    v_arr = np.asarray(vial_marker).squeeze()
    if v_arr.size == 0:
        return []

    # numeric if possible
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


def _infer_vial_windows_from_mat(
    format_name: str,
    time: np.ndarray,
    vial_marker: Optional[np.ndarray],
    data_stru_windows: Optional[List[Tuple[int, int]]] = None,
) -> List[Tuple[int, int]]:
    """
    Mirror the “same window logic” concept:

    - If data_stru: vial boundaries are the canonical windows.
    - Else: vial marker changes, else single window.
    """
    if format_name == "data_stru" and data_stru_windows:
        return data_stru_windows

    w = _infer_windows_from_vial_marker(vial_marker)
    if w:
        return w
    return _infer_windows_default(time)


# =============================================================================
# Per-vial slicing + mass cleaning (same as Excel script)
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
# Flat-array parsing
# =============================================================================

def _parse_flat_arrays(mat: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """
    Parse a generic MAT file with flat arrays.
    Requires Time and Mass. Conductivity and vial_marker optional.
    """
    time = _grab_case_insensitive(mat, "Time", "time", "t")
    mass = _grab_case_insensitive(mat, "Mass", "mass", "m")
    if time is None or mass is None:
        raise ValueError("MAT file missing required arrays: Time and Mass.")

    time = _as_1d_float_array(time)
    mass = _as_1d_float_array(mass)

    rcon = _grab_case_insensitive(mat, "RetentateCond", "Retentate_Cond", "retentate_cond", "ret_cond")
    pcon = _grab_case_insensitive(mat, "PermeateCond", "Permeate_Cond", "permeate_cond", "perm_cond")

    rcon_arr = _as_1d_float_array(rcon) if rcon is not None else np.full_like(time, np.nan, dtype=float)
    pcon_arr = _as_1d_float_array(pcon) if pcon is not None else np.full_like(time, np.nan, dtype=float)

    vial = _grab_case_insensitive(mat, "Vial", "VialSwap", "vial", "vial_swap", "vialswap")
    vial_arr = _as_1d_float_array(vial) if vial is not None else None

    # ensure same length as time
    n = time.size
    def _fit(x):
        if x is None:
            return None
        x = np.asarray(x).squeeze()
        if x.size >= n:
            return x[:n]
        return np.pad(x, (0, n - x.size), constant_values=np.nan)

    rcon_arr = _fit(rcon_arr)
    pcon_arr = _fit(pcon_arr)
    vial_arr = _fit(vial_arr)

    return time, mass, rcon_arr, pcon_arr, vial_arr


# =============================================================================
# Public API
# =============================================================================

def extract_mat_variables(path: str) -> Dict[str, Any]:
    """
    High-level functional entry point for MAT files.
    """
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

        conductivity_cF = getattr(ds, "conductivity_cF", None)
        if conductivity_cF is not None:
            try:
                conductivity_cF = bool(conductivity_cF)
            except Exception:
                conductivity_cF = None

        n_obj = getattr(ds, "n_obj", None)
        try:
            n_obj = int(n_obj) if n_obj is not None else None
        except Exception:
            n_obj = None

        data_config = _parse_data_config(ds)
        vials = _parse_data_raw_vials(ds)

        # Concatenate to full-run arrays + canonical windows + vial_marker
        time, mass, rcon, pcon, vial_marker, ds_windows = _concatenate_vials_to_full_run(vials)

        format_name = "data_stru"
        windows = _infer_vial_windows_from_mat(format_name, time, vial_marker, data_stru_windows=ds_windows)

    # ----------------------------------------------------------
    # Fallback: flat arrays
    # ----------------------------------------------------------
    else:
        dataset = None
        filename = None
        mode = None
        conductivity_cF = None
        n_obj = None
        data_config = {"extras": {}}

        time, mass, rcon, pcon, vial_marker = _parse_flat_arrays(mat)

        # build a single "vial" object for completeness
        vials = [{
            "number": 1,
            "time": time,
            "mass": mass,
            "retentate_cond": rcon,
            "permeate_cond": pcon,
            "nr": int(np.isfinite(mass).sum()),
            "extras": {},
        }]

        format_name = "flat_arrays"
        windows = _infer_vial_windows_from_mat(format_name, time, vial_marker, data_stru_windows=None)

    # Per-vial slices with cleaning (same as Excel script)
    per_t, per_m, per_r, per_p = _slice_per_vial_signals(time, mass, rcon, pcon, windows)

    variables: Dict[str, Any] = {
        # bookkeeping
        "path": path,
        "format": format_name,
        # data_stru core
        "dataset": dataset,
        "filename": filename,
        "mode": mode,
        "conductivity_cF": conductivity_cF,
        "n_obj": n_obj,
        # config + vials
        "data_config": data_config,
        "vials": vials,
        # full-run signals
        "time": time,
        "mass": mass,
        "retentate_cond": rcon,
        "permeate_cond": pcon,
        "vial_marker": vial_marker,
        # vial segmentation
        "windows": windows,
        "per_vial_time": per_t,
        "per_vial_mass": per_m,
        "per_vial_retentate_cond": per_r,
        "per_vial_permeate_cond": per_p,
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
