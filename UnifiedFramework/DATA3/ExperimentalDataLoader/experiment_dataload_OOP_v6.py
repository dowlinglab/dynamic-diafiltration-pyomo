#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unified experimental data loader for diafiltration projects.
Supports both:
    (1) MATLAB .mat files containing data_stru (load_data.m output)
    (2) Excel .xlsx files where each sheet contains an experiment

This version:
    - Adds PyDataStru (Python representation of MATLAB data_stru)
    - Adds robust constructors from_mat() and from_excel()
    - Keeps existing DataLoader fully intact
"""

from __future__ import annotations

# ================================================================
# Imports
# ================================================================
import os
import numpy as np
import pandas as pd
import scipy.io as sio
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict, Any, Protocol, runtime_checkable

# ================================================================
# SECTION 1 — PyDataStru (NEW)
# ================================================================

@dataclass
class ExperimentMetadata:
    """ This contains metadata about the experiment.
    It mimics the MATLAB data_stru.config. 
    """
    M_F0: float | None = None
    C_F0: np.ndarray | None = None
    M_O: float | None = None
    C_D: np.ndarray | None = None
    nc: int | None = None
    namec: Any | None = None
    ni: int | None = None
    delP: float | None = None
    Temp: float | None = None
    Am: float | None = None
    rho: float | None = None
    Lp0: float | None = None
    B0: np.ndarray | None = None
    sigma0: float | None = None
    theta0: np.ndarray | None = None
    n: int | None = None
    nr: int | None = None
    n_extra: int | None = None
    extras: Dict[str, Any] = field(default_factory=dict)

@dataclass
class VialData:
    """ This contains experimental data for a single vial.
    It mimics the MATLAB data_stru.data_raw(i).
    """
    number: int
    time: np.ndarray
    mass: np.ndarray
    cV_avg: np.ndarray | None = None
    cF_exp: np.ndarray | None = None
    nr: int | None = None
    extras: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ExperimentData:
    """
    This contains all of the data for an experiment.
    It mimics the MATLAB data_stru."""
    dataset: float | None
    filename: str | None
    mode: str | None
    conductivity_cF: bool | None
    n_obj: int | None
    data_config: ExperimentMetadata
    data_raw: List[VialData]

    # ------------------------------------------------------------
    # Constructor from .mat (load_data.m output)
    # ------------------------------------------------------------
    # Q: Why is this a static method?
    # A: Work around for an error?
    @staticmethod
    def from_mat(path: str):
        mat = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
        if "data_stru" not in mat:
            raise ValueError(f"{path} does not contain data_stru")

        ds = mat["data_stru"]

        dataset = float(getattr(ds, "dataset", np.nan))
        filename = str(getattr(ds, "filename", None))
        mode = str(getattr(ds, "mode", None))
        conductivity_cF = bool(getattr(ds, "conductivity_cF", False))
        n_obj = int(getattr(ds, "n_obj", 0))

        # ---- data_config ----
        dc = getattr(ds, "data_config", None)
        cfg_kwargs = {}
        extras = {}

        if dc is not None and hasattr(dc, "_fieldnames"):
            for fn in dc._fieldnames:
                val = getattr(dc, fn)
                if isinstance(val, (list, np.ndarray)):
                    val = np.asarray(val).squeeze()
                if fn in .__dataclass_fields__:
                    cfg_kwargs[fn] = val
                else:
                    extras[fn] = val

        data_config = (**cfg_kwargs, extras=extras)

        # ---- data_raw (vials) ----
        dr = getattr(ds, "data_raw", [])
        if isinstance(dr, np.ndarray):
            dr = dr.ravel()

        vials = []
        for i, v in enumerate(dr):
            t = np.asarray(getattr(v, "time", []), dtype=float).squeeze()
            m = np.asarray(getattr(v, "mass", []), dtype=float).squeeze()

            cV = getattr(v, "cV_avg", None)
            if cV is not None:
                cV = np.asarray(cV).squeeze()

            cF = getattr(v, "cF_exp", None)
            if cF is not None:
                cF = np.asarray(cF).squeeze()

            nr = getattr(v, "nr", None)
            num = getattr(v, "number", i + 1)

            vial_extras = {}
            if hasattr(v, "_fieldnames"):
                for fn in v._fieldnames:
                    if fn not in ["time", "mass", "cV_avg", "cF_exp", "nr", "number"]:
                        vial_extras[fn] = getattr(v, fn)

            vials.append(VialData(num, t, m, cV, cF, nr, vial_extras))

        return ExperimentData(dataset, filename, mode, conductivity_cF, n_obj, data_config, vials)

    # ------------------------------------------------------------
    # Constructor from Excel (robust header detection)
    # ------------------------------------------------------------
    @staticmethod
    def from_excel(path: str, sheet: str):
        raw = pd.read_excel(path, sheet_name=sheet, header=None)

        # ---- Header row detection (same logic as ExcelSource) ----
        def detect_header_row(df, max_rows=20):
            for i in range(min(max_rows, len(df))):
                row = df.iloc[i].astype(str).str.lower().tolist()
                has_time = any("time" in s for s in row)
                has_mass_or_cond = any(("mass" in s) or ("cond" in s) or ("pressure" in s) for s in row)
                if has_time and has_mass_or_cond:
                    return i
            return 0

        hdr = detect_header_row(raw)

        # Build measurement table
        df = raw.iloc[hdr + 1:].copy()
        df.columns = raw.iloc[hdr].astype(str).str.strip().tolist()
        df = df.reset_index(drop=True)

        # ---- Locate time and mass columns ----
        def find_col(df, keywords):
            for col in df.columns:
                low = str(col).lower()
                if any(k in low for k in keywords):
                    return col
            return None

        tcol = find_col(df, ["time", "sec", "seconds", "record"])
        mcol = find_col(df, ["mass", "weight", "collection", "permeate"])

        if tcol is None:
            raise ValueError(f"No time column found. Columns = {list(df.columns)}")
        if mcol is None:
            raise ValueError(f"No mass column found. Columns = {list(df.columns)}")

        # Convert
        time = pd.to_numeric(df[tcol], errors="coerce").to_numpy()
        mass = pd.to_numeric(df[mcol], errors="coerce").to_numpy()

        cF = df.filter(regex="retent.*cond", case=False)
        cV = df.filter(regex="permeat.*cond", case=False)

        vial = VialData(
            number=1,
            time=time,
            mass=mass,
            cF_exp=cF.to_numpy().squeeze() if not cF.empty else None,
            cV_avg=cV.to_numpy().squeeze() if not cV.empty else None,
            nr=np.isfinite(mass).sum(),
            extras={}
        )

        return ExperimentData(
            dataset=None,
            filename=os.path.basename(path),
            mode=None,
            conductivity_cF=False,
            n_obj=1,
            data_config=(),
            data_raw=[vial]
        )


# ================================================================
# SECTION 2 — DataRun (unchanged)
# ================================================================

# How does this relate to our new plan with two dataclasses (ExperimentData and VialData)?

# This is our first pass at unifying the data structure from >1 month ago. It only grabs minimal data.
# We will likely want to update the code below to use the more robust ExperimentData and VialData classes above.
@dataclass
class DataRun:
    time: np.ndarray
    mass: np.ndarray
    retentate_cond: np.ndarray
    permeate_cond: np.ndarray
    vial_marker: Optional[np.ndarray] = None
    extra_signals: Dict[str, np.ndarray] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        n = len(self.time)

        # This looks like a helper function that can be moved
        # to the global scope in case we want to reuse it elsewhere.
        def pad(x):
            if x is None:
                return np.full(n, np.nan)
            arr = np.asarray(x).squeeze()
            if len(arr) < n:
                arr = np.pad(arr, (0, n - len(arr)), constant_values=np.nan)
            return arr[:n]

        self.time = pad(self.time)
        self.mass = pad(self.mass)
        self.retentate_cond = pad(self.retentate_cond)
        self.permeate_cond = pad(self.permeate_cond)

        if self.vial_marker is not None:
            self.vial_marker = pad(self.vial_marker)

        # Normalize extra signals
        for k, v in list(self.extra_signals.items()):
            self.extra_signals[k] = pad(v)


# ================================================================
# SECTION 3 — ExcelSource (unchanged)
# ================================================================

@runtime_checkable
class DataSource(Protocol):
    def sniff(self) -> bool: ...
    def parse(self) -> DataRun: ...


# This gets simplified into a function
class ExcelSource:
    """Original Excel parser used by DataLoader."""
    def __init__(self, path: str, measure_sheet: Optional[str] = None, index_sheet: Optional[str] = None):
        self.path = path
        self.measure_sheet_override = measure_sheet
        self.index_sheet_override = index_sheet
        self.measure_sheet = None
        self.index_sheet = None

    def sniff(self):
        return os.path.splitext(self.path)[1].lower() in [".xlsx", ".xls", ".xlsm"]

    # --- same detection logic used earlier ---
    @staticmethod
    def detect_header_row(df):
        for i in range(min(20, len(df))):
            row = df.iloc[i].astype(str).str.lower().tolist()
            if "time" in " ".join(row) and ("mass" in " ".join(row) or "cond" in " ".join(row)):
                return i
        return 0

    def parse(self) -> DataRun:
        xl = pd.ExcelFile(self.path)
        sheet = self.measure_sheet_override or xl.sheet_names[0]

        raw = xl.parse(sheet, header=None)
        hdr = self.detect_header_row(raw)

        df = raw.iloc[hdr + 1:].copy()
        df.columns = raw.iloc[hdr].astype(str).str.strip().tolist()
        df = df.reset_index(drop=True)

        # Extract canonical fields
        tcol = [c for c in df.columns if "time" in str(c).lower()][0]
        mcol = [c for c in df.columns if "mass" in str(c).lower()][0]

        time = pd.to_numeric(df[tcol], errors="coerce").to_numpy()
        mass = pd.to_numeric(df[mcol], errors="coerce").to_numpy()

        r = df.filter(regex="retent.*cond", case=False)
        p = df.filter(regex="permeat.*cond", case=False)

        return DataRun(
            time=time,
            mass=mass,
            retentate_cond=r.to_numpy().squeeze() if not r.empty else None,
            permeate_cond=p.to_numpy().squeeze() if not p.empty else None,
            vial_marker=None,
        )


# ================================================================
# SECTION 4 — MatSource (unchanged)
# ================================================================

class MatSource:
    """Original MATLAB parser for DataLoader."""
    def __init__(self, path: str):
        self.path = path

    def sniff(self):
        return self.path.lower().endswith(".mat")

    def parse(self) -> DataRun:
        mat = sio.loadmat(self.path, squeeze_me=True, struct_as_record=False)

        if "data_stru" in mat:
            ds = mat["data_stru"]

            # Concatenate all vials into one DataRun
            times, masses, rcons, pcons, vm = [], [], [], [], []

            raw = ds.data_raw
            if isinstance(raw, np.ndarray):
                raw = raw.ravel()

            for vial in raw:
                t = np.asarray(getattr(vial, "time", []), dtype=float)
                m = np.asarray(getattr(vial, "mass", []), dtype=float)
                r = np.asarray(getattr(vial, "cF_exp", []), dtype=float)
                p = np.asarray(getattr(vial, "cV_avg", []), dtype=float)
                num = getattr(vial, "number", 1)

                times.append(t)
                masses.append(m)
                rcons.append(r)
                pcons.append(p)
                vm.append(np.full(len(t), num))

            return DataRun(
                time=np.concatenate(times),
                mass=np.concatenate(masses),
                retentate_cond=np.concatenate(rcons),
                permeate_cond=np.concatenate(pcons),
                vial_marker=np.concatenate(vm),
            )

        # If no data_stru, fallback
        raise ValueError("MAT file does not contain data_stru")


# ================================================================
# SECTION 5 — DataLoader (unchanged)
# ================================================================

class DataLoader:
    """Uses ExcelSource or MatSource to build DataRun, then slices into windows."""
    def __init__(self, path: str, measure_sheet: Optional[str] = None):
        self.path = path
        self.measure_sheet = measure_sheet
        self.source = None
        self.windows = []
        self.times = []
        self.masses = []
        self.retentate_cond = []
        self.permeate_cond = []

    def load(self):
        excel = ExcelSource(self.path, self.measure_sheet)
        mat = MatSource(self.path)

        if excel.sniff():
            self.source = excel
        elif mat.sniff():
            self.source = mat
        else:
            raise ValueError("Unknown file type")

        run = self.source.parse()

        # One window for now
        n = len(run.time)
        self.windows = [(0, n - 1)]

        self.times = [run.time]
        self.masses = [run.mass]
        self.retentate_cond = [run.retentate_cond]
        self.permeate_cond = [run.permeate_cond]


# ================================================================
# SECTION 6 — CLI (unchanged)
# ================================================================
def list_sheets(path: str):
    if not path.lower().endswith(".xlsx"):
        return []
    return pd.ExcelFile(path).sheet_names


def _main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("path")
    parser.add_argument("--quicklook", action="store_true")
    args = parser.parse_args()

    dl = DataLoader(args.path)
    dl.load()
    print("Windows:", dl.windows)


if __name__ == "__main__":
    _main()
