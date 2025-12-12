"""Contract validator for the Python Experimental IR.

Use this to assert that loader outputs obey the invariants in SCHEMA.md.

Typical usage:
    from experiment_ir_validator import validate_experimental_ir
    validate_experimental_ir(vars_)
"""

from __future__ import annotations

from typing import Any, Dict, Iterable

import numpy as np


class IRValidationError(AssertionError):
    """Raised when the Experimental IR violates the schema contract."""


def _require_keys(d: Dict[str, Any], keys: Iterable[str], ctx: str) -> None:
    missing = [k for k in keys if k not in d]
    if missing:
        raise IRValidationError(f"Missing required keys in {ctx}: {missing}")


def validate_experimental_ir(vars_: Dict[str, Any]) -> None:
    """Validate the canonical Experimental IR.

    Raises:
        IRValidationError: if any invariant is violated.
    """
    _require_keys(
        vars_,
        [
            "path",
            "format",
            "dataset",
            "filename",
            "mode",
            "continuous_cF",
            "data_config",
            "vials",
            "time",
            "mass",
            "cV_avg",
            "cF_exp",
            "vial_ranges",
            "per_vial_time",
            "per_vial_mass",
            "per_vial_cV_avg",
            "per_vial_cF_exp",
        ],
        ctx="vars_",
    )

    if vars_["format"] not in {"data_stru", "flat_arrays"}:
        raise IRValidationError(f"format must be 'data_stru' or 'flat_arrays', got {vars_['format']!r}")

    cfg = vars_["data_config"]
    if not isinstance(cfg, dict):
        raise IRValidationError("data_config must be a dict")
    _require_keys(cfg, ["n", "nr"], ctx="data_config")

    vials = vars_["vials"]
    if not isinstance(vials, list):
        raise IRValidationError("vials must be a list")

    n = int(cfg["n"])
    if len(vials) != n:
        raise IRValidationError(f"len(vials)={len(vials)} does not match data_config['n']={n}")

    time = np.asarray(vars_["time"]).squeeze()
    mass = np.asarray(vars_["mass"]).squeeze()
    cV = np.asarray(vars_["cV_avg"]).squeeze()
    cF = np.asarray(vars_["cF_exp"]).squeeze()

    N = len(time)
    if not (len(mass) == len(cV) == len(cF) == N):
        raise IRValidationError(
            f"Full-run arrays must share length: len(time)={len(time)}, len(mass)={len(mass)}, len(cV_avg)={len(cV)}, len(cF_exp)={len(cF)}"
        )

    vial_ranges = vars_["vial_ranges"]
    if not isinstance(vial_ranges, list) or not all(isinstance(t, tuple) and len(t) == 2 for t in vial_ranges):
        raise IRValidationError("vial_ranges must be a list of (start,end) tuples")

    if len(vial_ranges) != n:
        raise IRValidationError(f"len(vial_ranges)={len(vial_ranges)} does not match n={n}")

    if N > 0:
        if vial_ranges[0][0] != 0:
            raise IRValidationError("vial_ranges must start at index 0")
        if vial_ranges[-1][1] != N - 1:
            raise IRValidationError("vial_ranges must end at index N-1")

    for i, (s, e) in enumerate(vial_ranges):
        if not (0 <= s <= e < N):
            raise IRValidationError(f"vial_ranges[{i}] out of bounds: ({s},{e}) for N={N}")
        if i < len(vial_ranges) - 1:
            s2, _ = vial_ranges[i + 1]
            if e + 1 != s2:
                raise IRValidationError(
                    f"vial_ranges must be contiguous: vial {i} ends at {e}, next vial starts at {s2}"
                )

    per_lists = [
        ("per_vial_time", vars_["per_vial_time"]),
        ("per_vial_mass", vars_["per_vial_mass"]),
        ("per_vial_cV_avg", vars_["per_vial_cV_avg"]),
        ("per_vial_cF_exp", vars_["per_vial_cF_exp"]),
    ]
    for name, lst in per_lists:
        if not isinstance(lst, list):
            raise IRValidationError(f"{name} must be a list")
        if len(lst) != n:
            raise IRValidationError(f"len({name})={len(lst)} does not match n={n}")

    for i, (s, e) in enumerate(vial_ranges):
        expected = e - s + 1
        for name, lst in per_lists:
            arr = np.asarray(lst[i]).squeeze()
            if len(arr) != expected:
                raise IRValidationError(
                    f"{name}[{i}] length {len(arr)} != expected {expected} from vial_ranges[{i}]=({s},{e})"
                )

    for i, vial in enumerate(vials):
        if not isinstance(vial, dict):
            raise IRValidationError(f"vials[{i}] must be a dict")
        _require_keys(vial, ["number", "time", "mass"], ctx=f"vials[{i}]")
