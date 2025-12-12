# Experimental Data IR Schema (MATLAB ↔ Python)

This document defines the **canonical Python Experimental IR** produced by `extract_mat_variables(path)` and its **mapping** from MATLAB `data_stru` as produced by `load_data.m`.


## 1. Canonical Python Experimental IR

`extract_mat_variables(path) -> dict[str, Any]`

### 1.1 Top-level keys

#### Bookkeeping (required)
- `path: str`  
  Path used for loading.
- `format: "data_stru" | "flat_arrays"`  
  Which parsing branch was used.

#### Experiment metadata (required keys; values may be `None`)
- `dataset: float | None`
- `filename: str | None`
- `mode: str | None`
- `continuous_cF: bool | None`  
  **Important:** `None` means the `.mat` file did not store this flag (no inference under Suggestion A).

#### Configuration (required)
- `data_config: dict[str, Any]`  
  Contains configuration and parameter initialization data.
  Must include:
  - `n: int` total number of vials (computed if missing)
  - `nr: int` total number of residuals (computed if missing)
  May include:
  - `extras: dict[str, Any]` for non-canonical / unknown config fields.

#### Vials (required)
- `vils: list[dict[str, Any]]` length `n`

Each vial dict uses **MATLAB-native field names**:

Required fields per vial:
- `number: int`
- `time: np.ndarray` shape `(Ti,)`
- `mass: np.ndarray` shape `(Ti,)`

Optional per-vial measurement fields:
- `cV_avg: np.ndarray | None` shape `(Ti,)`
- `cF_exp: np.ndarray | None` shape `(Ti,)`
- `nr: int | Any`

Optional extras:
- `extras: dict[str, Any]` preserves additional vial fields.

#### Full-run concatenated arrays (required)
These are derived by concatenating vial arrays in order:

- `time: np.ndarray` shape `(N,)`
- `mass: np.ndarray` shape `(N,)`
- `cV_avg: np.ndarray` shape `(N,)` (NaN-filled if missing)
- `cF_exp: np.ndarray` shape `(N,)` (NaN-filled if missing)

Optional:
- `vial_marker: np.ndarray | None` shape `(N,)` indicating vial membership.

#### Vial segmentation index ranges (required)
- `vial_ranges: list[tuple[int, int]]` length `n`  
  Each tuple is inclusive `(start_idx, end_idx)` and partitions `0..N-1`.

#### Per-vial derived slices (required)
- `per_vial_time: list[np.ndarray]`
- `per_vial_mass: list[np.ndarray]`  
  **Semantics:** re-baselined per vial and spike-filtered per the loader policy.
- `per_vial_cV_avg: list[np.ndarray]`
- `per_vial_cF_exp: list[np.ndarray]`

---

## 2. Global invariants (must hold)

Let `N = len(time)`.

1) Length alignment:
- `len(time) == len(mass) == len(cV_avg) == len(cF_exp) == N`

2) Vial list alignment:
- `len(vials) == data_config["n"] == len(vial_ranges) == len(per_vial_time) == len(per_vial_mass) == len(per_vial_cV_avg) == len(per_vial_cF_exp)`

3) Vial range coverage and contiguity:
- `vial_ranges[0][0] == 0`
- `vial_ranges[-1][1] == N - 1`
- For each `i`:
  - `0 <= start_i <= end_i < N`
- For each adjacent pair:
  - `end_i + 1 == start_{i+1}`

4) Per-vial slice consistency:
For each vial `i` with range `(s, e)`:
- `per_vial_time[i].shape[0] == e - s + 1`  
(and same for mass/cV_avg/cF_exp slices)

---

## 3. MATLAB ↔ Python mapping

### 3.1 MATLAB source schema

MATLAB struct:
- `data_stru.dataset`
- `data_stru.filename`
- `data_stru.mode`
- `data_stru.continuous_cF` (may be absent in stored `.mat` files)
- `data_stru.data_config` struct
- `data_stru.data_raw(i)` array of structs, one per vial:
  - `.number`, `.time`, `.mass`, `.cV_avg`, `.cF_exp`, `.nr`, and possibly additional fields.

### 3.2 Field mapping rules

Metadata:
- `data_stru.dataset`  → `dataset`
- `data_stru.filename` → `filename`
- `data_stru.mode`     → `mode`
- `data_stru.continuous_cF` → `continuous_cF`  
  If absent in `.mat`, Python returns `None` (no inference).

Config:
- Each `data_stru.data_config.<field>`:
  - canonical fields → `data_config["<field>"]`
  - unknown fields   → `data_config["extras"]["<field>"]`
- If `n` or `nr` are missing, Python computes them.

Vials:
For each MATLAB vial `data_raw(i)`:
- `.number` → `vials[i]["number"]`
- `.time`   → `vials[i]["time"]`
- `.mass`   → `vials[i]["mass"]`
- `.cV_avg` → `vials[i]["cV_avg"]` (may be None/NaN-filled later)
- `.cF_exp` → `vials[i]["cF_exp"]` (may be None/NaN-filled later)
- `.nr`     → `vials[i]["nr"]`
- other fields → `vials[i]["extras"][field]`

Derived objects:
- Full-run arrays are concatenations of vial arrays.
- `vial_ranges` are computed from vial lengths (canonical for `data_stru` inputs).

---

## 4. Design note: Suggestion A (no inference for continuous_cF)

If the `.mat` file does not explicitly store `continuous_cF` (or a recognized equivalent), the loader returns:

- `continuous_cF = None`

Downstream code should treat this as “unknown / not recorded”.

---

## 5. References (engineering rationale)

- Wilson et al., *Best Practices for Scientific Computing*, PLoS Biology (2014).
- Sculley et al., *Hidden Technical Debt in Machine Learning Systems*, NeurIPS (2015).
