# DATA1-matlab

This folder contains MATLAB scripts and supporting functions for analyzing diafiltration and filtration experiments, performing parameter estimation, and generating heatmaps for Model-Based Design of Experiments (MBDoE).

## Dependencies

- **MATLAB**: Tested with versions R2020a and R2022a.

## Folder Organization

- **`data/`**  
  Contains experimental data files in `.mat` format used for simulations and analyses.

- **`functions/`**  
  Includes MATLAB functions for model configuration, data processing, and calculations.

## Key Files

### `run_data_analysis.m`
- Performs data analysis for filtration and diafiltration experiments.
- Key features:
  - Configures datasets and model structures.
  - Fits experimental data to the model using a regression method.
  - Simulates and visualizes results for different datasets.
  - Outputs parameter estimates for permeability (`Lp`), solute permeability (`B`), and reflection coefficient (`sigma`).
  - Supports optional contour plotting for parameter exploration.

### `doe_heatmap_diafiltration.m`
- Generates heatmaps for diafiltration experiments.
- Analyzes key metrics like A-optimality, D-optimality, E-optimality, and modified E-optimality.
- Simulates the effect of varying diafiltrate concentrations, applied pressures, and the number of end vials.

### `doe_heatmap_filtration.m`
- Generates heatmaps for filtration experiments.
- Focuses on the initial feed concentrations and applied pressures.
- Evaluates similar optimality metrics to identify the best experimental conditions.

## Usage

1. **Setup**  
   - Ensure MATLAB is installed and include the `functions/` folder.

2. **Run Scripts**
   - For comprehensive data analysis: Execute `run_data_analysis.m`.
   - For diafiltration DoE analysis: Execute `doe_heatmap_diafiltration.m`.
   - For filtration DoE analysis: Execute `doe_heatmap_filtration.m`.

3. **Customize Parameters**  
   - Modify model parameters (e.g., model_stru.concpolar for inclusion of concentration polarization) directly in the scripts to explore different setups.

## Outputs

- **Parameter Fits**  
  Parameter estimates and fitting results displayed in the console and saved to `.mat` files.

- **Simulation Results**  
  Processed simulation data saved as `.png` and `.csv` files.

- **Heatmaps**  
  Visualizations saved as `.png` files in the `doe_FIM/` folder, showcasing optimality metrics for experimental conditions.

---

For additional details, please refer to the inline comments in the MATLAB scripts.
