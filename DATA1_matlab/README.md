# DATA1-matlab

This folder contains MATLAB scripts and supporting functions for analyzing diafiltration and filtration experiments, performing parameter estimation, and generating heatmaps for Model-Based Design of Experiments (MBDoE).

## Dependencies

- **MATLAB**: Tested with versions R2020a and R2022a.

## Folder Organization

- **`data/`**  
  Contains experimental data files in `.mat` format used for simulations and analyses, also includes result files for paper figure preparation.

- **`functions/`**  
  Includes MATLAB functions for model configuration, data processing, and calculations.

## Key Files

### Simulation and Analysis Scripts

1. **`run_data_analysis.m`**  
    - Performs data analysis for filtration and diafiltration experiments.
    - Key features:
      - Configures datasets and model structures.
      - Fits experimental data to the model using a regression method.
      - Simulates and visualizes results for different datasets.
      - Outputs parameter estimates for permeability (`Lp`), solute permeability (`B`), and reflection coefficient (`sigma`).
      - Supports optional contour plotting for parameter exploration.
        
2. **`doe_heatmap_diafiltration.m`**  
    - Generates heatmaps for diafiltration experiments.
    - Analyzes key metrics like A-optimality, D-optimality, E-optimality, and modified E-optimality.
    - Simulates the effect of varying diafiltrate concentrations, applied pressures, and the number of end vials.
      
3. **`doe_heatmap_filtration.m`**  
    - Generates heatmaps for filtration experiments.
    - Focuses on the initial feed concentrations and applied pressures.
    - Evaluates similar optimality metrics to identify the best experimental conditions.

### Sigma Sensitivity Analysis

4. **`run_sigma_sensitivity.m`**  
   - Simulates sigma sensitivity for filtration and diafiltration experiments and generates comparative visualizations.
     
5. **`heatmap_sigma_sensitivity_diafiltration.m`**  
   - Generates heatmaps to study sigma sensitivity for diafiltration experiments.
   - Evaluates variations in predictions for mass, retentate, and permeate.

6. **`heatmap_sigma_sensitivity_filtration.m`**  
   - Similar to the diafiltration script but focuses on filtration experiments.

### Visualization Utilities

7. **`diafiltration_plots.py`**  
   - Python script stores functions generating enhanced plots and visualizations from simulation outputs.
   - Includes functions for simulation comparisons, contour plots, and sensitivity visualization.
   
7. **`DiafiltrationPaperPlots.ipynb`**
   - Figure reproduce for DATA1 paper.

9. **`PSE2021Plots.ipynb`**
   - Figure reproduce for DATA-MBDoE paper.
---

## Usage

1. **Setup**  
   - Ensure MATLAB is installed and include the `functions/` folder.
   - For Python-based visualization, set up a Python environment with the necessary libraries (`numpy`, `matplotlib`, `pandas`, `scipy`).

2. **Run Scripts**
   - For comprehensive data analysis: `run_data_analysis.m`.
   - For DOE heatmaps:
     - Diafiltration: `doe_heatmap_diafiltration.m`
     - Filtration: `doe_heatmap_filtration.m`
   - For sigma sensitivity:
     - Diafiltration: `heatmap_sigma_sensitivity_diafiltration.m`
     - Filtration: `heatmap_sigma_sensitivity_filtration.m`
     - General analysis: `run_sigma_sensitivity.m`

3. **Customize Parameters**  
   - Modify model parameters (e.g., model_stru.concpolar for inclusion of concentration polarization) directly in the scripts to explore different setups.
  
4. **Visualizations**  
   - Move results to `data/`
   - Refer to Jupyter notebooks to generate high-quality visualizations from simulation results.

## Outputs

- **Parameter Fits**  
  Parameter estimates and fitting results displayed in the console and saved to `.mat` files.

- **Simulation Results**  
  Processed simulation data and figure saved as `.csv` and `.png` files.

- **Heatmaps**  
  Visualizations saved as `.png` files in the `doe_FIM/` and `sigma_sensitivity/` folders.

---

For additional details, please refer to the inline comments in the MATLAB scripts.
