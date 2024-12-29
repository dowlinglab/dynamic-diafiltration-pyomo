# Dynamic Diafiltration Model in Pyomo

**Characterizing Transport Properties of Nanofiltration Membranes**

This repository contains a Pyomo-based dynamic diafiltration model used to analyze the transport properties of nanofiltration (NF) membranes. It supports robust parameter estimation and data analytics workflows based on dynamic experiments.

---

## Environment Setup

This project is developed with the following dependencies:
- **Python**: Version 3.11.6
- **Pyomo**: Version 6.8.0
- **MATLAB**: Implemented with R2020a and R2022a (for additional MATLAB analysis in `DATA1-matlab`)

### Suggested Conda Configuration

To set up the Python environment, recommend the following Conda configuration:
```bash
conda create -n dynamic-diafiltration -c anaconda -c conda-forge -c IDAES-PSE python=3.11 numpy matplotlib pandas scipy idaes-pse scikit-learn
```
Activate the environment:
```bash
conda activate dynamic-diafiltration
```
In case CasADi is not installed (needed for model initialization)
```bash
pip install casadi
```
---

## Repository Organization

### Root Directory

1. **`data_library/`**  
   Contains formatted experimental data (nested `.mat` files) used for the model and analysis.

2. **`DATA1-matlab/`**
   Original MATLAB scripts and supporting functions for analyzing diafiltration and filtration experiments, performing parameter estimation, and generating heatmaps for Model-Based Design of Experiments (MBDoE) in DATA1 and DATA-MBDoE papers. For usage details, please refer to the README file in this folder.
   
3. **`utility.py`**  
   Core library of Python functions for:
   - Loading and storing nested MATLAB `.mat` data.
   - Performing parameter estimation and providing information for model ranking.
   - Performing Fisher Information Matrix (FIM) analysis.
   - Visualization of experimental and simulation results.

4. **`DATA1_model_demo.ipynb`**  
   A Jupyter Notebook demonstrating:
   - Dynamic diafiltration Pyomo modeling.
   - Reproduce DATA1 analysis for NF90 membranes.
   - [WIP] Generating contours for sensitivity analysis/DoE in the DATA1 paper.
  
5. **`DATA2_model_demo.ipynb`**  
   A Jupyter Notebook demonstrating:
   - Dynamic diafiltration Pyomo modeling.
   - Analysis of NF270 membranes for DATA2 paper.

6. **`DATA2_visualization.ipynb`**  
   Notebook for generating publication-ready plots and visualizations for the DATA2 paper.

7. **`run_cross_verification.py`**  
   Python script for cross-verifying empirical solute permeability coefficient (`B`) models with additional datasets, and store visualizations.

8. **`run_DATA2_model_variations.py`**  
   Script for exploring model variations (e.g., including startup dynamics and/or time correction) and performing Fisher Information Matrix (FIM) calculations to evaluate information gain for the DATA2 paper.

9. **`run_pre_B_dependence.py`**  
   Preliminarily investigates the dependence of solute permeability coefficient (`B`) on interface concentrations. Outputs visualizations and insights for modeling and optimization.
   
---

## Key Features

- **Dynamic Diafiltration Modeling**: Captures startup dynamics, concentration-dependent transport properties, and solute flux mechanisms.
- **Data Analytics**: Implements advanced techniques like weighted least squares (WLS) and Akaike Information Criterion (AIC) for parameter estimation and model discrimination.
- **Support for Multiple Membranes**: Includes NF90 (near neutral) and NF270 (negative-charged) membrane analysis for comparing transport properties.
- **Visualization Tools**: Built-in scripts for generating model diagnostics and experiment comparisons.
- **FIM Analysis**: Calculates the Fisher Information Matrix for experimental design optimization.

---

## Usage

### Running the Model
1. **Prepare Data**: Place experimental `.mat` files in the `data/` directory.
2. **Configure Parameters**: Modify initial conditions and parameters in `utility.py` or the notebook as needed.
3. **Run the Notebook**: Execute `DATA2_model_demo.ipynb` to:
   - Load data.
   - Run the Pyomo model.
   - Visualize results.

### Visualization
Plots include:
- Time evolution of mass and concentration in vials.
- Comparison of experimental and simulated results.

### FIM Analysis
Use `calc_FIM` to optimize experimental design:
- Configure step size and finite difference scheme.
- Analyze parameter sensitivity printout.

---

## References

This repository accompanies the studies: 

**DATA2**

*"Characterizing Transport Properties of Surface-Charged Nanofiltration Membranes via Model-based Data Analytics"*, authored by [Xinhong Liu et al.](mailto:xliu27@alumni.nd.edu), [William A. Phillip](mailto:wphillip@nd.edu), and [Alexander W. Dowling](mailto:adowling@nd.edu).

**DATA-MBDoE**

*"Membrane Characterization with Model-Based Design of Experiments"*, authored by [Xinhong Liu et al.](mailto:xliu27@alumni.nd.edu), [William A. Phillip](mailto:wphillip@nd.edu), and [Alexander W. Dowling](mailto:adowling@nd.edu).

**DATA1**

*"DATA: Diafiltration Apparatus for high-Throughput Analysis"*, authored by [Jonathan A Ouimet](mailto:ouimetja@gmail.com), [Xinhong Liu et al.](mailto:xliu27@alumni.nd.edu), [William A. Phillip](mailto:wphillip@nd.edu), and [Alexander W. Dowling](mailto:adowling@nd.edu).


For details on the methodology, refer to the manuscript or contact the corresponding authors.
