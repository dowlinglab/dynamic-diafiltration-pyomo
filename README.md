# Dynamic Diafiltration Model in Pyomo

**Characterizing Transport Properties of Nanofiltration Membranes**

This repository contains a Pyomo-based dynamic diafiltration model used to analyze the transport properties of nanofiltration (NF) membranes. It supports robust parameter estimation and data analytics workflows based on dynamic experiments.

---

## Environment Setup

This project is developed with the following dependencies:
- **Python**: Version 3.11.6
- **Pyomo**: Version 6.8.0

### Suggested Conda Configuration
Create an isolated environment with the following command:
```bash
conda create -n dynamic-diafiltration -c anaconda -c conda-forge -c IDAES-PSE python=3.11 numpy matplotlib pandas scipy idaes-pse
```

---

## Repository Organization

### Root Directory
1. **`utility.py`**  
   Library of functions for:
   - Loading and storing nested MATLAB `.mat` data.
   - Performing parameter estimation and providing information for model ranking.
   - Performing Fisher Information Matrix (FIM) analysis.
   - Visualization of experimental and simulation results.

2. **`DATA_model_demo.ipynb`**  
   A Jupyter Notebook demonstrating:
   - Dynamic diafiltration Pyomo modeling.
   - Analysis of NF90 and NF270 membranes.

3. **`data/`**  
   Contains formatted experimental data (nested `.mat` files) used for the model and analysis.

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
3. **Run the Notebook**: Execute `DATA_model_demo.ipynb` to:
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

*"Characterizing Transport Properties of Surface-Charged Nanofiltration Membranes via Model-based Data Analytics"*, authored by [Xinhong Liu et al.](mailto:wphillip@nd.edu), [William A. Phillip](mailto:wphillip@nd.edu), and [Alexander W. Dowling](mailto:adowling@nd.edu).

*"Membrane Characterization with Model-Based Design of Experiments"*, authored by [Xinhong Liu et al.](mailto:wphillip@nd.edu), [William A. Phillip](mailto:wphillip@nd.edu), and [Alexander W. Dowling](mailto:adowling@nd.edu).

*"DATA: Diafiltration Apparatus for high-Throughput Analysis"*, authored by [Jonathan A Ouimet](mailto:ouimetja@gmail.com), [Xinhong Liu et al.](mailto:wphillip@nd.edu), [William A. Phillip](mailto:wphillip@nd.edu), and [Alexander W. Dowling](mailto:adowling@nd.edu).


For details on the methodology, see the manuscript or contact the corresponding authors.
