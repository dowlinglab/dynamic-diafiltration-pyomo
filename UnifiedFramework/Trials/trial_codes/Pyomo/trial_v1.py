#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 21:09:57 2025

@author: kkasturi
"""

# -*- coding: utf-8 -*-
"""
Reactor DoE with Pyomo.DAE + pyomo.contrib.doe
A --> B --> C with Arrhenius kinetics; temperature is piecewise-constant.
"""
import json
from pathlib import Path

import numpy as np
import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.contrib.doe import DesignOfExperiments


class ReactorExperiment:
    """
    Minimal experiment wrapper providing:
      - create_model()
      - finalize_model()
      - label_experiment()
      - get_labeled_model()
    Compatible with pyomo.contrib.doe.DesignOfExperiments.
    """

    def __init__(self, data, nfe=10, ncp=3):
        self.data = data
        self.nfe = nfe
        self.ncp = ncp
        self.model = None
        self.t_control = None

    # === Build the continuous-time model ===
    def create_model(self):
        m = self.model = pyo.ConcreteModel()

        # Parameters
        m.R = pyo.Param(mutable=False, initialize=8.314)  # J/mol-K

        # Time
        m.t = ContinuousSet(bounds=(0.0, 1.0))

        # States
        m.CA = pyo.Var(m.t, within=pyo.NonNegativeReals)
        m.CB = pyo.Var(m.t, within=pyo.NonNegativeReals)
        m.CC = pyo.Var(m.t, within=pyo.NonNegativeReals)
        m.T = pyo.Var(m.t, within=pyo.NonNegativeReals)

        # Kinetic parameters (unknowns to estimate)
        m.A1 = pyo.Var(within=pyo.NonNegativeReals)
        m.E1 = pyo.Var(within=pyo.NonNegativeReals)
        m.A2 = pyo.Var(within=pyo.NonNegativeReals)
        m.E2 = pyo.Var(within=pyo.NonNegativeReals)

        # Derivatives
        m.dCAdt = DerivativeVar(m.CA, wrt=m.t)
        m.dCBdt = DerivativeVar(m.CB, wrt=m.t)

        # Rate constants (Arrhenius)
        def k1_rule(m, t):
            return m.A1 * pyo.exp(-m.E1 * 1000.0 / (m.R * m.T[t]))

        def k2_rule(m, t):
            return m.A2 * pyo.exp(-m.E2 * 1000.0 / (m.R * m.T[t]))

        m.k1 = pyo.Expression(m.t, rule=k1_rule)
        m.k2 = pyo.Expression(m.t, rule=k2_rule)

        # ODEs
        def ca_ode(m, t):
            return m.dCAdt[t] == -m.k1[t] * m.CA[t]

        def cb_ode(m, t):
            return m.dCBdt[t] == m.k1[t] * m.CA[t] - m.k2[t] * m.CB[t]

        m.CA_rxn_ode = pyo.Constraint(m.t, rule=ca_ode)
        m.CB_rxn_ode = pyo.Constraint(m.t, rule=cb_ode)

        # Equimolar mass balance: CA0 = CA + CB + CC
        def cc_balance(m, t):
            return m.CA[m.t.first()] == m.CA[t] + m.CB[t] + m.CC[t]

        m.CC_balance = pyo.Constraint(m.t, rule=cc_balance)

        return m

    # === Discretize + impose experiment structure ===
    def finalize_model(self):
        m = self.model
        data = self.data

        # Control points (dict: time -> temperature)
        control_points = data["control_points"]  # {0.0: 300, 0.125: 300, ...}

        # Initial conditions and bounds
        m.CA[m.t.first()].value = data["CA0"]
        m.CA[m.t.first()].setlb(data["CA_bounds"][0])
        m.CA[m.t.first()].setub(data["CA_bounds"][1])

        m.CB[m.t.first()].fix(data["CB0"])  # allows nonzero CB0 if desired

        # Update time set with experiment horizon + control times
        m.t.update(data["t_range"])
        m.t.update(control_points.keys())

        # Fix unknown parameters to the provided "true" (or nominal) values for FIM
        m.A1.fix(data["A1"])
        m.A2.fix(data["A2"])
        m.E1.fix(data["E1"])
        m.E2.fix(data["E2"])

        # Discretize
        discr = pyo.TransformationFactory("dae.collocation")
        discr.apply_to(m, nfe=self.nfe, ncp=self.ncp, wrt=m.t)

        # Temperature bounds and piecewise-constant control:
        T_lo, T_hi = data["T_bounds"]
        for t in m.t:
            m.T[t].setlb(T_lo)
            m.T[t].setub(T_hi)

        # Fix T at control points to specified values
        for tc, val in control_points.items():
            m.T[tc].fix(val)

        # Between control points, hold T equal to the most recent control value
        control_times = sorted(control_points.keys())
        m.T_control = pyo.ConstraintList()
        for t in m.t:
            if t in control_points:
                continue
            # find neighbor control time to the left
            left = max(ct for ct in control_times if ct <= t)
            m.T_control.add(m.T[t] == m.T[left])

        # Store for labeling
        self.t_control = control_times

    # === Tag inputs/outputs/errors/unknowns for DoE ===
    def label_experiment(self):
        m = self.model

        # Measurement outputs at control times
        m.experiment_outputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        for t in self.t_control:
            m.experiment_outputs[m.CA[t]] = None
            m.experiment_outputs[m.CB[t]] = None
            m.experiment_outputs[m.CC[t]] = None

        # Measurement errors (diagonal, constant)
        m.measurement_error = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        concentration_error = 1e-2
        for t in self.t_control:
            m.measurement_error[m.CA[t]] = concentration_error
            m.measurement_error[m.CB[t]] = concentration_error
            m.measurement_error[m.CC[t]] = concentration_error

        # Design variables (inputs): initial CA and T at control points
        m.experiment_inputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        m.experiment_inputs[m.CA[m.t.first()]] = None
        for t in self.t_control:
            m.experiment_inputs[m.T[t]] = None

        # Unknown parameters to estimate
        m.unknown_parameters = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        for k in [m.A1, m.A2, m.E1, m.E2]:
            m.unknown_parameters[k] = pyo.value(k)

    def get_labeled_model(self):
        if self.model is None:
            self.create_model()
            self.finalize_model()
            self.label_experiment()
        return self.model


if __name__ == "__main__":
    # === Load experiment data ===
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "result.json"
    with open(file_path) as f:
        data_ex = json.load(f)

    # Ensure control points are float keys
    data_ex["control_points"] = {float(k): v for k, v in data_ex["control_points"].items()}

    # Build experiment
    experiment = ReactorExperiment(data=data_ex, nfe=10, ncp=3)

    # === Design of Experiments ===
    # Finite-difference sensitivity options
    fd_formula = "central"
    step_size = 1e-3

    # Objective: log-det(FIM) (a.k.a. D-optimal); uses Cholesky with lower bound
    doe_obj = DesignOfExperiments(
        experiment,
        fd_formula=fd_formula,
        step=step_size,
        objective_option="determinant",
        scale_constant_value=1.0,
        scale_nominal_param_value=True,
        prior_FIM=None,
        jac_initial=None,
        fim_initial=None,
        L_diagonal_lower_bound=1e-7,
        solver=None,
        tee=False,
        get_labeled_model_args=None,
        _Cholesky_option=True,
        _only_compute_fim_lower=True,
    )

    # Full-factorial ranges:
    #   "min, max, npoints" per variable
    design_ranges = {
        "CA[0.0]": [1.0, 9.0, 3],      # initial CA at t=0.0
        "T[0.0]": [300.0, 700.0, 3],   # initial temperature at first control point
    }

    # Compute FIM over factorial design (sequential simulation for speed)
    doe_obj.compute_FIM_full_factorial(design_ranges=design_ranges, method="sequential")

    # Plot (saved as example_reactor_compute_FIM.png/pdf)
    doe_obj.draw_factorial_figure(
        sensitivity_design_variables=["CA[0.0]", "T[0.0]"],
        fixed_design_variables={
            # hold all other control temperatures at 300 K
            "T[0.125]": 300,
            "T[0.25]": 300,
            "T[0.375]": 300,
            "T[0.5]": 300,
            "T[0.625]": 300,
            "T[0.75]": 300,
            "T[0.875]": 300,
            "T[1.0]": 300,
        },
        title_text="Reactor Example",
        xlabel_text="Concentration of A (M)",
        ylabel_text="Initial Temperature (K)",
        figure_file_name="example_reactor_compute_FIM",
        log_scale=False,
    )

    # Optional: run a single DoE optimization (depends on your objective settings)
    doe_obj.run_doe()
