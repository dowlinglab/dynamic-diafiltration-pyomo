#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 17 15:03:12 2025
@author: kkasturi
"""
"""
Modified library of functions for diafiltration experiment modeling DATA1
by Xinhong Liu
"""

###############################################################################
# Libraries to import
###############################################################################

# Scientific computing and data handling
import numpy as np                  # For numerical operations and array handling
import pandas as pd                # For structured data manipulation and exporting (e.g., to CSV)
import scipy.io as spio            # For reading MATLAB .mat files
from scipy import interpolate      # For 1D interpolation of simulation results to match experiment times

# Plotting and visualization
import matplotlib.pyplot as plt    # For generating plots and visualizing results
from matplotlib.lines import Line2D # Optional: allows customizing legends or manual line drawing

# IDAES framework (used in some advanced Pyomo modeling, optional here)
import idaes                       # Provides access to advanced modeling tools in Pyomo (not directly used here)
import time                        # For tracking performance or timestamps (not yet used in this script)
import copy                        # For deep copying nested structures, used in sensitivity/FIM analysis
import os                          # For path and file operations, e.g., creating folders or saving files
import json                        # For exporting results as JSON files

# Evaluation and metrics
from sklearn.metrics import r2_score  # For computing R² between predicted and experimental values (optional/extra analysis)

# Core Pyomo modeling tools
from pyomo.environ import *       # Core modeling components: Var, Param, Constraint, Objective, etc.
from pyomo.dae import *           # Tools for differential algebraic equations (DAEs)
import idaes.core.util.scaling as iscale  # For variable scaling (optional, not currently used but useful in large models)

###############################################################################

###############################################################################

def loadmat(filename): # Defines the main function loadmat which takes one 
                       # argument: the path to a .mat file.

    '''
    Read in nested structure(mat file) generated from MATLAB and output dictionaries.

    This function is called instead of using spio.loadmat directly as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries which are still mat-objects.
    Adapted from https://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries

    Arguments:
        filename: str, location + filename

    Returns:
        dictionary contains structured data
    '''
    print("\nLoading data file =",filename,"\n") # Prints the filename to the 
                                                 # console for tracking/debugging

# ⬇ Helper Function: _check_keys(d)
    def _check_keys(d): # Defines a helper function to recursively clean a 
                        # dictionary loaded from MATLAB.
        '''
        checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        '''

        for key in d: # Loops over each key in the dictionary d.
            if isinstance(d[key], spio.matlab.mat_struct): # Checks if the value 
                                # is a mat_struct (used by scipy.io to represent 
                                # MATLAB structs).
                d[key] = _todict(d[key]) # If it is, convert it to a native Python dictionary using '_todict'

            elif isinstance(d[key], np.ndarray): # If it’s a NumPy array (possibly from a MATLAB cell array),
                d[key] = _tolist(d[key]) # ... convert it recursively into a Python list using _tolist.
        return d # Returns the cleaned dictionary.

# ⬇ Helper Function: _todict(matobj)
    def _todict(matobj): # Defines a recursive function to convert a MATLAB object
                         # into a nested Python dictionary.
        '''
        A recursive function which constructs from matobjects nested dictionaries
        '''
        d = {} # Creates an empty Python dictionary.

        for strg in matobj._fieldnames: # Loops through the field names (keys) of the MATLAB struct.
            elem = matobj.__dict__[strg] # Fetches the corresponding value for that key.

            if isinstance(elem, spio.matlab.mat_struct):
                d[strg] = _todict(elem) # If the value is another MATLAB struct, recurse into _todict.

            elif isinstance(elem, np.ndarray):
                d[strg] = _tolist(elem) # If the value is a NumPy array, process it via _tolist.

            else:
                d[strg] = elem # Otherwise, assign the value directly.

        return d # Returns the constructed nested dictionary.

# ⬇ Helper Function: _tolist(ndarray)
    def _tolist(ndarray): # Defines a recursive function to convert MATLAB cell arrays (NumPy ndarrays) into Python lists.
        '''
        A recursive function which constructs lists from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        '''
        elem_list = [] # Initialize an empty list.
        for sub_elem in ndarray: # Iterate through each element in the NumPy array.

            if isinstance(sub_elem, spio.matlab.mat_struct):
                elem_list.append(_todict(sub_elem)) # If it’s a MATLAB struct, convert to dict.

            elif isinstance(sub_elem, np.ndarray):
                elem_list.append(_tolist(sub_elem))# If it’s a nested array (e.g. cell arrays inside cell arrays), recurse.

            elif isinstance(sub_elem, int):
                elem_list.append(float(sub_elem)) # Convert integers to float for consistency (MATLAB doesn't distinguish well).

            else:
                elem_list.append(sub_elem) # Otherwise, append the raw value.

        return elem_list # Returns the final list.

# ⬇ Load .mat file and convert
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True) # Loads the .mat file using scipy.io.loadmat with two important options:
                                            # struct_as_record=False: Load MATLAB structs as Python classes with attributes
                                            # squeeze_me=True: Remove singleton dimensions (e.g. [1x1] arrays → scalar)
    return _check_keys(data) # Final step: clean up all mat_struct and array remnants, then return the nested dictionary.