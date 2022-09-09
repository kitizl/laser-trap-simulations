import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess

# =============================================================================
# PARAMETRIC RESONANCE
# Two main adjustment parameters:
#	1. Modulation depth
#	2. Phase-shift
# Plot the effects they have on the amplitude (as a ratio of the control
# amplitude) to see what the ideal parameters, mathematically, are.

# 100 values linear spaced between 0.01 and 1.0
mod_depth_array = np.linspace(0.2, 0.3, 100)

# run the control experiment
print("Running experiments")


proc = subprocess.Popen(['wsl', './control', "data/control.csv", 
						"71.5e-9", "5.0", "150e3","1"])
proc.wait()


# run the experiments sweeping through modulation depths

# TODO : Set the number of "positions" in a variable
lead_zeros = int(np.log10(10))+1 # get the correct number of leading 0s

"""
for i, md in enumerate(mod_depth_array):
	filename = f"data/001-hard-coded-4/001.02-02_to_03_fine/eps-{str(i).zfill(lead_zeros)}.csv"
	proc = subprocess.Popen(['wsl', './parametric-modulation', filename, 
						"71.5e-9", "5.0", "150e3","1", f"{md}", "0.0"])
	proc.wait()
"""

# redo the phase-shift experiment, but across two different modulation depths
# 24 values linear spaced between 0 and 2*np.pi
phase_shift_array = np.linspace(0, 2*np.pi, 100)

# low modulation depth : 0.24
for i, phi in enumerate(phase_shift_array):
	filename = f"data/001-hard-coded-4/001.02-02_to_03_fine/lowmod-phase-{i}-rad.csv"
	proc = subprocess.Popen(['wsl', './parametric-modulation', filename, 
						"71.5e-9", "5.0", "150e3","1", "0.24", f"{phi}"])
	proc.wait()

# high modulation depth : 0.28
for i, phi in enumerate(phase_shift_array):
	filename = f"data/001-hard-coded-4/001.02-02_to_03_fine/himod-phase-{i}-rad.csv"
	proc = subprocess.Popen(['wsl', './parametric-modulation', filename, 
						"71.5e-9", "5.0", "150e3","1", "0.28", f"{phi}"])
	proc.wait()
