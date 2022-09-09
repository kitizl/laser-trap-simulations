import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess

# =============================================================================
# PARAMETRIC COOLING
# The parameters we can adjust so far :
# 1. Feedback gain
# 2. Pressure
# Let's first find "ideal feedback gain" and then find the ideal pressure to start, k?

# we need _much_ larger values mate.
fb_gain_array = np.linspace(1e4,+1e5,20)
data_dir = "data/002-xvproduct/002.01-fb-sweep-3/"

lead_zeros = int(np.log10(20))+1 # get the correct number of leading 0s

# setting pressure to be 5 mbar, because particle is _definitely_ ballistic in that regime

for i, fb_gain in enumerate(fb_gain_array):
	filename = f"{data_dir}/fb-{str(i).zfill(lead_zeros)}.csv"
	proc = subprocess.Popen(['wsl', './parametric-cooling', filename, 
						"71.5e-9", "5.0", "150e3","1", f"{fb_gain}"])
	proc.wait()
