import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess

# pressure check
pressures = [1000,300,100,30,10,3,1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001]
# fixing the frequency at 150 kHz


def data_stats(df, title, filename, N=100):
	average_pos = np.mean(df.x)
	average_vel = np.mean(df.v)
	pos_var = np.std(df.x)**2
	vel_var = np.std(df.v)**2

	fig, axs = plt.subplots(1, 2, figsize=(18,9))
	fig.suptitle(title)
	axs[0].plot(df.t[:N], df.x[:N]/pos_var, label="Position")
	axs[0].plot(df.t[:N], df.v[:N]/vel_var, label="Velocity")
	axs[0].legend(loc="upper right")
	axs[0].set_ylabel("Position/Velocity [a.u.]")
	axs[0].set_xlabel("Time [s]")

	axs[1].psd(df.x, label="Position", Fs=1e6, NFFT=2**12)
	axs[1].psd(df.v, label="Velocity", Fs=1e6, NFFT=2**12)
	axs[1].legend(loc="upper right")
	axs[1].set_xlabel("Frequency [kHz]")

	plt.savefig(filename[:-4]+".png")

	return (average_pos, pos_var)


# First running the control
"""
compile_proc = subprocess.Popen(["wsl", "g++", "control.cpp", "-o", "control"])
compile_proc.wait()

for p in pressures:
	filename = f"pressure-lfc/control_{p}mbar.csv"
	print(f"Simulating {p} mbar dataset...")
	proc = subprocess.Popen(["wsl","./control", filename, '71.5e-9', f'{p}', '150e3', '1'])
	proc.wait()

## Performing the analysis

control_df = pd.DataFrame()
control_means = []
control_var = []

for p in pressures:
	print(f"Analysing {p} mbar dataset...")
	filename = f"pressure-lfc/control_{p}mbar.csv"
	data = pd.read_csv(filename)
	mean, var = data_stats(data, f"Control : {p} mbar dataset", filename=filename)
	control_means += [mean]
	control_var += [var]

control_df["Mean"] = control_means
control_df["Variance"] = control_var

control_df.to_csv("control-summary.csv", index=False)

# Now running the LFC

compile_proc = subprocess.Popen(["wsl", "g++", "lfc.cpp", "-o", "lfc"])
compile_proc.wait()

for p in pressures:
	filename = f"pressure-lfc/lfc_{p}mbar.csv"
	print(f"Simulating {p} mbar dataset...")
	proc = subprocess.Popen(["wsl","./lfc", filename, '-7', '1000', '71.5e-9', f'{p}', '150e3', '1'])
	proc.wait()

lfc_df = pd.DataFrame()
lfc_means = []
lfc_var = []

for p in pressures:
	print(f"Analysing {p} mbar dataset...")
	filename = f"pressure-lfc/lfc_{p}mbar.csv"
	data = pd.read_csv(filename)
	mean, var = data_stats(data, f"LFC : {p} mbar dataset", filename=filename)
	lfc_means += [mean]
	lfc_var += [var]

lfc_df["Mean"] = lfc_means
lfc_df["Variance"] = lfc_var

lfc_df.to_csv("lfc-summary.csv", index=False)
"""

# Now we are going to sweep through Gamma_FB by changing the electric field stranges

compile_proc = subprocess.Popen(["wsl", "g++", "lfc.cpp", "-o", "lfc"])
compile_proc.wait()

e_fields = [700,800,900,1000,1100,1200,1300]

# choosing 1 mbar, because that's what we can access right now

for E in e_fields:
	filename = f"pressure-lfc/lfc_fb_sweep_{E}voltsmeter.csv"
	print(f"Simulating {E} V/m dataset...")
	proc = subprocess.Popen(["wsl","./lfc", filename, '-7', f'{E}', '71.5e-9', '1.0', '150e3', '1'])
	proc.wait()

lfc_fb_df = pd.DataFrame()
lfc_fb_means = []
lfc_fb_var = []

for E in e_fields:
	print(f"Analysing {E} V/m dataset...")
	filename = f"pressure-lfc/lfc_fb_sweep_{E}voltsmeter.csv"
	data = pd.read_csv(filename)
	mean, var = data_stats(data, f"LFC : {E} volts/meter dataset", filename=filename)
	lfc_fb_means += [mean]
	lfc_fb_var += [var]

lfc_fb_df["Mean"] = lfc_fb_means
lfc_fb_df["Variance"] = lfc_fb_var

lfc_fb_df.to_csv("lfc-fb-sweep-summary.csv", index=False)


