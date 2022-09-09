#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob
#%% Parametric Resonance Tests

mod_depth_array = np.linspace(0.2, 0.3, 100)
phase_shift_array = np.linspace(0, 2*np.pi, 10)

data_dir = "data/001-hard-coded-4/001.02-02_to_03_fine/"
control_data = pd.read_csv("data/control.csv")
eps_files = glob(data_dir + "eps-*.csv")
phase_files = glob(data_dir + "*phase-*.csv")
print(phase_files)
eps_data = {md : pd.read_csv(filepath) for md, filepath in zip(mod_depth_array, eps_files)}

phase_data = {phi : pd.read_csv(filepath) for phi, filepath in zip(phase_shift_array, phase_files)}

# the files are ordered in the same way the arrays are.

# %% Analysing Control Data

control_amp = np.std(control_data.x)
print(control_amp*1e9)

# %%

eps_amps = np.asarray([np.std(eps_data[dat].x) for dat in eps_data])
phase_amps = np.asarray([np.std(phase_data[dat].x) for dat in phase_data])
print(eps_amps)

# %%
fig, axs = plt.subplots(1, 2, figsize=(10,5))
fig.suptitle("Effect of Parametric Feedback As a Function of Modulation Depth", fontweight="bold")

axs[0].set_title("Linear plot")
axs[0].plot(mod_depth_array, eps_amps, "o-", label="Parametrically modulated amplitude")
axs[0].axhline(control_amp, color="r", label="Control amplitude")
axs[0].set_xlabel("Modulation depth parameter")
axs[0].set_ylabel("Particle amplitude [m]")

axs[1].set_title("Log plot")
axs[1].semilogy(mod_depth_array, eps_amps, "o-", label="Parametrically modulated amplitude")
axs[1].axhline(control_amp, color="r", label="Control amplitude")
axs[1].set_xlabel("Modulation depth parameter")
axs[1].set_ylabel("Particle amplitude [m]")
axs[1].legend(loc="best", bbox_to_anchor=(2.0,0.5))

plt.savefig(data_dir + "/mod-depth-sweep.png", bbox_inches="tight")
plt.show()

# %%

fig, axs = plt.subplots(1, figsize=(5,5))
fig.suptitle("Effect of Parametric Feedback As a Function of Phase Shift", fontweight="bold")

axs.set_title("Linear plot")
axs.plot(phase_shift_array*180/np.pi, phase_amps*1e9, label="Parametrically modulated amplitude")
axs.axhline(control_amp*1e9, color="r", label="Control amplitude")
axs.set_xlabel("Phase shift [degrees]")
axs.set_ylabel("Particle amplitude [nm]")
axs.legend(bbox_to_anchor=(1.9,0.5))
plt.savefig(data_dir + "phase-shift-sweep.png", bbox_inches="tight")
plt.show()
# %% Interpolate and find the maxima of the phase-shift thing

from scipy.interpolate import interp1d

f = interp1d(phase_shift_array, phase_amps, kind="cubic")
th_phase = np.linspace(0,np.pi, 1000)
angle = th_phase[np.where(f(th_phase) == np.max(f(th_phase)))[0]]
print(angle*180/np.pi)
# %% PSD of all the modulation data!!! WATERFALL PLOT, SOMETHING???

plt.figure(figsize=(5,10))
plt.title("PSD as a function of modulation depth (down to up)", fontweight="bold")
for i, eps in enumerate(eps_data):
	plt.psd(eps_data[eps].x*10**i, Fs=1e3, NFFT=2**10, label=f"Mod depth : {mod_depth_array[i]*100:.1f}%")
plt.xlabel("Frequency [kHz]")
plt.legend(loc='best', bbox_to_anchor=(1.8,0.75))
plt.savefig(data_dir + "psd-mod-sweep.png", bbox_inches="tight")
plt.show()

plt.figure(figsize=(5,10))
plt.title("PSD as a function of phase shift (down to up)", fontweight="bold")
for i, phi in enumerate(phase_data):
	plt.psd(phase_data[phi].x*10**i, Fs=1e3, NFFT=2**10, label=f"Phase shift : {phase_shift_array[i]*180/np.pi:.0f}")
plt.xlabel("Frequency [kHz]")
plt.legend(loc='best', bbox_to_anchor=(1.8,0.75))
plt.savefig(data_dir + "psd-phase-sweep.png", bbox_inches="tight")
plt.show()

# %% FINE-SWEEP PHASE, FOR TWO MODULATIONS


data_dir = r"data\001-hard-coded-4\001.02-02_to_03_fine"
phase_shift_array = np.linspace(0, 2*np.pi, 100)

himod_phase_files = glob(data_dir + "\\himod-*.csv")
lomod_phase_files = glob(data_dir + "\\lowmod-*.csv")

himod_phase_data = {phi : pd.read_csv(filepath) for phi, filepath in zip(phase_shift_array, himod_phase_files)}
lomod_phase_data = {phi : pd.read_csv(filepath) for phi, filepath in zip(phase_shift_array, lomod_phase_files)}

himod_phase_amps = np.asarray([np.std(himod_phase_data[dat].x) for dat in himod_phase_data])
lomod_phase_amps = np.asarray([np.std(lomod_phase_data[dat].x) for dat in lomod_phase_data])

# %% PLOTTING THE ABOVE

fig, axs = plt.subplots(1, 2, figsize=(10,5))
fig.suptitle("Parametric response as a function of phase shift", fontweight="bold")

axs[0].set_title("24% modulation")
axs[1].set_title("28% modulation")

axs[0].plot(phase_shift_array*180/np.pi, lomod_phase_amps*1e9, "o-", label="Modulated amplitude")
axs[0].axhline(control_amp*1e9, label="Control amplitude", color="red")

axs[1].plot(phase_shift_array*180/np.pi, himod_phase_amps*1e9, label="Modulated amplitude")
axs[1].axhline(control_amp*1e9, label="Control amplitude", color="red")

axs[0].set_ylabel("Amplitude [nm]")

for i in [0,1]:
	axs[i].set_xlabel("Phase shift [degrees]")

axs[1].legend(loc="best", bbox_to_anchor=(1.6,0.75))

plt.savefig(data_dir + "phase-amp-sweep.png", bbox_inches="tight")
plt.show()
# %% WATERFALL PLOT TO COMPARE THE TWO


fig, axs = plt.subplots(1, 2, figsize=(8,30))


axs[0].set_title("24% modulation")
axs[1].set_title("28% modulation")

fig.suptitle("PSD as a function of phase shift (down to up)", fontweight="bold")

for i, phi in enumerate(himod_phase_data):
	hi_data = himod_phase_data[phi].x.to_numpy()
	lo_data = lomod_phase_data[phi].x.to_numpy()
	axs[0].psd(hi_data*10**i, Fs=1e3, NFFT=2**10, 
	label=f"Phase shift : {phase_shift_array[i]*180/np.pi:.2f}")
	axs[1].psd(lo_data*10**i, Fs=1e3, NFFT=2**10, 
	label=f"Phase shift : {phase_shift_array[i]*180/np.pi:.2f}")
	

plt.xlabel("Frequency [kHz]")
plt.legend(loc='best', bbox_to_anchor=(1.4,1.0))
plt.savefig(data_dir + "psd-phase-sweep.png", bbox_inches="tight")
plt.show()

# %%
