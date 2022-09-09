#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob
#%% PFC : GAIN SWEEPS

fb_gain_array = np.linspace(1e4,+1e5,20)
data_dir = "data/002-xvproduct/002.01-fb-sweep-3/"

control_data = pd.read_csv("data/control.csv")
control_amp = np.std(control_data.x)

fb_files = glob(data_dir + "fb*.csv")

fb_data = {fb : pd.read_csv(filename) for fb, filename in zip(fb_gain_array, fb_files)}

#%% AMPLITUDE AS A FUNCTION OF GAIN

fb_amps = np.asarray([np.std(fb_data[dat].x) for dat in fb_data.keys()])

fig, axs = plt.subplots(1,2, figsize=(14,5))

fig.suptitle("Amplitude as a function of feedback gain", fontweight="bold")
axs[0].plot(fb_gain_array, fb_amps*1e9,"o-", label="Modulated amplitude")
axs[0].axhline(control_amp*1e9, color="red", label="Control amplitude")
axs[1].semilogy(fb_gain_array, fb_amps*1e9, label="Modulated amplitude")
axs[1].axhline(control_amp*1e9, color="red", label="Control amplitude")

axs[1].legend(loc="best", bbox_to_anchor=(1.45,0.75))
axs[0].set_ylabel("Amplitude [nm]")

for i in [0,1]:
	axs[i].set_xlabel(r"Feedback gain $g_\mathrm{fb}$")
plt.savefig(data_dir + "fb-sweep-summary.png", bbox_inches="tight")
plt.show()
# %% WATERFALL PLOTS
plt.figure(figsize=(10,10))
i = 0
plt.psd(control_data.x*10, NFFT=2**15, Fs=1e3, label="Control", lw=3)
for fb, fb_dat in fb_data.items():
	print(fb_dat.x.std())
	plt.psd(fb_dat.x*10**(-i), NFFT=2**15, Fs=1e3, label=f"{fb:.2e}", alpha=0.1)

plt.legend(loc="upper right", bbox_to_anchor=(1.2,1.0))
plt.title("PSD as a function of feedback gain", fontweight="bold")
plt.xlim((125,200))
plt.xlabel("Frequency [kHz]")
plt.ylim((-269,-168))
plt.savefig(data_dir + "psd_fb_sweep-overlap.png", bbox_inches="tight")
plt.show()

# %%
