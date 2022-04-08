#!python3
#%%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

test_data = pd.read_csv("test/sine-driving.csv")
print(test_data.head())
# %%

plt.psd(test_data.x, Fs=1e6, NFFT=2**15)
plt.axvline(150e3, color="r")
plt.show()
# %%
plt.plot(test_data.t)
plt.show()
# %%
