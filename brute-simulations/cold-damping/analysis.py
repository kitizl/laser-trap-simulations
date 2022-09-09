import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


pressures = [1000,300,100,30,10,3,1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001]
control_data = pd.read_csv("pressure-lfc-2/control-summary.csv")
lfc_data = pd.read_csv("pressure-lfc-2/lfc-summary.csv")

#%%

print(control_data.head())
print(lfc_data.head())