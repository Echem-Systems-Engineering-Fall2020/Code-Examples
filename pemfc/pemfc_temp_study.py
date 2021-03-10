""" This file runs simulated polarization curves for a range of temperatures,   
    as specified in the array `T_range`, below.

    At present, it calls the function 'polarization' for the specified temperature, then plots the curves, one at a time.
"""

# pemfc_temp_study.py
import numpy as np
from pemfc_polarization import polarization
from matplotlib import pyplot as plt

# Define temperatures to test
T_lower = 298.15 # Kelvin
T_upper = 368.15 # Kelvin
n_Temps = 5 # Number of temperatures to simulate
T_range = np.linspace(T_lower, T_upper, n_Temps)

for T in T_range:
    i, V = polarization(T)
    plt.plot(i, V)

plt.legend(T_range-273.15)
plt.xlabel('Current density (A/m$^2$)',fontsize=14)
plt.ylabel('Cell Potential (V)',fontsize=14)
plt.savefig('polarization_vs_temperature.pdf',dpi=350)
plt.show()