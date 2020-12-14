# pemfc_temp_study.py
import numpy as np
from pemfc_polarization import polarization
from matplotlib import pyplot as plt

# Define temperatures to test
T_range = np.linspace(298.15,368.15,5)

for T in T_range:
    print('T = ',T)
    i, V = polarization(T)
    plt.plot(i, V)

plt.legend(T_range-273.15)
plt.xlabel('Current density (A/m$^2$)',fontsize=14)
plt.ylabel('Cell Potential (V)',fontsize=14)
plt.savefig('polarization_vs_temperature.pdf',dpi=350)
plt.show()