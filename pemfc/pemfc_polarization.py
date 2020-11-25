import numpy as np
from pemfc_model import pemfc_model


i_array = np.linspace(0.00000001,100000,50)
V_cell = np.zeros_like(i_array)

for j, current in enumerate(i_array):
    print(current)
    solution = pemfc_model(current)
    V_cell[j] = solution.y[-1,-1] - solution.y[0,-1]
    print(V_cell[j])


from matplotlib import pyplot as plt

plt.plot(i_array,V_cell,'.')
plt.savefig('results/polarization.png',dpi=350)
plt.show()