import numpy as np
from pemfc_model import pemfc_model

def polarization(T=None):
    print("temperature = ", T)
    i_array = np.linspace(1,20000,25)
    V_cell = np.zeros_like(i_array)

    for j, current in enumerate(i_array):
        print('    i_ext = ',current)
        solution = pemfc_model(current,T)
        V_cell[j] = solution.y[-1,-1] - solution.y[0,-1]
        print('        V_cell =', V_cell[j])

    return i_array, V_cell

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    
    i, V = polarization()

    plt.plot(i/10000, V)
    plt.show()
