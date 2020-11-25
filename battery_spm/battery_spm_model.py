""" Battery single particle model."""
from scipy.integrate import solve_ivp
import numpy as np
from battery_spm_function import residual

from battery_spm_init import SV_0, t_final, pars

time_span = np.array([0,t_final])

solution = solve_ivp(lambda t, y: residual(t, y, pars), time_span, SV_0,
    rtol=1e-6, atol=1e-8)

from matplotlib import pyplot as plt
for var in solution.y:
    plt.plot(solution.t,var)
    
plt.legend(['Anode double layer','Cathode double layer'])
plt.show()