"""
    This file runs and executes a simple PEMFC Model, calculating the cell potential as a function of a user-specified current density.

    The code structure at this level is meant to be very simple, and delegates most functions to lower-level modules, which can be updated with new capabilties, over time.  The code:

        1 - Initializes the model by calling pemfc_init.py
            a - pemfc_init.py reads in the user inputs from pemfc_inputs.py, 
            b - pemfc_init.py then creates an initial solution vector 'SV_0'
            c - pemfc_init.py returns SV_0 and a class 'pars' holding all simulation parameters.
        2 - Call the function pemfc_function.py and integrate over the 
            user-defined time span from pemfc_intputs.py.  For this simulation, we want steady-state behvaior, so we simply integrate over a sufficiently long time span (100 s), with a fixed boundary conition.
        3 - The simulation then returns the solution vector at the end of the 
            integration, which can be processed as needed to plot or analyze any quantities of interest.
"""

# Import necessary modules:
from scipy.integrate import solve_ivp #integration function for ODE system.
from pemfc_function import residual # point the model to the residual function
from pemfc_init import pars, SV_0

solution = solve_ivp(residual,pars.time_span,SV_0,rtol=1e-4, atol=1e-6)


# Some initial plotting to make sure it works:
from matplotlib import pyplot as plt
for var in solution.y:
    plt.plot(solution.t,var)
    
plt.legend(['Anode double layer','Cathode double layer'])

plt.show()