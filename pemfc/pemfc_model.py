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

def pemfc_model(i_ext=None, temp=None):
    # Import necessary modules:
    from scipy.integrate import solve_ivp #integration function for ODE system.
    from pemfc_function import residual # point to the residual function

    # Import the initial solution, parameters, and pointer structures:
    from pemfc_init import SV_0, pars, ptr

    # Parse and overwrite any variables passed to the function call:
    if i_ext:
        pars.i_ext = i_ext
    if temp:
        # Adjust ideal gas concentrations in SV_0:
        SV_0[ptr.C_k_an_GDL] *= pars.T/temp
        SV_0[ptr.C_k_an_CL] *= pars.T/temp
        # Overwrite temperature:
        pars.T = temp

    # The use of the 'lambda' function is required here so that we can pass the 
    #   class variablels 'pars' and 'ptr.'  Otherwise, we can only pass the 
    #   time span and our initial solution SV_0:
    solution = solve_ivp(lambda t, y: residual(t, y, pars, ptr),
        pars.time_span, SV_0, rtol=1e-9, atol=1e-7, method='BDF')

    # Return the solution results to whatever routine called the function:
    return solution


# If you want to run this as a script, doing one-off simulations, this will 
#    call the pemfc_model function, above
if __name__ == '__main__':
    pemfc_model()