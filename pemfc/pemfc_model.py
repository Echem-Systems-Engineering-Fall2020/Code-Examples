"""
    This file runs and executes a simple PEMFC Model, calculating the cell potential as a function of a user-specified current density.

    The code structure at this level is meant to be very simple, and delegates most functions to lower-level modules, which can be updated with new capabilties, over time.  The code:

        1 - Initializes the model by calling pemfc_init.py
            a - pemfc_init.py reads in the user inputs from pemfc_inputs.py, 
            b - pemfc_init.py then creates an initial solution vector 'SV_0'
            c - pemfc_init.py returns SV_0 a class 'inputs' holding all input values, and a class 'pars' holding all parameters.
        2 - Call the function pemfc_function.py and integrate over the 
            user-defined time span from pemfc_intputs.py.  For this simulation, we want steady-state behvaior, so we simply integrate over a sufficiently long time span (100 s), with a fixed boundary conition.
        3 - The simulation then returns the solution vector at the end of the 
            integration, which can be processed as needed to plot or analyze any quantities of interest.
"""

