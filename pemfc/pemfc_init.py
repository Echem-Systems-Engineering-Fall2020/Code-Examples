""" Initialize PEMFC model.

    For now, we hard-code all inputs.  Eventually, we may move them to a 'pemfc_inputs.py' file.
    
    Read in user inputs, then construct initial solution vector SV_0.  Create a class for all other needed parameters, called 'params.'
"""
import numpy as np

""" BEGIN USER INPUTS """

i_ext = 20

t_final = 100

T = 298   # Simulation temperature, K
P_an_0 = 101325 #Pa
X_k_an_0 = np.array([0.97, 0.03])

" Electric potentials"
phi_an_0 = 0
phi_elyte_0 = 0.6
phi_ca_0 = 1.1

C_dl_an = 1e2 # F/m2
C_dl_ca = 1e2 # F/m2

i_o_an = 2.5
i_o_ca = 1
n_an = 2
n_ca = 4

F = 96485
R = 8.3145

beta_ca = 0.5
beta_an = 0.5

delta_Phi_eq_an = 0.61
delta_Phi_eq_ca = 0.55
""" END USER INPUTS """


"============ INITIALIZE SOLUTION VECTOR ============"
C_k_an_CL_0 = P_an_0*X_k_an_0/R/T


# Construct initial solution vector
SV_0 = np.hstack((np.array([phi_elyte_0 - phi_an_0]), C_k_an_CL_0, 
    np.array([phi_ca_0 - phi_elyte_0])))

# Construct class containing all parameters.  Mostly we are just copying the 
#   variable names from the user inputs into the class structure:
class pars:
    time_span = np.array([0,t_final])

    T = T

    i_ext = i_ext

    # Anode
    delta_Phi_eq_an = delta_Phi_eq_an
    i_o_an = i_o_an
    n_an = n_an
    beta_an = beta_an

    C_dl_an = C_dl_an

    # Cathode
    delta_Phi_eq_ca = delta_Phi_eq_ca
    i_o_ca = i_o_ca
    n_ca = n_ca
    beta_ca = beta_ca

    C_dl_ca = C_dl_ca

# Create a 'pointer' class that specifies where in SV certain variables are 
#    stored:
class ptr:
    phi_dl_an = 0

    # C_k in anode CL: starts just after phi_dl, is same size as X_k_an:
    C_k_an_CL = np.array([phi_dl_an+1, phi_dl_an+X_k_an_0.shape[0]])

    # phi_dl_ca: starts just after C_k_an_CL:
    phi_dl_ca = C_k_an_CL[-1] + 1