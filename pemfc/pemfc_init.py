""" Initialize PEMFC model.

    1. Read in user inputs.
    2. Construct initial solution vector SV_0.
    3. Create a class for all other needed parameters, called 'pars.' 
    4. Create a class containing all variable pointers, called 'ptr.'
"""
import numpy as np

"============ READ IN USER INPUTS ============"
from pemfc_inputs import *

"============ INITIALIZE SOLUTION VECTOR ============"
C_k_an_0 = P_an_0*X_k_an_0/R/T


# Construct initial solution vector
SV_0_GDL_an = C_k_an_0
SV_0_CL_an = np.tile(np.hstack((np.array([phi_an_0 - phi_elyte_0]),C_k_an_0)),
    npoints_CL)
SV_0_ca = np.array([phi_ca_0 - phi_elyte_0])
SV_0 = np.hstack((SV_0_GDL_an, SV_0_CL_an, SV_0_ca))

"============ CREATE PARAMETERS CLASS ============"
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

    dy_GDL = dy_GDL
    dy_CL = dy_CL/npoints_CL
    npoints_CL = npoints_CL

    eps_g_dy_Inv_CL = 1/dy_CL/eps_gas_CL
    eps_g_dy_Inv_GDL = 1/dy_GDL/eps_gas_GDL

    eps_g_GDL = eps_gas_GDL
    eps_g_CL = eps_gas_CL
    n_Brugg_GDL = alpha_Brugg_GDL
    n_Brugg_CL = alpha_Brugg_CL

    nu_k_an = nu_k_an
    nu_k_ca = nu_k_ca

    X_k_GDL = X_k_an_0

    D_k_g_an = D_k_g_an
    mu_g_an = mu_g_an

    # Pt surface area per unit geometric area:
    A_fac_Pt = 0.5*eps_solid_CL*3.*dy_CL*Pt_surf_frac/d_part_CL
    # Geometric area per unit double layer area:
    A_fac_dl = Pt_surf_frac/A_fac_Pt

    # Fraction of Carbon surface covered by Pt:
    f_Pt = Pt_surf_frac

    d_solid_GDL = d_part_GDL
    d_solid_CL = d_part_CL
    
    # Cathode
    delta_Phi_eq_ca = delta_Phi_eq_ca
    i_o_ca = i_o_ca
    n_ca = n_ca
    beta_ca = beta_ca

    C_dl_ca = C_dl_ca

"============ CREATE CLASS OF POINTERS ============"
# Create a 'pointer' class that specifies where in SV certain variables are 
#    stored:
class ptr:
    
    # C_k in anode GDL: starts just after phi_dl, is same size as X_k_an:
    C_k_an_GDL = np.arange(0, X_k_an_0.shape[0])
    
    # C_k in anode CL: starts just after GDL, is same size as X_k_an:
    phi_dl_an = np.empty((npoints_CL),int)
    C_k_an_CL = np.empty((npoints_CL, X_k_an_0.shape[0]),int)
    nvars_an_CL = 1+X_k_an_0.shape[0]
    
    for j in np.arange(npoints_CL):
        phi_dl_an[j] = C_k_an_GDL[-1] + 1 + j*nvars_an_CL
        C_k_an_CL[j,:] = np.arange(phi_dl_an[j]+1, phi_dl_an[j]+1+X_k_an_0.shape[0])
    
    # phi_dl_ca: starts just after C_k_an_CL:
    phi_dl_ca = C_k_an_CL[-1,-1] + 1
