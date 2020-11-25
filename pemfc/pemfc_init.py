""" Initialize PEMFC model.

    For now, we hard-code all inputs.  Eventually, we may move them to a 'pemfc_inputs.py' file.
    
    Read in user inputs, then construct initial solution vector SV_0.  Create a class for all other needed parameters, called 'params.'
"""
import numpy as np

""" BEGIN USER INPUTS """
save_tag = 'porosity_study'
i_ext = 10000  # External current (A/m2)

t_final = 100000

" Thermo-chemical inputs "
T = 298   # Simulation temperature, K
P_an_0 = 101325 #Pa
X_k_an_0 = np.array([0.97, 0.03])

" Microstructure and geometry "
eps_gas_GDL = 0.7 # Gas phase volume fraction in GDL

eps_solid_CL = 0.6 # Volume fraction of solids (Pt +C) in CL
eps_gas_CL = 0.28  # Volume fraction of gas in CL

# Exponent in the Bruggeman correlatin: tau_fac = eps^alpha
alpha_Brugg_GDL = -1
alpha_Brugg_CL = -0.5

Pt_surf_frac = 0.1 # Fraction of carbon surface covered by Pt in CL.

dy_GDL = 100e-6 # GDL thickness (m)
dy_CL = 20e-6   # CL thickness (m)

# Carbon particle diameter:
d_part_GDL = 5e-6 # m
d_part_CL = 100e-9 # m

" Transport properties"
D_k_g_an = np.array([5.48e-4, 6.13e-5]) #Mixture-averaged diffusion coeffs, m2/s
mu_g_an = 9.54e-6                       # Dynamic viscsity, Pa-s

" Initial electric potential values "
phi_an_0 = 0
phi_elyte_0 = 0.6
phi_ca_0 = 1.1

" Charge transfer inputs "
C_dl_an = 1e2 # F/m2
C_dl_ca = 1e2 # F/m2

i_o_an = 2.5e-3
i_o_ca = 1e-3

n_an = -2.
n_ca = 4

F = 96485
R = 8.3145

# Charge transfer stoichiometric coefficients:
# Must be in same order as X_k_an:
nu_k_an = np.array([-1., 0.])

beta_ca = 0.5
beta_an = 0.5

delta_Phi_eq_an = -0.61
delta_Phi_eq_ca = 0.55
""" END USER INPUTS """


"============ INITIALIZE SOLUTION VECTOR ============"
C_k_an_CL_0 = P_an_0*X_k_an_0/R/T


# Construct initial solution vector
SV_0 = np.hstack((np.array([phi_an_0 - phi_elyte_0]), C_k_an_CL_0, C_k_an_CL_0,
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

    dy_GDL = dy_GDL
    dy_CL = dy_CL

    eps_g_dy_Inv_CL = 1/dy_CL/eps_gas_CL
    eps_g_dy_Inv_GDL = 1/dy_GDL/eps_gas_GDL

    eps_g_GDL = eps_gas_GDL
    eps_g_CL = eps_gas_CL
    n_Brugg_GDL = alpha_Brugg_GDL
    n_Brugg_CL = alpha_Brugg_CL

    nu_k_an = nu_k_an

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

# Create a 'pointer' class that specifies where in SV certain variables are 
#    stored:
class ptr:
    phi_dl_an = 0
    
    # C_k in anode GDL: starts just after phi_dl, is same size as X_k_an:
    C_k_an_GDL = np.arange(phi_dl_an+1, phi_dl_an+1+X_k_an_0.shape[0])
    
    # C_k in anode CL: starts just after GDL, is same size as X_k_an:
    C_k_an_CL = np.arange(C_k_an_GDL[-1]+1, C_k_an_GDL[-1]+1+X_k_an_0.shape[0])
    
    # phi_dl_ca: starts just after C_k_an_CL:
    phi_dl_ca = C_k_an_CL[-1] + 1