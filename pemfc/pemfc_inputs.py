""" 
    Specify all user inputs for the pemfc model.

    Some of these can be overwritten by user inputs at run-time.
"""
import numpy as np

""" BEGIN USER INPUTS """
# This tag will be appended to the folder where outputs are saved.
save_tag = 'porosity_study'

i_ext = 10000  # External current (A/m2)

t_final = 100000 # How many seconds of operation to simulate. Here, we pick a suficiently long time to reach steady state for a given current density.

" Thermo-chemical inputs "
T = 298   # Simulation temperature, K
P_an_0 = 101325 #Pa
"All anode species arrays assume the order of H2, H2O"
X_k_an_0 = np.array([0.97, 0.03])
"All cathode species arrays assume the order of O2, N2"
X_k_ca_0 = np.array([0.21, 0.79])

" Microstructure and geometry--Anode and Cathode are assumed equal"
eps_gas_GDL = 0.7 # Gas phase volume fraction in GDL

eps_solid_CL = 0.6 # Volume fraction of solids (Pt +C) in CL
eps_gas_CL = 0.28  # Volume fraction of gas in CL
# Remainder is ionomer

# Exponent in the Bruggeman correlatin: tau_fac = eps^alpha
alpha_Brugg_GDL = -1
alpha_Brugg_CL = -0.5

Pt_surf_frac = 0.1 # Fraction of carbon surface covered by Pt in CL.

dy_GDL = 100e-6 # GDL thickness (m)
dy_CL = 20e-6   # CL thickness (m)
npoints_CL = 3  # Number of finite volumes/mesh points in CL (anode and cathode)

# Carbon particle diameter:
d_part_GDL = 5e-6 # m
d_part_CL = 100e-9 # m

" Transport properties"
" Properties are assumed/made up, for now (12 March 2021)"
D_k_g_an = np.array([5.48e-4, 6.13e-5]) #Mixture-averaged diffusion coeffs, m2/s
D_k_g_ca = np.array([1.48e-4, 9.13e-5]) #Mixture-averaged diffusion coeffs, m2/s
mu_g_an = 9.54e-6                       # Dynamic viscsity, Pa-s
mu_g_ca = 6.54e-6                       # Dynamic viscsity, Pa-s

" Initial electric potential values "
phi_an_0 = 0
phi_elyte_0 = 0.6
phi_ca_0 = 1.1

" Charge transfer inputs "
C_dl_an = 1e2 # F/m2
C_dl_ca = 1e2 # F/m2

i_o_an = 2.5e-1 
i_o_ca = 1e-1

n_an = -2.
n_ca = 4

F = 96485
R = 8.3145

# Charge transfer stoichiometric coefficients:
# Must be in same order as X_k_an:
nu_k_an = np.array([-1., 0.])
nu_k_ca = np.array([-0.5, 0])

beta_ca = 0.5
beta_an = 0.5

delta_Phi_eq_an = -0.61
delta_Phi_eq_ca = 0.55