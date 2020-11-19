""" Residual for 1D pemfc model.
        Inputs:
            t: current simulation time (s).  Required by the integrator, but not used.
            SV: current state of the PEMFC simulation domain.
        Returns:   
            dSV_dt: time derivative of solution vector variables SV.
"""
import numpy as np
from pemfc_init import pars, ptr
from math import exp

# Constants:
F = 96485    # Faraday's constant, C/mol of equivalent
R = 8.3145   # Universal gas constant, J/mol-K

def residual(t,SV):
    
    # Initialize the residual: must be the exact same size as SV.  Initialize 
    #   as all zeros; that way, a variable will not change with time, if we do 
    #   not specify a residual.
    dSV_dt = np.zeros_like(SV)

    "========ANODE==========="
    # Calculate overpotential.
    # SHOULD BE A FUNCTION OF SPECIES ACTIVITIES:
    eta_an = SV[ptr.phi_dl_an] - pars.delta_Phi_eq_an

    # Butler-Vollmer equation:
    i_Far_an = pars.i_o_an*(exp(-pars.n_an*F*pars.beta_an*eta_an/R/pars.T)
                      - exp(pars.n_an*F*(1-pars.beta_an)*eta_an/R/pars.T))

    # Double layer current density (per unit area surface)
    i_dl_an = pars.i_ext - i_Far_an

    # Change in double layer potential per time:
    dSV_dt[ptr.phi_dl_an] = -i_dl_an/pars.C_dl_an

    # Catalyst layer gas phase
    " THESE ARE TEMPORARY "
    C_k_an_CL = SV[ptr.C_k_an_CL]
    N_k_i = np.zeros_like(C_k_an_CL)
    sdot_k = np.zeros_like(N_k_i)
    pars.A_fac = 1.
    pars.eps_dy_Inv_CL = 1.
    
    # Change in catalyst layer gas phase mole fractions:
    dCk_dt = (N_k_i + sdot_k/pars.A_fac)*pars.eps_dy_Inv_CL
    dSV_dt[ptr.C_k_an_CL] = dCk_dt
    




    "========CATHODE==========="
    eta_ca = SV[3] - pars.delta_Phi_eq_ca
    i_Far_ca = pars.i_o_ca*(exp(-pars.n_ca*F*pars.beta_ca*eta_ca/R/pars.T)
                      - exp(pars.n_ca*F*(1-pars.beta_ca)*eta_ca/R/pars.T))
    i_dl_ca = pars.i_ext - i_Far_ca
    
    dSV_dt[3] = -i_dl_ca/pars.C_dl_ca
    return dSV_dt