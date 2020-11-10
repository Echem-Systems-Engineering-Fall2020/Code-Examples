""" Residual for 1D pemfc model.
        Inputs:
            t: current simulation time (s).  Required by the integrator, but not used.
            SV: current state of the PEMFC simulation domain.
        Returns:   
            dSV_dt: time derivative of solution vector variables SV.
"""
import numpy as np
from pemfc_init import pars
from math import exp

# Constants:
F = 96485    # Faraday's constant, C/mol of equivalent
R = 8.3145   # Universal gas constant, J/mol-K

def residual(t,SV):
    dSV_dt = np.zeros_like(SV)
    
    eta_an = SV[0] - pars.delta_Phi_eq_an
    i_Far_an = pars.i_o_an*(exp(-pars.n_an*F*pars.beta_an*eta_an/R/pars.T)
                      - exp(pars.n_an*F*(1-pars.beta_an)*eta_an/R/pars.T))
    i_dl_an = pars.i_ext - i_Far_an
    dSV_dt[0] = -i_dl_an/pars.C_dl_an
    
    
    eta_ca = SV[1] - pars.delta_Phi_eq_ca
    i_Far_ca = pars.i_o_ca*(exp(-pars.n_ca*F*pars.beta_ca*eta_ca/R/pars.T)
                      - exp(pars.n_ca*F*(1-pars.beta_ca)*eta_ca/R/pars.T))
    i_dl_ca = pars.i_ext - i_Far_ca
    
    dSV_dt[1] = -i_dl_ca/pars.C_dl_ca
    return dSV_dt