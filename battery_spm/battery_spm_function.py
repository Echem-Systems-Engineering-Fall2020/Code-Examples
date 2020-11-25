import numpy as np
from math import exp

# Constants
F = 96485
R = 8.3145

def residual(t, SV, pars):

    RTinv = 1/R/pars.T
    dSV_dt = np.zeros_like(SV)
    
    eta_an = SV[0] - pars.dPhi_eq_an
    i_Far_an = pars.i_o_an*(exp(-pars.n_an*F*pars.beta_an*eta_an*RTinv)
                      - exp(pars.n_an*F*(1-pars.beta_an)*eta_an*RTinv))
    i_dl_an = pars.i_ext*pars.A_fac_an - i_Far_an
    dSV_dt[0] = i_dl_an*pars.C_dl_an_inv
    
    
    eta_ca = SV[1] - pars.dPhi_eq_ca
    i_Far_ca = pars.i_o_ca*(exp(-pars.n_ca*F*pars.beta_ca*eta_ca*RTinv)
                      - exp(pars.n_ca*F*(1-pars.beta_ca)*eta_ca*RTinv))
    i_dl_ca = -pars.i_ext*pars.A_fac_ca - i_Far_ca
    
    
    dSV_dt[1] = i_dl_ca*pars.C_dl_ca_inv
    
    return dSV_dt