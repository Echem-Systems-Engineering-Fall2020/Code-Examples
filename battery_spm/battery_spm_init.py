from battery_spm_inputs import *
import numpy as np

phi_dl_an_0 = phi_an_0 - phi_sep_0
phi_dl_ca_0 = phi_ca_0 - phi_sep_0


capacity_anode = capacity_graphite*H_an*eps_graphite*density_graphite
capacity_cathode = capacity_LCO*H_ca*eps_LCO*density_LCO
capacity_area = min(capacity_anode,capacity_cathode)


t_final = charge_frac*3600./C_rate



SV_0 = np.array([phi_dl_an_0, phi_dl_ca_0])


# Load inputs and other parameters into 'pars' class:
class pars:
    T = T

    dPhi_eq_an = dPhi_eq_an
    dPhi_eq_ca = dPhi_eq_ca

    i_o_an = i_o_an
    n_an = n_an
    beta_an = beta_an

    i_o_ca = i_o_ca
    n_ca = n_ca
    beta_ca = beta_ca

    C_dl_an_inv = 1/C_dl_an
    C_dl_ca_inv = 1/C_dl_ca

    i_ext = C_rate*capacity_area

    A_fac_an = r_p_an/3/H_an/eps_graphite
    A_fac_ca = r_p_ca/3/H_ca/eps_LCO