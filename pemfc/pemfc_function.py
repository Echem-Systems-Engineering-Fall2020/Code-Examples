""" Residual for 1D pemfc model.
        Inputs:
            t: current simulation time (s).  Required by the integrator, but not used.
            SV: current state of the PEMFC simulation domain.
        Returns:   
            dSV_dt: time derivative of solution vector variables SV.
"""
import numpy as np
from math import exp

# Constants:
F = 96485    # Faraday's constant, C/mol of equivalent
R = 8.3145   # Universal gas constant, J/mol-K

def residual(t, SV, pars, ptr):
    # Initialize the residual: must be the exact same size as SV.  Initialize 
    #   as all zeros; that way, a variable will not change with time, if we do 
    #   not specify a residual.
    dSV_dt = np.zeros_like(SV)

    "========ANODE==========="
    # Calculate overpotential.
    # SHOULD BE A FUNCTION OF SPECIES ACTIVITIES:
    eta_an = SV[ptr.phi_dl_an] - pars.delta_Phi_eq_an 

    # Butler-Vollmer equation:
    i_Far_an = -pars.i_o_an*(exp(-pars.n_an*F*pars.beta_an*eta_an/R/pars.T)
                      - exp(pars.n_an*F*(1-pars.beta_an)*eta_an/R/pars.T))

    # Double layer current density (per unit area surface)
    i_dl_an = -pars.i_ext*pars.A_fac_dl - i_Far_an*pars.f_Pt

    # Change in double layer potential per time:
    dSV_dt[ptr.phi_dl_an] = -i_dl_an/pars.C_dl_an

    # GDL gas phase:
    C_k_an_GDL = SV[ptr.C_k_an_GDL]
    # Catalyst layer gas phase
    C_k_an_CL = SV[ptr.C_k_an_CL]
   
    s1 = {'C_k': C_k_an_GDL, 'dy':pars.dy_GDL, 'eps_g':pars.eps_g_GDL, 
        'n_Brugg':pars.n_Brugg_GDL, 'd_solid':pars.d_solid_GDL}
    s2 = {'C_k': C_k_an_CL, 'dy':pars.dy_CL, 'eps_g':pars.eps_g_CL,
        'n_Brugg':pars.n_Brugg_CL, 'd_solid':pars.d_solid_CL}
    gas_props = {'T':pars.T, 'D_k':pars.D_k_g_an, 'mu':pars.mu_g_an}
    
    N_k_i = pemfc_gas_flux(s1, s2, gas_props)
    # print(N_k_i)
    
    # Molar production rates due to Faradaic current:
    sdot_k = i_Far_an*pars.nu_k_an/pars.n_an/F
    
    # Change in catalyst layer gas phase mole fractions:
    dCk_dt = (N_k_i + sdot_k*pars.A_fac_Pt)*pars.eps_g_dy_Inv_CL
    dSV_dt[ptr.C_k_an_CL] = dCk_dt
    


    "========CATHODE==========="
    " THIS SHOULD EVENTUALLY DEPEND ON SPECIES ACTIVITIES. "
    eta_ca = SV[ptr.phi_dl_ca] - pars.delta_Phi_eq_ca
    i_Far_ca = pars.i_o_ca*(exp(-pars.n_ca*F*pars.beta_ca*eta_ca/R/pars.T)
                      - exp(pars.n_ca*F*(1-pars.beta_ca)*eta_ca/R/pars.T))
    i_dl_ca = pars.i_ext*pars.A_fac_dl - i_Far_ca*pars.f_Pt
    
    dSV_dt[ptr.phi_dl_ca] = -i_dl_ca/pars.C_dl_ca
    return dSV_dt


def pemfc_gas_flux(node1, node2, gas_props):
    N_k  = np.zeros_like(node1['C_k'])

    f1 = node1['dy']/(node1['dy'] + node2['dy'])
    f2 = 1-f1

    C_int = f1*node1['C_k'] + f2*node2['C_k']

    X_k_1 = node1['C_k']/np.sum(node1['C_k'])
    X_k_2 = node2['C_k']/np.sum(node2['C_k'])
    X_k_int = f1*X_k_1 + f2*X_k_2

    P_1 = np.sum(node1['C_k'])*R*gas_props['T']
    P_2 = np.sum(node2['C_k'])*R*gas_props['T']
    # print(P_1, P_2)

    eps_g = f1*node1['eps_g'] + f2*node2['eps_g']
    tau_fac = (f1*node1['eps_g']**node1['n_Brugg'] 
        + f2*node2['eps_g']**node2['n_Brugg'])
    D_k_eff = eps_g*gas_props['D_k']/tau_fac
    # print(f2, f1)
    # print(tau_fac)
    # print(f1*node1['eps_g']**node1['n_Brugg'])
    
    d_part = f1*node1['d_solid'] + f2*node2['d_solid']
    K_g = eps_g**3*d_part**2*tau_fac**(-2)*(1-eps_g)**(-2)/72

    dY = 0.5*(node1['dy'] + node2['dy'])

    V_conv = -K_g*(P_2 - P_1)/dY/gas_props['mu']
    V_k_diff = -D_k_eff*(X_k_2 - X_k_1)/dY/X_k_int

    V_k  = V_conv + V_k_diff

    N_k = C_int*X_k_int*V_k

    return N_k