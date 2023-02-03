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

    # GDL gas phase:
    j = 0
    s1, _  = read_state(SV, ptr, pars, 'gdl_an', j)
    
    # Read the state in the first CL node, save as state s2:
    s2, eta_2 = read_state(SV, ptr, pars, 'cl_an', j)

    gas_props = {'T':pars.T, 'D_k':pars.D_k_g_an, 'mu':pars.mu_g_an}

    # Calculate the flux between states 1 (GDL) and 2 (first CL node)
    #   The subscript 'o' refers to the flux 'out' of the s2 node.
    N_k_o = pemfc_gas_flux(s1, s2, gas_props)
    "========ANODE CATALYST LAYER==========="
    for j in np.arange(pars.npoints_CL-1):
        # Re-assign state 2 's2' as 's1', previously 'out' flux as 'in' flux.
        s1 = s2
        N_k_i = N_k_o
        eta = eta_2
    
        # Calculate currents and chemical production rates:
        i_dl, sdot_k = chem_calcs(eta, pars, 'an')
       
        # Change in double layer potential per time:
        dSV_dt[ptr.phi_dl_an[j]] = -i_dl/pars.C_dl_an

        # Gas properties used in flux calculations:
        gas_props = {'T':pars.T, 'D_k':pars.D_k_g_an, 'mu':pars.mu_g_an}

        # Read the state in the next CL node, save as state s2:
        s2, eta_2 = read_state(SV, ptr, pars, 'cl_an', j+1)
        
        # Calculate fluxes between nodes 1 and 2:
        N_k_o = pemfc_gas_flux(s1, s2, gas_props)
        
        # Change in catalyst layer gas phase mole fractions:
        dCk_dt = (N_k_i - N_k_o + sdot_k*pars.A_fac_Pt)*pars.eps_g_dy_Inv_CL
        dSV_dt[ptr.C_k_an_CL[j,:]] = dCk_dt
    
    # Re-assign state 2 's2' as 's1'
    s1 = s2
    N_k_i = N_k_o
    eta = eta_2

    # Look at the final CL node:
    j = pars.npoints_CL-1

    # Calculate currents and chemical production rates:
    i_dl, sdot_k = chem_calcs(eta, pars, 'an')

    # Change in double layer potential per time:
    dSV_dt[ptr.phi_dl_an[j]] = -i_dl/pars.C_dl_an
    
    # Change in catalyst layer gas phase mole fractions:
    dCk_dt = (N_k_i + sdot_k*pars.A_fac_Pt)*pars.eps_g_dy_Inv_CL
    dSV_dt[ptr.C_k_an_CL[j,:]] = dCk_dt

    "========CATHODE==========="
    " THIS SHOULD EVENTUALLY DEPEND ON SPECIES ACTIVITIES. "
    eta = SV[ptr.phi_dl_ca] - pars.delta_Phi_eq_ca
    
    i_dl, _ = chem_calcs(eta, pars, 'ca')

    dSV_dt[ptr.phi_dl_ca] = -i_dl/pars.C_dl_ca

    # Return the SV derivatives w/r/t time:
    return dSV_dt

def read_state(SV, ptr, pars, domain, j):
    if domain=='gdl_an':
        # No ionomer, so no eta:
        eta = 0

        # Read out gas phase molar concentrations:
        C_k = SV[ptr.C_k_an_GDL[j]]
        # Load the GDL state as s1
        state = {'C_k': C_k, 'dy':pars.dy_GDL, 'eps_g':pars.eps_g_GDL, 
            'n_Brugg':pars.n_Brugg_GDL, 'd_solid':pars.d_solid_GDL}
        
    elif domain=='cl_an':
        # Calculate overpotential.
        # SHOULD BE A FUNCTION OF SPECIES ACTIVITIES:
        eta = SV[ptr.phi_dl_an[j]] - pars.delta_Phi_eq_an 

        # Read out gas phase molar concentrations:
        C_k = SV[ptr.C_k_an_CL[j]]
        # Load the GDL state as s1
        state = {'C_k': C_k, 'dy':pars.dy_CL, 'eps_g':pars.eps_g_CL, 
            'n_Brugg':pars.n_Brugg_CL, 'd_solid':pars.d_solid_CL}

    return state, eta

def chem_calcs(eta, pars, domain):
    if domain == 'an':
        # Butler-Vollmer equation:
        # i_o SHOULD VARY WITH SPECIES ACTIVITIES:
        i_Far = -pars.i_o_an*(exp(-pars.n_an*F*pars.beta_an*eta/R/pars.T)
                        - exp(pars.n_an*F*(1-pars.beta_an)*eta/R/pars.T))

        # Molar production rates due to Faradaic current:
        sdot_k = i_Far*pars.nu_k_an/pars.n_an/F

        # Double layer current density (per unit area surface)
        i_dl = -pars.i_ext*pars.A_fac_dl - i_Far*pars.f_Pt
    elif domain=='ca':
        i_Far = pars.i_o_ca*(exp(-pars.n_ca*F*pars.beta_ca*eta/R/pars.T)
                          - exp(pars.n_ca*F*(1-pars.beta_ca)*eta/R/pars.T))
        # Molar production rates due to Faradaic current:
        sdot_k = i_Far*pars.nu_k_ca/pars.n_ca/F

        # Double layer current density (per unit area surface)
        i_dl = pars.i_ext*pars.A_fac_dl - i_Far*pars.f_Pt

    return i_dl, sdot_k

def pemfc_gas_flux(node1, node2, gas_props):
    # Initialize molar fluxes:
    N_k  = np.zeros_like(node1['C_k'])

    # Weighting fractions between current (GDL) and next (CL) node:
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
    
    d_part = f1*node1['d_solid'] + f2*node2['d_solid']
    K_g = eps_g**3*d_part**2*tau_fac**(-2)*(1-eps_g)**(-2)/72

    dY = 0.5*(node1['dy'] + node2['dy'])

    V_conv = -K_g*(P_2 - P_1)/dY/gas_props['mu']
    V_k_diff = -D_k_eff*(X_k_2 - X_k_1)/dY/X_k_int

    V_k  = V_conv + V_k_diff

    N_k = C_int*X_k_int*V_k

    return N_k