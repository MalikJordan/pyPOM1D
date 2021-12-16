from bfm.bfm50.Functions.other_functions import eTq_vector, get_concentration_ratio

def mesozoo_eqns(conc, mesozoo_parameters, constant_parameters, environmental_parameters, zc, zn, zp, i_c, i_n, i_p, temp):
    """ Calculates the mesozooplankton (Z3 & Z4) terms needed for the zooplankton biological rate equations
        Equations come from the BFM user manual and the fortran code MesoZoo.F90
    """         
    
    # Dissolved oxygen concentration (mg O_2 m^-3)
    disOxygen_IO_O = conc[0]              # Dissolved oxygen (mg O_2 m^-3)

    # Concentration ratios
    zn_zc = get_concentration_ratio(zn, zc, constant_parameters["p_small"])
    zp_zc = get_concentration_ratio(zp, zc, constant_parameters["p_small"])
    
    # Temperature regulating factor
    fTZ = eTq_vector(temp, environmental_parameters["basetemp"], environmental_parameters["q10z"])
    
    # Oxygen dependent regulation factor
    fZO = (max(constant_parameters["p_small"],disOxygen_IO_O)**3)/(max(constant_parameters["p_small"],disOxygen_IO_O)**3 + mesozoo_parameters["z_disOxygen_IO_O"]**3)
    
    # energy cost of ingestion
    prI = 1.0 - mesozoo_parameters["etaZ"] - mesozoo_parameters["betaZ"]
    
    # Zooplankton total repiration rate (from 'MesoZoo.F90' line 343)
    dZcdt_rsp_disInorgCarbon_IO_C = prI*i_c + mesozoo_parameters["bZ"]*fTZ*zc
    
    # Specific rates of low oxygen mortality and Density dependent mortality
    # from fortran code MesoZoo.F90 lines 343-344
    rdo_c = mesozoo_parameters["d_Zdns"]*(1.0 - fZO)*fTZ*zc
    rd_c = mesozoo_parameters["d_Z"]*zc**mesozoo_parameters["gammaZ"]
    
    # Total egestion including pellet production (from MesoZoo.F90 line 359 - 361)
    dZcdt_rel_particOrganDetritus_NO_C = mesozoo_parameters["betaZ"]*i_c + rdo_c + rd_c
    dZndt_rel_particOrganDetritus_NO_N = mesozoo_parameters["betaZ"]*i_n + zn_zc*(rdo_c + rd_c)
    dZpdt_rel_particOrganDetritus_NO_P = mesozoo_parameters["betaZ"]*i_p + zp_zc*(rdo_c + rd_c)
    
    # Check the assimilation rate for Carbon, Nitrogen and Phosphorus
    # compute P:C and N:C ratios in the assimilation rate
    # from MesoZoo.F90 lines 371-375
    ru_c = mesozoo_parameters["etaZ"]*i_c
    ru_n = (mesozoo_parameters["etaZ"] + prI)*i_n
    ru_p = (mesozoo_parameters["etaZ"] + prI)*i_p
    pu_e_n = ru_n/(constant_parameters["p_small"] + ru_c)
    pu_e_p = ru_p/(constant_parameters["p_small"] + ru_c)
    
    # Eliminate the excess of the non-limiting constituent
    # Determine whether C, P or N is the limiting element and assign the value to variable limiting_nutrient
    # from MesoZoo.F90 lines 
    limiting_nutrient = 'carbon'
    temp_p = pu_e_p/(zp_zc + constant_parameters["p_small"])
    temp_n = pu_e_n/(zn_zc + constant_parameters["p_small"])
    
    if temp_p<temp_n or abs(temp_p - temp_n)<constant_parameters["p_small"]:
        if pu_e_p<zp_zc:
            limiting_nutrient = 'phosphorus'
    else:
        if pu_e_n<zn_zc:
            limiting_nutrient = 'nitrogen'
    
    # Compute the correction terms depending on the limiting constituent
    if limiting_nutrient == 'carbon':
        q_Zc = 0.0
        q_Zp = max(0.0, (1.0 - mesozoo_parameters["betaZ"])*i_p - mesozoo_parameters["p_Zopt"]*ru_c)
        q_Zn = max(0.0, (1.0 - mesozoo_parameters["betaZ"])*i_n - mesozoo_parameters["n_Zopt"]*ru_c)
    elif limiting_nutrient == 'phosphorus':
        q_Zp = 0.0
        q_Zc = max(0.0, ru_c - (1.0 - mesozoo_parameters["betaZ"])*i_p/mesozoo_parameters["p_Zopt"])
        q_Zn = max(0.0, (1.0 - mesozoo_parameters["betaZ"])*i_n - mesozoo_parameters["n_Zopt"]*(ru_c - q_Zc))
    elif limiting_nutrient == 'nitrogen':
        q_Zn = 0.0
        q_Zc = max(0.0, ru_c - (1.0 - mesozoo_parameters["betaZ"])*i_n/mesozoo_parameters["n_Zopt"])
        q_Zp = max(0.0, (1.0 - mesozoo_parameters["betaZ"])*i_p - mesozoo_parameters["p_Zopt"]*(ru_c - q_Zc))

    # Nutrient remineralization basal metabolism + excess of non-limiting nutrients
    dZpdt_rel_phospate_IO_P = mesozoo_parameters["bZ"]*fZO*fTZ*zp + q_Zp
    dZndt_rel_ammonium_IO_N = mesozoo_parameters["bZ"]*fZO*fTZ*zn + q_Zn
    
    # Fluxes to particulate organic matter 
    # Add the correction term for organic carbon release based on the limiting constituent
    dZcdt_rel_particOrganDetritus_NO_C += q_Zc
    
    # mesozooplankton are assumed to have no dissolved products
    dZcdt_rel_labileDOM_NO_C = 0.0
    dZndt_rel_labileDOM_NO_N = 0.0
    dZpdt_rel_labileDOM_NO_P = 0.0
    
    
    return dZcdt_rel_labileDOM_NO_C, dZcdt_rel_particOrganDetritus_NO_C, dZcdt_rsp_disInorgCarbon_IO_C, dZndt_rel_labileDOM_NO_N, dZndt_rel_particOrganDetritus_NO_N, dZpdt_rel_labileDOM_NO_P, dZpdt_rel_particOrganDetritus_NO_P, dZpdt_rel_phospate_IO_P, dZndt_rel_ammonium_IO_N
