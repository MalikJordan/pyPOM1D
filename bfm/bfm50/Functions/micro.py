from bfm.bfm50.Functions.other_functions import eTq_vector, get_concentration_ratio

def microzoo_eqns(conc, microzoo_parameters, constant_parameters, environmental_parameters, zc, zn, zp, i_c, i_n, i_p, temp):
    """ Calculates the micorzooplankton (Z5 & Z6) terms needed for the zooplankton biological rate equations
        Equations come from the BFM user manual and the fortran code MicroZoo.F90
    """

    # Dissolved oxygen concentration (mg O_2 m^-3)
    disOxygen_IO_O = conc[0]
    
    # Concentration ratios
    zn_zc = get_concentration_ratio(zn, zc, constant_parameters["p_small"])
    zp_zc = get_concentration_ratio(zp, zc, constant_parameters["p_small"])

    # Temperature regulating factor
    fTZ = eTq_vector(temp, environmental_parameters["basetemp"], environmental_parameters["q10z"])
    
    # Oxygen dependent regulation factor
    fZO = min(1.0, (disOxygen_IO_O/(disOxygen_IO_O + microzoo_parameters["z_disOxygen_IO_O"])))
    
    #---------------------- Microzooplankton Respiration ----------------------
    # Zooplankton total repiration rate (eqn. 2.4.8, and matches fortran code)
    rrac = i_c*(1.0 - microzoo_parameters["etaZ"] - microzoo_parameters["betaZ"])
    rrsc = microzoo_parameters["bZ"]*fTZ*zc
    dZcdt_rsp_disInorgCarbon_IO_C = rrac + rrsc
    
    #------------- Microzooplankton mortality and activity excretion ----------
    # From fortran code MesoZoo.F90 lines 327-331
    rdc = ((1.0 - fZO)*microzoo_parameters["d_ZO"] + microzoo_parameters["d_Z"])*zc
    reac = i_c*(1.0 - microzoo_parameters["etaZ"])*microzoo_parameters["betaZ"]
    rric = reac + rdc
    dZcdt_rel_labileDOM_NO_C = rric*constant_parameters["epsilon_c"]
    dZcdt_rel_particOrganDetritus_NO_C = rric*(1.0 - constant_parameters["epsilon_c"])    

    #------------------- Microzooplankton nutrient dynamics -------------------
    # Organic Nitrogen dynamics (from fortran code) [mmol N m^-3 s^-1]
    rrin = i_n*microzoo_parameters["betaZ"] + rdc*zn_zc
    dZndt_rel_labileDOM_NO_N = rrin*constant_parameters["epsilon_n"]
    dZndt_rel_particOrganDetritus_NO_N = rrin - dZndt_rel_labileDOM_NO_N

    # Organic Phosphorus dynamics (from fortran code) [mmol P m^-3 s^-1]
    rrip = i_p*microzoo_parameters["betaZ"] + rdc*zp_zc
    dZpdt_rel_labileDOM_NO_P = rrip*constant_parameters["epsilon_p"]
    dZpdt_rel_particOrganDetritus_NO_P = rrip - dZpdt_rel_labileDOM_NO_P

    #--------------- Microzooplankton Dissolved nutrient dynamics -------------     
    # Equations from fortran code (MicroZoo.F90 line 368-371)
    runc = max(0.0, i_c*(1.0 - microzoo_parameters["betaZ"])-rrac)
    runn = max(0.0, i_n*(1.0 - microzoo_parameters["betaZ"]) + rrsc*zn_zc)
    runp = max(0.0, i_p*(1.0 - microzoo_parameters["betaZ"]) + rrsc*zp_zc)
    dZpdt_rel_phospate_IO_P = max(0.0, runp/(constant_parameters["p_small"] + runc) - microzoo_parameters["p_Zopt"])*runc
    dZndt_rel_ammonium_IO_N = max(0.0, runn/(constant_parameters["p_small"] + runc) - microzoo_parameters["n_Zopt"])*runc
    
    return dZcdt_rel_labileDOM_NO_C, dZcdt_rel_particOrganDetritus_NO_C, dZcdt_rsp_disInorgCarbon_IO_C, dZndt_rel_labileDOM_NO_N, dZndt_rel_particOrganDetritus_NO_N, dZpdt_rel_labileDOM_NO_P, dZpdt_rel_particOrganDetritus_NO_P, dZpdt_rel_phospate_IO_P, dZndt_rel_ammonium_IO_N
    
