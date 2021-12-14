from Functions.other_functions import eTq_vector

def pel_chem_eqns(pel_chem_parameters, environmental_parameters, constant_parameters, temper, conc, flPTreductEquiv_IO_R):
    """ calculates the non-living equations for DOM, POM, and nutrients """
    
    # State variables
    disOxygen_IO_O = conc[0]              # Dissolved oxygen (mg O_2 m^-3)
    nitrate_IO_N = conc[2]              # Nitrate (mmol N m^-3)
    ammonium_IO_N = conc[3]              # Ammonium (mmol N m^-3)
    reductEquiv_IO_R = conc[6]              # Reduction equivalents (mmol S m^-3)
    particOrganDetritus_NO_Si = conc[47]             # Particulate organic silicate (mmol Si m^-3)
    
    # Regulating factors
    eo = max(constant_parameters["p_small"], disOxygen_IO_O)/(max(constant_parameters["p_small"], disOxygen_IO_O)+ pel_chem_parameters["h_o"])
    er = reductEquiv_IO_R/(reductEquiv_IO_R + pel_chem_parameters["h_r"])
    
    # Temperature regulating factors
    fTn = eTq_vector(temper, environmental_parameters["basetemp"], environmental_parameters["q10n"])
    fTr6 = eTq_vector(temper, environmental_parameters["basetemp"], environmental_parameters["q10n5"])
    
    # Nitrification in the water  [mmol N m^-3 s^-1]   
    dammonium_IO_Ndt_nit_nitrate_IO_N = max(0.0, pel_chem_parameters["lambda_ammonium_IO_Nit"]*ammonium_IO_N*fTn*eo)

    # Denitrification flux [mmol N m^-3 s^-1] from PelChem.F90 line 134
    rPAo = flPTreductEquiv_IO_R/constant_parameters["omega_r"]
    dnitrate_IO_Ndt_denit = max(0.0, pel_chem_parameters["lambda_N3denit"]*fTn*er*rPAo/pel_chem_parameters["m_o"]*nitrate_IO_N)
    
    # Reoxidation of reduction equivalents [mmol S m^-3 s^-1]
    dreductEquiv_IO_Rdt_reox = pel_chem_parameters["lambda_reductEquiv_IO_Reox"]*eo*reductEquiv_IO_R
    
    # Dissolution of biogenic silicate [mmol Si m^-3 s^-1]
    dparticOrganDetritus_NO_Sidt_rmn_silicate_IO_Si = pel_chem_parameters["lambda_srmn"]*fTr6*particOrganDetritus_NO_Si
    
    return (dammonium_IO_Ndt_nit_nitrate_IO_N, dnitrate_IO_Ndt_denit, dreductEquiv_IO_Rdt_reox, dparticOrganDetritus_NO_Sidt_rmn_silicate_IO_Si)
