import sys
from bfm.bfm50.Functions.other_functions import eTq_vector, get_concentration_ratio

def bacteria_eqns(conc, bacteria_parameters, constant_parameters, environmental_parameters, temper):
    """ Calculates the terms needed for the bacteria biological rate equations
        Equations come from the BFM user manual
    """
    
    # State variables
    disOxygen_IO_O = conc[0]              # Dissolved oxygen (mg O_2 m^-3)
    phospate_IO_P = conc[1]              # Phosphate (mmol P m^-3)
    ammonium_IO_N = conc[3]              # Ammonium (mmol N m^-3)
    pelBacteria_LO_C = conc[7]              # Pelagic bacteria carbon (mg C m^-3)
    pelBacteria_LO_N = conc[8]              # Pelagic bacteria nitrogen (mmol N m^-3)
    pelBacteria_LO_P = conc[9]              # Pelagic bacteria phosphate (mmol P m^-3)
    labileDOM_NO_C = conc[39]             # Labile dissolved organic carbon (mg C m^-3)
    labileDOM_NO_N = conc[40]             # Labile dissolved organic nitrogen (mmol N m^-3)
    labileDOM_NO_P = conc[41]             # Labile dissolved organic phosphate (mmol P m^-3)
    semilabileDOC_NO_C = conc[42]             # Semi-labile dissolved organic carbon (mg C m^-3)
    semirefractDOC_NO_C = conc[43]             # Semi-refractory Dissolved Organic Carbon (mg C m^-3)
    particOrganDetritus_NO_C = conc[44]             # Particulate organic carbon (mg C m^-3)
    particOrganDetritus_NO_N = conc[45]             # Particulate organic nitrogen (mmol N m^-3)
    particOrganDetritus_NO_P = conc[46]             # Particulate organic phosphate (mmol P m^-3)
    
    # concentration ratios   
    bp_bc = get_concentration_ratio(pelBacteria_LO_P, pelBacteria_LO_C, constant_parameters["p_small"])
    bn_bc = get_concentration_ratio(pelBacteria_LO_N, pelBacteria_LO_C, constant_parameters["p_small"])
    labileDOM_NO_P_labileDOM_NO_C = get_concentration_ratio(labileDOM_NO_P, labileDOM_NO_C, constant_parameters["p_small"])
    particOrganDetritus_NO_P_particOrganDetritus_NO_C = get_concentration_ratio(particOrganDetritus_NO_P, particOrganDetritus_NO_C, constant_parameters["p_small"])
    labileDOM_NO_N_labileDOM_NO_C = get_concentration_ratio(labileDOM_NO_N, labileDOM_NO_C, constant_parameters["p_small"])
    particOrganDetritus_NO_N_particOrganDetritus_NO_C = get_concentration_ratio(particOrganDetritus_NO_N, particOrganDetritus_NO_C, constant_parameters["p_small"])
    
    # Temperature effect on pelagic bacteria
    fTB = eTq_vector(temper, environmental_parameters["basetemp"], environmental_parameters["q10b"])

    # oxygen non-dimensional regulation factor[-]
    # Oxygen environment: bacteria are both aerobic and anaerobic
    f_B_O = max(constant_parameters["p_small"],disOxygen_IO_O)**3/(max(constant_parameters["p_small"],disOxygen_IO_O)**3 + bacteria_parameters["h_B_O"]**3)
    
    # external nutrient limitation
    f_B_n = ammonium_IO_N/(ammonium_IO_N + bacteria_parameters["h_B_n"])
    f_B_p = phospate_IO_P/(phospate_IO_P + bacteria_parameters["h_B_p"])
    
    # Bacteria mortality (lysis) process [mg C m^-3 s^-1]
    dBcdt_lys = (bacteria_parameters["d_0B"]*fTB + bacteria_parameters["d_B_d"]*pelBacteria_LO_C)*pelBacteria_LO_C
    dBcdt_lys_labileDOM_NO_C = dBcdt_lys*constant_parameters["epsilon_c"]
    dBcdt_lys_labileDOM_NO_N = dBcdt_lys*bn_bc*constant_parameters["epsilon_n"]
    dBcdt_lys_labileDOM_NO_P = dBcdt_lys*bp_bc*constant_parameters["epsilon_p"]
    dBcdt_lys_particOrganDetritus_NO_C = dBcdt_lys*(1.0 - constant_parameters["epsilon_c"])
    dBcdt_lys_particOrganDetritus_NO_N = dBcdt_lys*bn_bc*(1.0 - constant_parameters["epsilon_n"])
    dBcdt_lys_particOrganDetritus_NO_P = dBcdt_lys*bp_bc*(1.0 - constant_parameters["epsilon_p"])


    # Substrate availability
    if bacteria_parameters["bact_version"]==1 or bacteria_parameters["bact_version"]==2:
        # nutrient limitation (intracellular)
        nut_lim_n = min(1.0, max(0.0, bn_bc/bacteria_parameters["n_B_opt"]))         # Nitrogen
        nut_lim_p = min(1.0, max(0.0, bp_bc/bacteria_parameters["p_B_opt"]))         # Phosphorus
        f_B_n_P = min(nut_lim_n, nut_lim_p)
        
        # Potential uptake by bacteria
        potential_upt = f_B_n_P*fTB*bacteria_parameters["r_0B"]*pelBacteria_LO_C
        
        # correction of substrate quality depending on nutrient content
        f_r1_n_P = min(1.0, labileDOM_NO_P_labileDOM_NO_C/bacteria_parameters["p_B_opt"], labileDOM_NO_N_labileDOM_NO_C/bacteria_parameters["n_B_opt"])
        f_r6_n_P = min(1.0, particOrganDetritus_NO_P_particOrganDetritus_NO_C/bacteria_parameters["p_B_opt"], particOrganDetritus_NO_N_particOrganDetritus_NO_C/bacteria_parameters["n_B_opt"])
    else:
        sys.exit('This code does not support this parameterization option, only bact_version=1')
        
    # Calculate the realized substrate uptake rate depending on the type of detritus and quality
    upt_labileDOM_NO_C = (bacteria_parameters["v_B_r1"]*f_r1_n_P + bacteria_parameters["v_0B_r1"]*(1.0 - f_r1_n_P))*labileDOM_NO_C
    upt_semilabileDOC_NO_C = bacteria_parameters["v_B_r2"]*semilabileDOC_NO_C
    upt_semirefractDOC_NO_C = bacteria_parameters["v_B_r3"]*semirefractDOC_NO_C
    upt_particOrganDetritus_NO_C = bacteria_parameters["v_B_r6"]*f_r6_n_P*particOrganDetritus_NO_C
    realized_upt = constant_parameters["p_small"] + upt_labileDOM_NO_C + upt_semilabileDOC_NO_C + upt_semirefractDOC_NO_C + upt_particOrganDetritus_NO_C
    
    # Actual uptake by bacteria
    actual_upt = min(potential_upt, realized_upt)
    
    # Carbon fluxes into bacteria
    dBcdt_upt_labileDOM_NO_C = actual_upt*upt_labileDOM_NO_C/realized_upt
    dBcdt_upt_semilabileDOC_NO_C = actual_upt*upt_semilabileDOC_NO_C/realized_upt
    dBcdt_upt_semirefractDOC_NO_C = actual_upt*upt_semirefractDOC_NO_C/realized_upt
    dBcdt_upt_particOrganDetritus_NO_C = actual_upt*upt_particOrganDetritus_NO_C/realized_upt
    
    # Organic Nitrogen and Phosphrous uptake
    dBcdt_upt_labileDOM_NO_N = labileDOM_NO_N_labileDOM_NO_C*dBcdt_upt_labileDOM_NO_C
    dBcdt_upt_particOrganDetritus_NO_N = particOrganDetritus_NO_N_particOrganDetritus_NO_C*dBcdt_upt_particOrganDetritus_NO_C
    dBcdt_upt_labileDOM_NO_P = labileDOM_NO_P_labileDOM_NO_C*dBcdt_upt_labileDOM_NO_C
    dBcdt_upt_particOrganDetritus_NO_P = particOrganDetritus_NO_P_particOrganDetritus_NO_C*dBcdt_upt_particOrganDetritus_NO_C
    
    # Bacteria respiration [mc C m^-3 s^-1]
    dBcdt_rsp_disInorgCarbon_IO_C = (bacteria_parameters["gamma_B_a"] + bacteria_parameters["gamma_B_O"]*(1.0 - f_B_O))*actual_upt + bacteria_parameters["b_B"]*pelBacteria_LO_C*fTB

    # Fluxes from bacteria
    if bacteria_parameters["bact_version"]==1:
        
        # There is no Carbon excretion
        dBcdt_rel_semilabileDOC_NO_C = 0.0
        dBcdt_rel_semirefractDOC_NO_C = 0.0
        
        # Dissolved Nitrogen dynamics
        dBndt_upt_rel_ammonium_IO_N = (bn_bc - bacteria_parameters["n_B_opt"])*pelBacteria_LO_C*bacteria_parameters["v_B_n"]
            
        # Dissolved Phosphorus dynamics
        dBpdt_upt_rel_phospate_IO_P = (bp_bc - bacteria_parameters["p_B_opt"])*pelBacteria_LO_C*bacteria_parameters["v_B_p"]

    # BACT2 parameterization
    if bacteria_parameters["bact_version"]==2:
        print('This code does not support this parameterization option, only bact_version=1')
        
    # BACT3 parameterization
    if bacteria_parameters["bact_version"]==3:
        print('This code does not support this parameterization option, only bact_version=1')

    # Term needed for denitrification flux (dnitrate_IO_Ndt_denit) (from PelBac.F90 line 352)
    flPTreductEquiv_IO_R = (1.0 - f_B_O)*dBcdt_rsp_disInorgCarbon_IO_C*constant_parameters["omega_c"]*constant_parameters["omega_r"]
    
    return (dBcdt_lys_labileDOM_NO_C, dBcdt_lys_labileDOM_NO_N, dBcdt_lys_labileDOM_NO_P, dBcdt_lys_particOrganDetritus_NO_C, dBcdt_lys_particOrganDetritus_NO_N, dBcdt_lys_particOrganDetritus_NO_P, 
            dBcdt_upt_labileDOM_NO_C, dBcdt_upt_particOrganDetritus_NO_C, dBpdt_upt_rel_phospate_IO_P, dBndt_upt_rel_ammonium_IO_N, dBcdt_upt_semilabileDOC_NO_C, dBcdt_upt_semirefractDOC_NO_C, 
            dBcdt_rel_semilabileDOC_NO_C, dBcdt_rel_semirefractDOC_NO_C, dBcdt_rsp_disInorgCarbon_IO_C, flPTreductEquiv_IO_R, f_B_O, f_B_n, f_B_p)
    
