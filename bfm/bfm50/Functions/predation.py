from Functions.other_functions import eTq_vector, get_concentration_ratio

def get_mesozoo_predation_terms(conc, mesozoo3_parameters, mesozoo4_parameters, zoo_availability_parameters, environmental_parameters, constant_parameters, temp):
    """ Calculates the predation terms for mesozooplankton """
    
    # Species concentrations
    diatoms_LO_C = conc[10]             # Diatoms carbon (mg C m^-3)
    diatoms_LO_N = conc[11]             # Diatoms nitrogen (mmol N m^-3)
    diatoms_LO_P = conc[12]             # Diatoms phosphate (mmol P m^-3)
    diatoms_LO_Chl = conc[13]             # Diatoms chlorophyll (mg Chl-a m^-3)
    diatoms_LO_Si = conc[14]             # Diatoms silicate (mmol Si m^-3) 
    nanoflagellates_LO_C = conc[15]             # NanoFlagellates carbon (mg C m^-3)
    nanoflagellates_LO_N = conc[16]             # NanoFlagellates nitrogen (mmol N m^-3)
    nanoflagellates_LO_P = conc[17]             # NanoFlagellates phosphate (mmol P m^-3)
    nanoflagellates_LO_Chl = conc[18]             # NanoFlagellates chlorophyll (mg Chl-a m^-3)
    picophyto_LO_C = conc[19]             # Picophytoplankton carbon (mg C m^-3)
    picophyto_LO_N = conc[20]             # Picophytoplankton nitrogen (mmol N m^-3)
    picophyto_LO_P = conc[21]             # Picophytoplankton phosphate (mmol P m^-3)
    picophyto_LO_Chl = conc[22]             # Picophytoplankton chlorophyll (mg Chl-a m^-3)
    largephyto_LO_C = conc[23]             # Large phytoplankton carbon (mg C m^-3)
    largephyto_LO_N = conc[24]             # Large phytoplankton nitrogen (mmol N m^-3)
    largephyto_LO_P = conc[25]             # Large phytoplankton phosphate (mmol P m^-3) 
    largephyto_LO_Chl = conc[26]             # Large phytoplankton chlorophyll (mg Chl-a m^-3)
    carnivMesozoo_LO_C = conc[27]             # Carnivorous mesozooplankton carbon (mg C m^-3)
    carnivMesozoo_LO_N = conc[28]             # Carnivorous mesozooplankton nitrogen (mmol N m^-3)
    carnivMesozoo_LO_P = conc[29]             # Carnivorous mesozooplankton phosphate (mmol P m^-3)
    omnivMesozoo_LO_C = conc[30]             # Omnivorous mesozooplankton carbon (mg C m^-3)
    omnivMesozoo_LO_N = conc[31]             # Omnivorous mesozooplankton nitrogen (mmol N m^-3)
    omnivMesozoo_LO_P = conc[32]             # Omnivorous mesozooplankton phosphate (mmol P m^-3)
    microzoo_LO_C = conc[33]             # Microzooplankton carbon (mg C m^-3)
    microzoo_LO_N = conc[34]             # Microzooplankton nitrogen (mmol N m^-3)
    microzoo_LO_P = conc[35]             # Microzooplankton phosphate (mmol P m^-3)
    heteroFlagellates_LO_C = conc[36]             # Heterotrophic flagellates carbon (mg C m^-3)
    heteroFlagellates_LO_N = conc[37]             # Heterotrophic flagellates nitrogen (mmol N m^-3)
    heteroFlagellates_LO_P = conc[38]             # Heterotrophic flagellates phosphate (mmol P m^-3)
    
    # concentration ratios
    conc_ratio_n = {
            "p1": get_concentration_ratio(diatoms_LO_N, diatoms_LO_C, constant_parameters["p_small"]),
            "p2": get_concentration_ratio(nanoflagellates_LO_N, nanoflagellates_LO_C, constant_parameters["p_small"]),
            "p3": get_concentration_ratio(picophyto_LO_N, picophyto_LO_C, constant_parameters["p_small"]),
            "p4": get_concentration_ratio(largephyto_LO_N, largephyto_LO_C, constant_parameters["p_small"]),
            "z3": get_concentration_ratio(carnivMesozoo_LO_N, carnivMesozoo_LO_C, constant_parameters["p_small"]),
            "z4": get_concentration_ratio(omnivMesozoo_LO_N, omnivMesozoo_LO_C, constant_parameters["p_small"]),
            "z5": get_concentration_ratio(microzoo_LO_N, microzoo_LO_C, constant_parameters["p_small"]),
            "z6": get_concentration_ratio(heteroFlagellates_LO_N, heteroFlagellates_LO_C, constant_parameters["p_small"])
            }
    conc_ratio_p = {
            "p1": get_concentration_ratio(diatoms_LO_P, diatoms_LO_C, constant_parameters["p_small"]),
            "p2": get_concentration_ratio(nanoflagellates_LO_P, nanoflagellates_LO_C, constant_parameters["p_small"]),
            "p3": get_concentration_ratio(picophyto_LO_P, picophyto_LO_C, constant_parameters["p_small"]),
            "p4": get_concentration_ratio(largephyto_LO_P, largephyto_LO_C, constant_parameters["p_small"]),
            "z3": get_concentration_ratio(carnivMesozoo_LO_P, carnivMesozoo_LO_C, constant_parameters["p_small"]),
            "z4": get_concentration_ratio(omnivMesozoo_LO_P, omnivMesozoo_LO_C, constant_parameters["p_small"]),
            "z5": get_concentration_ratio(microzoo_LO_P, microzoo_LO_C, constant_parameters["p_small"]),
            "z6": get_concentration_ratio(heteroFlagellates_LO_P, heteroFlagellates_LO_C, constant_parameters["p_small"])
            }
    
    # Zooplankton temperature regulating factor
    fTZ = eTq_vector(temp, environmental_parameters["basetemp"], environmental_parameters["q10z"])
    
    # Calculate total potential food given the non-dim prey availability
    # There is no parameter for capture efficiency in mesozooplankton
    # From MesoZoo.F90 lines 247-259
    # Phytoplankton LFG: Food availability of prey Phytoplankton for predator Z3 and Z4
    available_phyto_c3 = (zoo_availability_parameters["del_carnivMesozoo_LO_P1"]*diatoms_LO_C) + (zoo_availability_parameters["del_carnivMesozoo_LO_P2"]*nanoflagellates_LO_C) + (zoo_availability_parameters["del_carnivMesozoo_LO_P3"]*picophyto_LO_C) + (zoo_availability_parameters["del_carnivMesozoo_LO_P4"]*largephyto_LO_C)
    available_phyto_c4 = (zoo_availability_parameters["del_omnivMesozoo_LO_P1"]*diatoms_LO_C) + (zoo_availability_parameters["del_omnivMesozoo_LO_P2"]*nanoflagellates_LO_C) + (zoo_availability_parameters["del_omnivMesozoo_LO_P3"]*picophyto_LO_C) + (zoo_availability_parameters["del_omnivMesozoo_LO_P4"]*largephyto_LO_C)
    
    # Mesozooplankton LFG
    available_mesozoo_c3 = (zoo_availability_parameters["del_z3z3"]*carnivMesozoo_LO_C) + (zoo_availability_parameters["del_z3z4"]*omnivMesozoo_LO_C)
    available_mesozoo_c4 = (zoo_availability_parameters["del_z4z3"]*carnivMesozoo_LO_C) + (zoo_availability_parameters["del_z4z4"]*omnivMesozoo_LO_C)
    
    # Microzooplankton LFG
    available_microzoo_c3 = (zoo_availability_parameters["del_z3z5"]*microzoo_LO_C) + (zoo_availability_parameters["del_z3z6"]*heteroFlagellates_LO_C)
    available_microzoo_c4 = (zoo_availability_parameters["del_z4z5"]*microzoo_LO_C) + (zoo_availability_parameters["del_z4z6"]*heteroFlagellates_LO_C)
#    sys.exit(available_microzoo_c4)
    
    # Total potential food (from Meso.F90 'rumc')
    f_c3 = available_phyto_c3 + available_mesozoo_c3 + available_microzoo_c3
    f_c4 = available_phyto_c4 + available_mesozoo_c4 + available_microzoo_c4
    
    # Calculate total food uptake rate (from Meso.F90 'rugc')
    total_uptake_rate_z3 = fTZ*mesozoo3_parameters["r_Z0"]*(mesozoo3_parameters["nu_z"]*f_c3/((mesozoo3_parameters["nu_z"]*f_c3) + mesozoo3_parameters["r_Z0"]))*carnivMesozoo_LO_C
    total_uptake_rate_z4 = fTZ*mesozoo4_parameters["r_Z0"]*(mesozoo4_parameters["nu_z"]*f_c4/((mesozoo4_parameters["nu_z"]*f_c4) + mesozoo4_parameters["r_Z0"]))*omnivMesozoo_LO_C
    
    # Calculate specific uptake rate considering potentially available food (from Meso.F90 'sut')
    specific_uptake_rate_z3 = total_uptake_rate_z3/(constant_parameters["p_small"] + f_c3)
    specific_uptake_rate_z4 = total_uptake_rate_z4/(constant_parameters["p_small"] + f_c4)

    # Total Gross Uptakes from every LFG
    dcarnivMesozoo_LO_Cdt_prd = {
            "p1": specific_uptake_rate_z3*zoo_availability_parameters["del_carnivMesozoo_LO_P1"]*diatoms_LO_C,
            "p2": specific_uptake_rate_z3*zoo_availability_parameters["del_carnivMesozoo_LO_P2"]*nanoflagellates_LO_C,
            "p3": specific_uptake_rate_z3*zoo_availability_parameters["del_carnivMesozoo_LO_P3"]*picophyto_LO_C,
            "p4": specific_uptake_rate_z3*zoo_availability_parameters["del_carnivMesozoo_LO_P4"]*largephyto_LO_C,
            "z3": specific_uptake_rate_z3*zoo_availability_parameters["del_z3z3"]*carnivMesozoo_LO_C,
            "z4": specific_uptake_rate_z3*zoo_availability_parameters["del_z3z4"]*omnivMesozoo_LO_C,
            "z5": specific_uptake_rate_z3*zoo_availability_parameters["del_z3z5"]*microzoo_LO_C,
            "z6": specific_uptake_rate_z3*zoo_availability_parameters["del_z3z6"]*heteroFlagellates_LO_C
            }

    domnivMesozoo_LO_Cdt_prd = {
            "p1": specific_uptake_rate_z4*zoo_availability_parameters["del_omnivMesozoo_LO_P1"]*diatoms_LO_C,
            "p2": specific_uptake_rate_z4*zoo_availability_parameters["del_omnivMesozoo_LO_P2"]*nanoflagellates_LO_C,
            "p3": specific_uptake_rate_z4*zoo_availability_parameters["del_omnivMesozoo_LO_P3"]*picophyto_LO_C,
            "p4": specific_uptake_rate_z4*zoo_availability_parameters["del_omnivMesozoo_LO_P4"]*largephyto_LO_C,
            "z3": specific_uptake_rate_z4*zoo_availability_parameters["del_z4z3"]*carnivMesozoo_LO_C,
            "z4": specific_uptake_rate_z4*zoo_availability_parameters["del_z4z4"]*omnivMesozoo_LO_C,
            "z5": specific_uptake_rate_z4*zoo_availability_parameters["del_z4z5"]*microzoo_LO_C,
            "z6": specific_uptake_rate_z4*zoo_availability_parameters["del_z4z6"]*heteroFlagellates_LO_C
            }
    # Total ingestion rate
    ic3 = 0.0
    in3 = 0.0
    ip3 = 0.0
    
    for key in dcarnivMesozoo_LO_Cdt_prd:
        ic3 += dcarnivMesozoo_LO_Cdt_prd[key]
        in3 += dcarnivMesozoo_LO_Cdt_prd[key]*conc_ratio_n[key]
        ip3 += dcarnivMesozoo_LO_Cdt_prd[key]*conc_ratio_p[key]
    
    ic4 = 0.0
    in4 = 0.0
    ip4 = 0.0
    
    for key in domnivMesozoo_LO_Cdt_prd:
        ic4 += domnivMesozoo_LO_Cdt_prd[key]
        in4 += domnivMesozoo_LO_Cdt_prd[key]*conc_ratio_n[key]
        ip4 += domnivMesozoo_LO_Cdt_prd[key]*conc_ratio_p[key]
  
    return dcarnivMesozoo_LO_Cdt_prd, domnivMesozoo_LO_Cdt_prd, ic3, in3, ip3, ic4, in4, ip4


def get_microzoo_predation_terms(conc, microzoo5_parameters, microzoo6_parameters, zoo_availability_parameters, environmental_parameters, constant_parameters, temp):
    """ Calculates the predation terms for microzooplankton """
    
    # Species concentrations
    pelBacteria_LO_C = conc[7]              # Pelagic bacteria carbon (mg C m^-3)
    pelBacteria_LO_N = conc[8]              # Pelagic bacteria nitrogen (mmol N m^-3)
    pelBacteria_LO_P = conc[9]              # Pelagic bacteria phosphate (mmol P m^-3)
    diatoms_LO_C = conc[10]             # Diatoms carbon (mg C m^-3)
    diatoms_LO_N = conc[11]             # Diatoms nitrogen (mmol N m^-3)
    diatoms_LO_P = conc[12]             # Diatoms phosphate (mmol P m^-3)
    diatoms_LO_Chl = conc[13]             # Diatoms chlorophyll (mg Chl-a m^-3)
    diatoms_LO_Si = conc[14]             # Diatoms silicate (mmol Si m^-3) 
    nanoflagellates_LO_C = conc[15]             # NanoFlagellates carbon (mg C m^-3)
    nanoflagellates_LO_N = conc[16]             # NanoFlagellates nitrogen (mmol N m^-3)
    nanoflagellates_LO_P = conc[17]             # NanoFlagellates phosphate (mmol P m^-3)
    nanoflagellates_LO_Chl = conc[18]             # NanoFlagellates chlorophyll (mg Chl-a m^-3)
    picophyto_LO_C = conc[19]             # Picophytoplankton carbon (mg C m^-3)
    picophyto_LO_N = conc[20]             # Picophytoplankton nitrogen (mmol N m^-3)
    picophyto_LO_P = conc[21]             # Picophytoplankton phosphate (mmol P m^-3)
    picophyto_LO_Chl = conc[22]             # Picophytoplankton chlorophyll (mg Chl-a m^-3)
    largephyto_LO_C = conc[23]             # Large phytoplankton carbon (mg C m^-3)
    largephyto_LO_N = conc[24]             # Large phytoplankton nitrogen (mmol N m^-3)
    largephyto_LO_P = conc[25]             # Large phytoplankton phosphate (mmol P m^-3) 
    largephyto_LO_Chl = conc[26]             # Large phytoplankton chlorophyll (mg Chl-a m^-3)
    carnivMesozoo_LO_C = conc[27]             # Carnivorous mesozooplankton carbon (mg C m^-3)
    carnivMesozoo_LO_N = conc[28]             # Carnivorous mesozooplankton nitrogen (mmol N m^-3)
    carnivMesozoo_LO_P = conc[29]             # Carnivorous mesozooplankton phosphate (mmol P m^-3)
    omnivMesozoo_LO_C = conc[30]             # Omnivorous mesozooplankton carbon (mg C m^-3)
    omnivMesozoo_LO_N = conc[31]             # Omnivorous mesozooplankton nitrogen (mmol N m^-3)
    omnivMesozoo_LO_P = conc[32]             # Omnivorous mesozooplankton phosphate (mmol P m^-3)
    microzoo_LO_C = conc[33]             # Microzooplankton carbon (mg C m^-3)
    microzoo_LO_N = conc[34]             # Microzooplankton nitrogen (mmol N m^-3)
    microzoo_LO_P = conc[35]             # Microzooplankton phosphate (mmol P m^-3)
    heteroFlagellates_LO_C = conc[36]             # Heterotrophic flagellates carbon (mg C m^-3)
    heteroFlagellates_LO_N = conc[37]             # Heterotrophic flagellates nitrogen (mmol N m^-3)
    heteroFlagellates_LO_P = conc[38]             # Heterotrophic flagellates phosphate (mmol P m^-3)
    
    # concentration ratios
    conc_ratio_n = {
            "b1": get_concentration_ratio(pelBacteria_LO_N, pelBacteria_LO_C, constant_parameters["p_small"]),
            "p1": get_concentration_ratio(diatoms_LO_N, diatoms_LO_C, constant_parameters["p_small"]),
            "p2": get_concentration_ratio(nanoflagellates_LO_N, nanoflagellates_LO_C, constant_parameters["p_small"]),
            "p3": get_concentration_ratio(picophyto_LO_N, picophyto_LO_C, constant_parameters["p_small"]),
            "p4": get_concentration_ratio(largephyto_LO_N, largephyto_LO_C, constant_parameters["p_small"]),
            "z5": get_concentration_ratio(microzoo_LO_N, microzoo_LO_C, constant_parameters["p_small"]),
            "z6": get_concentration_ratio(heteroFlagellates_LO_N, heteroFlagellates_LO_C, constant_parameters["p_small"])
            }
    
    conc_ratio_p = {
            "b1": get_concentration_ratio(pelBacteria_LO_P, pelBacteria_LO_C, constant_parameters["p_small"]),
            "p1": get_concentration_ratio(diatoms_LO_P, diatoms_LO_C, constant_parameters["p_small"]),
            "p2": get_concentration_ratio(nanoflagellates_LO_P, nanoflagellates_LO_C, constant_parameters["p_small"]),
            "p3": get_concentration_ratio(picophyto_LO_P, picophyto_LO_C, constant_parameters["p_small"]),
            "p4": get_concentration_ratio(largephyto_LO_P, largephyto_LO_C, constant_parameters["p_small"]),
            "z5": get_concentration_ratio(microzoo_LO_P, microzoo_LO_C, constant_parameters["p_small"]),
            "z6": get_concentration_ratio(heteroFlagellates_LO_P, heteroFlagellates_LO_C, constant_parameters["p_small"])
            }
    
    # Zooplankton temperature regulating factor
    fTZ = eTq_vector(temp, environmental_parameters["basetemp"], environmental_parameters["q10z"])
    
    # Capture efficiencies 
    capture_efficiencies_z5 = {
            "b1": pelBacteria_LO_C/(pelBacteria_LO_C + microzoo5_parameters["mu_z"]),
            "p1": diatoms_LO_C/(diatoms_LO_C + microzoo5_parameters["mu_z"]),
            "p2": nanoflagellates_LO_C/(nanoflagellates_LO_C + microzoo5_parameters["mu_z"]),
            "p3": picophyto_LO_C/(picophyto_LO_C + microzoo5_parameters["mu_z"]),
            "p4": largephyto_LO_C/(largephyto_LO_C + microzoo5_parameters["mu_z"]),
            "z5": microzoo_LO_C/(microzoo_LO_C + microzoo5_parameters["mu_z"]),
            "z6": heteroFlagellates_LO_C/(heteroFlagellates_LO_C + microzoo5_parameters["mu_z"])
            }
    
    capture_efficiencies_z6 = {
            "b1": pelBacteria_LO_C/(pelBacteria_LO_C + microzoo6_parameters["mu_z"]),
            "p1": diatoms_LO_C/(diatoms_LO_C + microzoo6_parameters["mu_z"]),
            "p2": nanoflagellates_LO_C/(nanoflagellates_LO_C + microzoo6_parameters["mu_z"]),
            "p3": picophyto_LO_C/(picophyto_LO_C + microzoo6_parameters["mu_z"]),
            "p4": largephyto_LO_C/(largephyto_LO_C + microzoo6_parameters["mu_z"]),
            "z5": microzoo_LO_C/(microzoo_LO_C + microzoo6_parameters["mu_z"]),
            "z6": heteroFlagellates_LO_C/(heteroFlagellates_LO_C + microzoo6_parameters["mu_z"])
            }
    
    # Calculate total potential food given the non-dim prey availability and capture efficiency
    # From MicroZoo.F90 lines 209-237
    # Bacteria LFG: Food availability of prey Bacteria for predator Z5 and Z6
    available_bact_c5 = zoo_availability_parameters["del_z5b1"]*pelBacteria_LO_C*capture_efficiencies_z5["b1"]
    available_bact_n5 = zoo_availability_parameters["del_z5b1"]*pelBacteria_LO_C*capture_efficiencies_z5["b1"]*conc_ratio_n["b1"]
    available_bact_p5 = zoo_availability_parameters["del_z5b1"]*pelBacteria_LO_C*capture_efficiencies_z5["b1"]*conc_ratio_p["b1"]
    available_bact_c6 = zoo_availability_parameters["del_z6b1"]*pelBacteria_LO_C*capture_efficiencies_z6["b1"]
    available_bact_n6 = zoo_availability_parameters["del_z6b1"]*pelBacteria_LO_C*capture_efficiencies_z6["b1"]*conc_ratio_n["b1"]
    available_bact_p6 = zoo_availability_parameters["del_z6b1"]*pelBacteria_LO_C*capture_efficiencies_z6["b1"]*conc_ratio_p["b1"]
    
    # Phytoplankton LFG: Food availability of prey Phytoplankton for predator Z5 and Z6
    available_phyto_c5 = ((zoo_availability_parameters["del_microzoo_LO_P1"]*diatoms_LO_C*capture_efficiencies_z5["p1"]) + (zoo_availability_parameters["del_carnivMesozoo_LO_P2"]*nanoflagellates_LO_C*capture_efficiencies_z5["p2"]) + 
                          (zoo_availability_parameters["del_carnivMesozoo_LO_P3"]*picophyto_LO_C*capture_efficiencies_z5["p3"]) + (zoo_availability_parameters["del_carnivMesozoo_LO_P4"]*largephyto_LO_C*capture_efficiencies_z5["p4"]))
    available_phyto_n5 = ((zoo_availability_parameters["del_microzoo_LO_P1"]*diatoms_LO_C*capture_efficiencies_z5["p1"]*conc_ratio_n["p1"]) + (zoo_availability_parameters["del_carnivMesozoo_LO_P2"]*nanoflagellates_LO_C*capture_efficiencies_z5["p2"]*conc_ratio_n["p2"]) + 
                          (zoo_availability_parameters["del_carnivMesozoo_LO_P3"]*picophyto_LO_C*capture_efficiencies_z5["p3"]*conc_ratio_n["p3"]) + (zoo_availability_parameters["del_carnivMesozoo_LO_P4"]*largephyto_LO_C*capture_efficiencies_z5["p4"]*conc_ratio_n["p4"]))
    available_phyto_p5 = ((zoo_availability_parameters["del_microzoo_LO_P1"]*diatoms_LO_C*capture_efficiencies_z5["p1"]*conc_ratio_p["p1"]) + (zoo_availability_parameters["del_carnivMesozoo_LO_P2"]*nanoflagellates_LO_C*capture_efficiencies_z5["p2"]*conc_ratio_p["p2"]) + 
                          (zoo_availability_parameters["del_carnivMesozoo_LO_P3"]*picophyto_LO_C*capture_efficiencies_z5["p3"]*conc_ratio_p["p3"]) + (zoo_availability_parameters["del_carnivMesozoo_LO_P4"]*largephyto_LO_C*capture_efficiencies_z5["p4"]*conc_ratio_p["p4"]))
    available_phyto_c6 = ((zoo_availability_parameters["del_heteroFlagellates_LO_P1"]*diatoms_LO_C*capture_efficiencies_z6["p1"]) + (zoo_availability_parameters["del_heteroFlagellates_LO_P2"]*nanoflagellates_LO_C*capture_efficiencies_z6["p2"]) + 
                          (zoo_availability_parameters["del_heteroFlagellates_LO_P3"]*picophyto_LO_C*capture_efficiencies_z6["p3"]) + (zoo_availability_parameters["del_heteroFlagellates_LO_P4"]*largephyto_LO_C*capture_efficiencies_z6["p4"]))
    available_phyto_n6 = ((zoo_availability_parameters["del_heteroFlagellates_LO_P1"]*diatoms_LO_C*capture_efficiencies_z6["p1"]*conc_ratio_n["p1"]) + (zoo_availability_parameters["del_heteroFlagellates_LO_P2"]*nanoflagellates_LO_C*capture_efficiencies_z6["p2"]*conc_ratio_n["p2"]) + 
                          (zoo_availability_parameters["del_heteroFlagellates_LO_P3"]*picophyto_LO_C*capture_efficiencies_z6["p3"]*conc_ratio_n["p3"]) + (zoo_availability_parameters["del_heteroFlagellates_LO_P4"]*largephyto_LO_C*capture_efficiencies_z6["p4"]*conc_ratio_n["p4"]))
    available_phyto_p6 = ((zoo_availability_parameters["del_heteroFlagellates_LO_P1"]*diatoms_LO_C*capture_efficiencies_z6["p1"]*conc_ratio_p["p1"]) + (zoo_availability_parameters["del_heteroFlagellates_LO_P2"]*nanoflagellates_LO_C*capture_efficiencies_z6["p2"]*conc_ratio_p["p2"]) + 
                          (zoo_availability_parameters["del_heteroFlagellates_LO_P3"]*picophyto_LO_C*capture_efficiencies_z6["p3"]*conc_ratio_p["p3"]) + (zoo_availability_parameters["del_heteroFlagellates_LO_P4"]*largephyto_LO_C*capture_efficiencies_z6["p4"]*conc_ratio_p["p4"]))
    
    # Phytoplankton LFG: Food availability of prey Microzooplankton for predator Z5 and Z6
    available_microzoo_c5 = (zoo_availability_parameters["del_z5z5"]*microzoo_LO_C*capture_efficiencies_z5["z5"]) + (zoo_availability_parameters["del_z5z6"]*heteroFlagellates_LO_C*capture_efficiencies_z5["z6"])
    available_microzoo_n5 = (zoo_availability_parameters["del_z5z5"]*microzoo_LO_C*capture_efficiencies_z5["z5"]*conc_ratio_n["z5"]) + (zoo_availability_parameters["del_z5z6"]*heteroFlagellates_LO_C*capture_efficiencies_z5["z6"]*conc_ratio_n["z6"])
    available_microzoo_p5 = (zoo_availability_parameters["del_z5z5"]*microzoo_LO_C*capture_efficiencies_z5["z5"]*conc_ratio_p["z5"]) + (zoo_availability_parameters["del_z5z6"]*heteroFlagellates_LO_C*capture_efficiencies_z5["z6"]*conc_ratio_p["z6"])
    
    available_microzoo_c6 = (zoo_availability_parameters["del_z6z5"]*microzoo_LO_C*capture_efficiencies_z6["z5"]) + (zoo_availability_parameters["del_z6z6"]*heteroFlagellates_LO_C*capture_efficiencies_z6["z6"])
    available_microzoo_n6 = (zoo_availability_parameters["del_z6z5"]*microzoo_LO_C*capture_efficiencies_z6["z5"]**conc_ratio_n["z5"]) + (zoo_availability_parameters["del_z6z6"]*heteroFlagellates_LO_C*capture_efficiencies_z6["z6"]**conc_ratio_n["z6"])
    available_microzoo_p6 = (zoo_availability_parameters["del_z6z5"]*microzoo_LO_C*capture_efficiencies_z6["z5"]**conc_ratio_p["z5"]) + (zoo_availability_parameters["del_z6z6"]*heteroFlagellates_LO_C*capture_efficiencies_z6["z6"]**conc_ratio_p["z6"])
    
    # Total potential food (from MicroZoo.F90 'rumc', 'rumn', 'rump')
    f_c5 = available_bact_c5 + available_phyto_c5 + available_microzoo_c5
    f_n5 = available_bact_n5 + available_phyto_n5 + available_microzoo_n5
    f_p5 = available_bact_p5 + available_phyto_p5 + available_microzoo_p5
    f_c6 = available_bact_c6 + available_phyto_c6 + available_microzoo_c6
    f_n6 = available_bact_n6 + available_phyto_n6 + available_microzoo_n6
    f_p6 = available_bact_p6 + available_phyto_p6 + available_microzoo_p6
    
    # Calculate total food uptake rate (from MicroZoo.F90 line 243 'rugc')
    total_uptake_rate_z5 = fTZ*microzoo5_parameters["r_Z0"]*(f_c5/(f_c5 + microzoo5_parameters["h_Z_F"]))*microzoo_LO_C
    total_uptake_rate_z6 = fTZ*microzoo6_parameters["r_Z0"]*(f_c6/(f_c6 + microzoo6_parameters["h_Z_F"]))*heteroFlagellates_LO_C
    
    # Calculate specific uptake rate considering potentially available food (from MicroZoo.F90 line 244 'sut')
    specific_uptake_rate_z5 = total_uptake_rate_z5/(f_c5 + constant_parameters["p_small"])
    specific_uptake_rate_z6 = total_uptake_rate_z6/(f_c6 + constant_parameters["p_small"])
    
    # Total Gross Uptakes from every LFG
    dmicrozoo_LO_Cdt_prd = {
            "b1": specific_uptake_rate_z5*zoo_availability_parameters["del_z5b1"]*pelBacteria_LO_C*capture_efficiencies_z5["b1"],
            "p1": specific_uptake_rate_z5*zoo_availability_parameters["del_microzoo_LO_P1"]*diatoms_LO_C*capture_efficiencies_z5["p1"],
            "p2": specific_uptake_rate_z5*zoo_availability_parameters["del_microzoo_LO_P2"]*nanoflagellates_LO_C*capture_efficiencies_z5["p2"],
            "p3": specific_uptake_rate_z5*zoo_availability_parameters["del_microzoo_LO_P3"]*picophyto_LO_C*capture_efficiencies_z5["p3"],
            "p4": specific_uptake_rate_z5*zoo_availability_parameters["del_microzoo_LO_P4"]*largephyto_LO_C*capture_efficiencies_z5["p4"],
            "z5": specific_uptake_rate_z5*zoo_availability_parameters["del_z5z5"]*microzoo_LO_C*capture_efficiencies_z5["z5"],
            "z6": specific_uptake_rate_z5*zoo_availability_parameters["del_z5z6"]*heteroFlagellates_LO_C*capture_efficiencies_z5["z6"]
            }
        
    dheteroFlagellates_LO_Cdt_prd = {
            "b1": specific_uptake_rate_z6*zoo_availability_parameters["del_z6b1"]*pelBacteria_LO_C*capture_efficiencies_z6["b1"],
            "p1": specific_uptake_rate_z6*zoo_availability_parameters["del_heteroFlagellates_LO_P1"]*diatoms_LO_C*capture_efficiencies_z6["p1"],
            "p2": specific_uptake_rate_z6*zoo_availability_parameters["del_heteroFlagellates_LO_P2"]*nanoflagellates_LO_C*capture_efficiencies_z6["p2"],
            "p3": specific_uptake_rate_z6*zoo_availability_parameters["del_heteroFlagellates_LO_P3"]*picophyto_LO_C*capture_efficiencies_z6["p3"],
            "p4": specific_uptake_rate_z6*zoo_availability_parameters["del_heteroFlagellates_LO_P4"]*largephyto_LO_C*capture_efficiencies_z6["p4"],
            "z5": specific_uptake_rate_z6*zoo_availability_parameters["del_z6z5"]*microzoo_LO_C*capture_efficiencies_z6["z5"],
            "z6": specific_uptake_rate_z6*zoo_availability_parameters["del_z6z6"]*heteroFlagellates_LO_C*capture_efficiencies_z6["z6"]
            }
    
    # Total ingestion rate
    ic5 = 0.0
    in5 = 0.0
    ip5 = 0.0
    
    for key in dmicrozoo_LO_Cdt_prd:
        ic5 += dmicrozoo_LO_Cdt_prd[key]
        in5 += dmicrozoo_LO_Cdt_prd[key]*conc_ratio_n[key]
        ip5 += dmicrozoo_LO_Cdt_prd[key]*conc_ratio_p[key]
    
    ic6 = 0.0
    in6 = 0.0
    ip6 = 0.0
    
    for key in dheteroFlagellates_LO_Cdt_prd:
        ic6 += dheteroFlagellates_LO_Cdt_prd[key]
        in6 += dheteroFlagellates_LO_Cdt_prd[key]*conc_ratio_n[key]
        ip6 += dheteroFlagellates_LO_Cdt_prd[key]*conc_ratio_p[key]
    
    return dmicrozoo_LO_Cdt_prd, dheteroFlagellates_LO_Cdt_prd, ic5, in5, ip5, ic6, in6, ip6
