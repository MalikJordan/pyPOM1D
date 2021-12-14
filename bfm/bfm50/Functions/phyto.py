import numpy
import sys
from Functions.other_functions import insw_vector, eTq_vector, get_concentration_ratio

def phyto_eqns(conc, phyto_parameters, env_parameters, constant_parameters, group, pc, pn, pp, pl, qs, temp, time):
    """ Calculates the terms needed for the phytoplnaktion biological rate equations
        Equations come from the BFM user manual
    """
    
    # Species concentrations
    phospate_IO_P = conc[1]               # Phosphate (mmol P m^-3)
    nitrate_IO_N = conc[2]               # Nitrate (mmol N m^-3)
    ammonium_IO_N = conc[3]               # Ammonium (mmol N m^-3)
    silicate_IO_Si = conc[5]               # Silicate (mmol Si m^-3)
    diatoms_LO_Chl = conc[13]              # Diatoms chlorophyll (mg Chl-a m^-3)
    diatoms_LO_Si = conc[14]              # Diatoms silicate (mmol Si m^-3) 
    nanoflagellates_LO_Chl = conc[18]              # NanoFlagellates chlorophyll (mg Chl-a m^-3)
    picophyto_LO_Chl = conc[22]              # Picophytoplankton chlorophyll (mg Chl-a m^-3)
    largephyto_LO_Chl = conc[26]              # Large phytoplankton chlorophyll (mg Chl-a m^-3)
    particOrganDetritus_NO_C = conc[44]              # Particulate organic carbon (mg C m^-3)
    
    # Concentration ratios (constituents quota in phytoplankton)
    pn_pc = get_concentration_ratio(pn, pc, constant_parameters["p_small"])
    pp_pc = get_concentration_ratio(pp, pc, constant_parameters["p_small"])
    pl_pc = get_concentration_ratio(pl, pc, constant_parameters["p_small"])

    #--------------------------------------------------------------------------
    # Temperature response of Phytoplankton Include cut-off at low temperature if p_temp>0
    et = eTq_vector(temp, env_parameters["basetemp"], env_parameters["q10z"])
    fTP = max(0.0, et - phyto_parameters["p_temp"])

    #--------------------------------------------------------------------------
    # Nutrient limitations (intracellular and extracellular) fpplim is the 
    # combined non-dimensional factor limiting photosynthesis
    # from Phyto.F90 lines 268-308
    iphospate_IO_P = min(1.0, max(constant_parameters["p_small"], (pp_pc - phyto_parameters["phi_Pmin"])/(phyto_parameters["phi_Popt"] - phyto_parameters["phi_Pmin"])))
    in1n = min(1.0, max(constant_parameters["p_small"], (pn_pc - phyto_parameters["phi_Nmin"])/(phyto_parameters["phi_Nopt"] - phyto_parameters["phi_Nmin"])))
    
    if group == 1:
        fpplim = min(1.0, silicate_IO_Si/(silicate_IO_Si + phyto_parameters["h_Ps"] + (phyto_parameters["rho_Ps"]*diatoms_LO_Si)))
    else:
        fpplim = 1.0

    #--------------------------------------------------------------------------
    # multiple nutrient limitation, Liebig rule (from Phyto.F90 line ~318, iN)
    multiple_nut_lim = min(iphospate_IO_P, in1n)
    
    #--------------------------------------------------------------------------    
    # Total extinction coef (m^-1)
    suspended_sediments = 0.0
    # from CalcVerticalExtinction.F90 line 82
    xEPS = env_parameters["p_eps0"] + env_parameters["p_epsESS"]*suspended_sediments + env_parameters["p_epsR6"]*particOrganDetritus_NO_C
    # from CalcVerticalExtinction.F90 line 101 (ChlAttenFlag=1, ChlDynamicsFlag=2)       
    xEPS = xEPS + diatoms_LO_Chl + nanoflagellates_LO_Chl*2 + picophyto_LO_Chl*3 + largephyto_LO_Chl*4
    
    #-------------------------------------------------------------------------- 
    # Light limitation with Chl dynamics
    # irradiance (uE m^-2 s^-1) from Phyto.F90 lines 353-355
    irradiance = qs*env_parameters["epsilon_PAR"]/constant_parameters["e2w"]
    r = xEPS * env_parameters["del_z"]
    r = irradiance/xEPS/env_parameters["del_z"]*(1.0 - numpy.exp(-r))
    irr = max(constant_parameters["p_small"], r)
    
    # Compute exponent E_PAR/E_K = alpha0/PBmax (part of eqn. 2.2.4)
    exponent = pl_pc*phyto_parameters["alpha_chl"]/phyto_parameters["rP0"]*irr

    # light limitation factor (from Phyto.f90 line 374, eiPPY)
    light_lim = (1.0 - numpy.exp(-exponent))

    #--------------------------------------------------------------------------
    # total photosynthesis (from Phyto.F90 line ~380, sum)
    photosynthesis = phyto_parameters["rP0"]*fTP*light_lim*fpplim

    #--------------------------------------------------------------------------
    # Lysis nad excretion
    # nutr. -stress lysis (from Phyto.F90 lines ~385-387, sdo)
    nut_stress_lysis = (phyto_parameters["h_Pnp"]/(multiple_nut_lim + phyto_parameters["h_Pnp"]))*phyto_parameters["d_P0"]
    nut_stress_lysis += phyto_parameters["p_seo"]*pc/(pc + phyto_parameters["p_sheo"] + constant_parameters["p_small"])

    # activity excretion (Phyto.F90 line 389)
    activity_excretion = photosynthesis*phyto_parameters["betaP"]

    # nutrient stress excretion from Phyto.F90 line 396
    nut_stress_excretion = photosynthesis*(1.0 - phyto_parameters["betaP"])*(1.0 - multiple_nut_lim)

    #--------------------------------------------------------------------------
    # Apportioning over R1 and R6: Cell lysis generates both DOM and POM
    pe_R6 = min(phyto_parameters["phi_Pmin"]/(pp_pc + constant_parameters["p_small"]), phyto_parameters["phi_Nmin"]/(pn_pc + constant_parameters["p_small"]))
    pe_R6 = min(1.0, pe_R6)
    rparticOrganDetritus_NO_C = pe_R6*nut_stress_lysis*pc
    rlabileDOM_NO_C = (1.0 - pe_R6)*nut_stress_lysis*pc

    #--------------------------------------------------------------------------
    # Respiration rate
    # activity (from Phyto.F90 line 416)
    activity_rsp = phyto_parameters["gammaP"]*(photosynthesis - activity_excretion - nut_stress_excretion)

    # basal (from Phyto.F90 line 417)
    basal_rsp = fTP*phyto_parameters["bP"]

    # total (from Phyto.F90 line 418)
    total_rsp = activity_rsp + basal_rsp

    # total actual respiration
    dPcdt_rsp_disInorgCarbon_IO_C = total_rsp*pc

    #--------------------------------------------------------------------------
    # Production, productivity and C flows
    # Phytoplankton gross primary production [mg C m^-3 s^-1]
    dPcdt_gpp_disInorgCarbon_IO_C = photosynthesis*pc

    # specific loss terms (from Phyto.F90 line 428)
    specific_loss_terms = activity_excretion + nut_stress_excretion + total_rsp + nut_stress_lysis

    # All activity excretions are assigned to R1
    # p_switchDOC=1 and P_netgrowth=FLASE: [mg C m^-3 s^-1]
    rlabileDOM_NO_C += activity_excretion*pc + nut_stress_excretion*pc
    dPcdt_exu_semilabileDOC_NO_C = 0.0

    # Phytoplankton DOM cell lysis- carbon lost to DOM [mg C m^-3 s^-1]
    dPcdt_lys_labileDOM_NO_C = rlabileDOM_NO_C

    # Phytoplankton POM cell lysis- carbon lost to POM (eqn. 2.2.9) [mg C m^-3 s^-1]
    dPcdt_lys_particOrganDetritus_NO_C = rparticOrganDetritus_NO_C

    #--------------------------------------------------------------------------
    # Potential-Net primary production
    # from Phyto.F90 line 455
    sadap = fTP*phyto_parameters["rP0"]
    
    # Net production (from Phyto.F90 line 457, 'run')
    net_production = max(0.0, (photosynthesis - specific_loss_terms)*pc)
    
    #--------------------------------------------------------------------------
    # Nutrient Uptake: calculate maximum uptake of N, P based on affinity

    cqun3 = phyto_parameters["h_Pn"]/(constant_parameters["p_small"] + phyto_parameters["h_Pn"] + ammonium_IO_N)

    # max potential uptake of N3 (from Phyto.F90 'rumn3')
    max_upt_nitrate_IO_N = phyto_parameters["a_N"]*nitrate_IO_N*pc*cqun3

    # max potential uptake of N4 (from Phyto.F90 'rumn4')
    max_upt_ammonium_IO_N = phyto_parameters["a_N"]*ammonium_IO_N*pc

    # max potential uptake of DIN (from Phyto.F90 'rumn')
    max_upt_DIN = max_upt_nitrate_IO_N + max_upt_ammonium_IO_N

    # max pot. uptake of PO4 (from Phyto.F90 line 468)
    rump = phyto_parameters["a_P"]*phospate_IO_P*pc

    #--------------------------------------------------------------------------
    # Nutrient dynamics: NITROGEN

    # Intracellular missing amount of N (from Phyto.F90)
    misn = sadap*(phyto_parameters["phi_Nmax"]*pc - pn)

    # N uptake based on net assimilat. C (from Phyto.F90)
    rupn = phyto_parameters["phi_Nmax"]*net_production

    # actual uptake of NI (from Phyto.F90, 'runn')
    dPndt_upt = min(max_upt_DIN, rupn + misn)

    # if nitrogen uptake rate is positive, then uptake is divided between coming from the nitrate and ammonium reservoir
    # if nitrogen uptake is negative, all nitrogen goes to the DOM pool
    upt_switch_n = insw_vector(dPndt_upt)

    # actual uptake of nitrate_IO_N (from Phyto.F90, 'runn3')
    dPndt_upt_nitrate_IO_N = upt_switch_n*dPndt_upt*max_upt_nitrate_IO_N/(constant_parameters["p_small"] + max_upt_DIN)

    # actual uptake of ammonium_IO_N (from Phyto.F90, 'runn4')
    dPndt_upt_ammonium_IO_N = upt_switch_n*dPndt_upt*max_upt_ammonium_IO_N/(constant_parameters["p_small"] + max_upt_DIN)

    extra_n = -dPndt_upt*(1.0 - upt_switch_n)

    #--------------------------------------------------------------------------
    # Nutrient dynamics: PHOSPHORUS

    # intracellular missing amount of P (from Phyto.F90 line 514)
    misp = sadap*(phyto_parameters["phi_Pmax"]*pc-pp)

    # P uptake based on C uptake (from Phyto.F90 line 517)
    rupp = phyto_parameters["phi_Pmax"]*net_production

    # Actual uptake
    runp = min(rump, rupp + misp)
    upt_switch_p = insw_vector(runp)
    dPpdt_upt_phospate_IO_P = runp*upt_switch_p

    # is uptake is negative flux goes to DIP (labileDOM_NO_P) pool
    dPpdt_upt_labileDOM_NO_P = -runp*(1.0 - upt_switch_p)

    #--------------------------------------------------------------------------
    # Excretion of N and P to PON and POP
    dPndt_lys_particOrganDetritus_NO_N = pe_R6*nut_stress_lysis*pn
    dPndt_lys_labileDOM_NO_N = nut_stress_lysis*pn - dPndt_lys_particOrganDetritus_NO_N

    dPpdt_lys_particOrganDetritus_NO_P = pe_R6*nut_stress_lysis*pp
    dPpdt_lys_labileDOM_NO_P = nut_stress_lysis*pp - dPpdt_lys_particOrganDetritus_NO_P

    #--------------------------------------------------------------------------
    # Nutrient dynamics: SILICATE
    if group == 1:
        # Gross uptake of silicate excluding respiratory costs (from Phyto.F90, 'runs')
        dPsdt_upt_silicate_IO_Si = max(0.0, phyto_parameters["phi_Sopt"]*pc*(photosynthesis - basal_rsp))

        # losses of Si (from Phyto.F90)
        dPsdt_lys_particOrganDetritus_NO_Si = nut_stress_lysis*diatoms_LO_Si
    else:
        dPsdt_upt_silicate_IO_Si = 0.0
        dPsdt_lys_particOrganDetritus_NO_Si = 0.0

    #--------------------------------------------------------------------------
    # Chl-a synthesis and photoacclimation
    if phyto_parameters["chl_switch"] == 1:
        # dynamical chl:c ratio from Fortran code Phyto.F90
        rho_chl = phyto_parameters["theta_chl0"]*min(1.0, phyto_parameters["rP0"]*light_lim*pc/(phyto_parameters["alpha_chl"]*(pl + constant_parameters["p_small"])*irr))
        
        # Chlorophyll synthesis from Fortran code Phyto.F90, rate_chl
        dPldt_syn = rho_chl*(photosynthesis - nut_stress_excretion - activity_excretion - activity_rsp)*pc - nut_stress_lysis*pl
    else:
        sys.exit("Warning: This code does not support other chl systhesis parameterizations")
    
    #--------------------------------------------------------------------------

    return (dPcdt_gpp_disInorgCarbon_IO_C, dPcdt_rsp_disInorgCarbon_IO_C, dPcdt_lys_labileDOM_NO_C, dPcdt_lys_particOrganDetritus_NO_C, dPcdt_exu_semilabileDOC_NO_C, 
            dPndt_upt_nitrate_IO_N, dPndt_upt_ammonium_IO_N, extra_n, dPndt_lys_labileDOM_NO_N, dPndt_lys_particOrganDetritus_NO_N, 
            dPpdt_upt_phospate_IO_P, dPpdt_upt_labileDOM_NO_P, dPpdt_lys_labileDOM_NO_P, dPpdt_lys_particOrganDetritus_NO_P, 
            dPldt_syn, dPsdt_upt_silicate_IO_Si, dPsdt_lys_particOrganDetritus_NO_Si)
