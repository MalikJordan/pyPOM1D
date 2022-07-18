from __future__ import division, print_function
import numpy as np
import json
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from bfm.bfm50.Functions.seasonal_cycling_functions import get_wind, get_salinity, get_sunlight, get_temperature, calculate_density
from bfm.bfm50.Functions.phyto import phyto_eqns
from bfm.bfm50.Functions.bacteria import bacteria_eqns
from bfm.bfm50.Functions.predation import get_mesozoo_predation_terms, get_microzoo_predation_terms
from bfm.bfm50.Functions.micro import microzoo_eqns
from bfm.bfm50.Functions.meso import mesozoo_eqns
from bfm.bfm50.Functions.pel_chem import pel_chem_eqns
from bfm.bfm50.Functions.oxygen import calculate_oxygen_reaeration
from bfm.bfm50.Functions.co2_flux_functions import calculate_co2_flux
from bfm.bfm50.Functions.other_functions import insw_vector, get_concentration_ratio
from pom.constants import current_path


# Names of species in the system
species_names = ['disOxygen_IO_O', 'phospate_IO_P', 'nitrate_IO_N', 'ammonium_IO_N', 'nitrogenSink', 'silicate_IO_Si', 'reductEquiv_IO_R', 'pelBacteria_LO_C', 'pelBacteria_LO_N', 'pelBacteria_LO_P', 
                 'diatoms_LO_C', 'diatoms_LO_N', 'diatoms_LO_P', 'diatoms_LO_Chl', 'diatoms_LO_Si', 'nanoflagellates_LO_C', 'nanoflagellates_LO_N', 'nanoflagellates_LO_P', 'nanoflagellates_LO_Chl',
                 'picophyto_LO_C', 'picophyto_LO_N', 'picophyto_LO_P', 'picophyto_LO_Chl', 'largephyto_LO_C', 'largephyto_LO_N', 'largephyto_LO_P', 'largephyto_LO_Chl',
                 'carnivMesozoo_LO_C', 'carnivMesozoo_LO_N', 'carnivMesozoo_LO_P', 'omnivMesozoo_LO_C', 'omnivMesozoo_LO_N', 'omnivMesozoo_LO_P', 'microzoo_LO_C', 'microzoo_LO_N', 'microzoo_LO_P',
                 'heteroFlagellates_LO_C', 'heteroFlagellates_LO_N', 'heteroFlagellates_LO_P', 'labileDOM_NO_C', 'labileDOM_NO_N', 'labileDOM_NO_P', 'semilabileDOC_NO_C', 'semirefractDOC_NO_C', 'particOrganDetritus_NO_C', 
                 'particOrganDetritus_NO_N', 'particOrganDetritus_NO_P', 'particOrganDetritus_NO_Si', 'disInorgCarbon_IO_C', 'totalAlkalinity_IO']


def bfm50_rate_eqns(count, bfm_phys_vars, time, conc, seasonal_cycle=False):
    """ Calculates the change in concentration for the 50 state variables
        NOTE: iron dynamics are not included, this is to compare to the standalone pelagic system
    """



    #--------------------------------------------------------------------------
    # # import parameters from json file
    # with open("bfm50_parameters.json", "r") as read_parameters:
    with open(current_path + "/bfm/bfm50/bfm50_parameters.json") as read_parameters:
        parameters = json.load(read_parameters)
    
    constant_parameters = parameters["constants"]
    environmental_parameters = parameters["environmental_parameters"]
    bacteria_parameters = parameters["bacteria_parameters"]
    phyto1_prameters = parameters["phyto1_parameters"]
    phyto2_prameters = parameters["phyto2_parameters"]
    phyto3_prameters = parameters["phyto3_parameters"]
    phyto4_prameters = parameters["phyto4_parameters"]
    mesozoo3_parameters = parameters["mesozoo3_parameters"]
    mesozoo4_parameters = parameters["mesozoo4_parameters"]
    microzoo5_parameters = parameters["microzoo5_parameters"]
    microzoo6_parameters = parameters["microzoo6_parameters"]
    zoo_availability_parameters = parameters["zoo_availability_parameters"]
    pel_chem_parameters = parameters["pel_chem_parameters"]
    co2_flux_parameters = parameters["co2_flux_parameters"]
    oxygen_reaeration_parameters = parameters["oxygen_reaeration_parameters"]
    
    #--------------------------------------------------------------------------
    # Seasonal wind, temp, salinity, and radiation values
    if seasonal_cycle:

        t = time

        # Wind
        w_win = environmental_parameters["w_win"]                               # Winter wind speed
        w_sum = environmental_parameters["w_sum"]                               # Summer wind speed
        wind = get_wind(t,w_win,w_sum)                                          # Yearly wind cylce

        # Temperature
        t_win = environmental_parameters["t_win"]                               # Winter temp value
        t_sum = environmental_parameters["t_sum"]                               # Summer temp value
        tde = environmental_parameters["tde"]                                   # Sinusoidal temperature daily excursion degC
        temper = get_temperature(t,t_win,t_sum, tde)                            # Yearly temp cycle

        # Salinity
        s_win = environmental_parameters["s_win"]                               # Winter salinity value
        s_sum = environmental_parameters["s_sum"]                               # Summer salinity value
        salt = get_salinity(t,s_win,s_sum)                                      # Yearly salinity cycle

        # Short wave irradiance flux (W/m^2)
        qs_win = environmental_parameters["qs_win"]                             # Winter irradiance value
        qs_sum = environmental_parameters["qs_sum"]                             # Summer irradiance value
        qs = get_sunlight(t,qs_win,qs_sum)                                      # Yearly irradiance cycle

    else:
        # wind = environmental_parameters["w_win"]
        # temper = environmental_parameters["t_win"]
        # salt = environmental_parameters["s_win"]
        # qs = environmental_parameters["qs_win"]
        wind = bfm_phys_vars.wind
        temper = bfm_phys_vars.temperature[count]
        salt = bfm_phys_vars.salinity[count]
        # qs = bfm_phys_vars.irradiance
        xEPS = bfm_phys_vars.vertical_extinction[count,:]
        irradiance = bfm_phys_vars.irradiance[count,:]
        del_z = bfm_phys_vars.depth[count]
        suspended_sediments = bfm_phys_vars.suspended_matter[count]
        rho = bfm_phys_vars.density[count]


    #--------------------------------------------------------------------------
    # State variables
    disOxygen_IO_O = conc[0]              # Dissolved oxygen (mg O_2 m^-3)
    phospate_IO_P = conc[1]              # Phosphate (mmol P m^-3)
    nitrate_IO_N = conc[2]              # Nitrate (mmol N m^-3)
    ammonium_IO_N = conc[3]              # Ammonium (mmol N m^-3)
    nitrogenSink = conc[4]              # Nitrogen sink (mmol N m^-3)
    silicate_IO_Si = conc[5]              # Silicate (mmol Si m^-3)
    reductEquiv_IO_R = conc[6]              # Reduction equivalents (mmol S m^-3)
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
    labileDOM_NO_C = conc[39]             # Labile dissolved organic carbon (mg C m^-3)
    labileDOM_NO_N = conc[40]             # Labile dissolved organic nitrogen (mmol N m^-3)
    labileDOM_NO_P = conc[41]             # Labile dissolved organic phosphate (mmol P m^-3)
    semilabileDOC_NO_C = conc[42]             # Semi-labile dissolved organic carbon (mg C m^-3)
    semirefractDOC_NO_C = conc[43]             # Semi-refractory Dissolved Organic Carbon (mg C m^-3)
    particOrganDetritus_NO_C = conc[44]             # Particulate organic carbon (mg C m^-3)
    particOrganDetritus_NO_N = conc[45]             # Particulate organic nitrogen (mmol N m^-3)
    particOrganDetritus_NO_P = conc[46]             # Particulate organic phosphate (mmol P m^-3)
    particOrganDetritus_NO_Si = conc[47]             # Particulate organic silicate (mmol Si m^-3)
    disInorgCarbon_IO_C = conc[48]             # Dissolved inorganic carbon(mg C m^-3)
    totalAlkalinity_IO = conc[49]             # Total alkalinity (mmol Eq m^-3)

    #--------------------------------------------------------------------------
    # concentration ratios
    diatoms_LO_N_diatoms_LO_C = get_concentration_ratio(diatoms_LO_N, diatoms_LO_C, constant_parameters["p_small"])
    diatoms_LO_P_diatoms_LO_C = get_concentration_ratio(diatoms_LO_P, diatoms_LO_C, constant_parameters["p_small"])
    diatoms_LO_Chl_diatoms_LO_C = get_concentration_ratio(diatoms_LO_Chl, diatoms_LO_C, constant_parameters["p_small"])
    diatoms_LO_Si_diatoms_LO_C = get_concentration_ratio(diatoms_LO_Si, diatoms_LO_C, constant_parameters["p_small"])
    nanoflagellates_LO_N_nanoflagellates_LO_C = get_concentration_ratio(nanoflagellates_LO_N, nanoflagellates_LO_C, constant_parameters["p_small"])
    nanoflagellates_LO_P_nanoflagellates_LO_C = get_concentration_ratio(nanoflagellates_LO_P, nanoflagellates_LO_C, constant_parameters["p_small"])
    nanoflagellates_LO_Chl_nanoflagellates_LO_C = get_concentration_ratio(nanoflagellates_LO_Chl, nanoflagellates_LO_C, constant_parameters["p_small"])
    picophyto_LO_N_picophyto_LO_C = get_concentration_ratio(picophyto_LO_N, picophyto_LO_C, constant_parameters["p_small"])
    picophyto_LO_P_picophyto_LO_C = get_concentration_ratio(picophyto_LO_P, picophyto_LO_C, constant_parameters["p_small"])
    picophyto_LO_Chl_picophyto_LO_C = get_concentration_ratio(picophyto_LO_Chl, picophyto_LO_C, constant_parameters["p_small"])
    largephyto_LO_N_largephyto_LO_C = get_concentration_ratio(largephyto_LO_N, largephyto_LO_C, constant_parameters["p_small"])
    largephyto_LO_P_largephyto_LO_C = get_concentration_ratio(largephyto_LO_P, largephyto_LO_C, constant_parameters["p_small"])
    largephyto_LO_Chl_largephyto_LO_C = get_concentration_ratio(largephyto_LO_Chl, largephyto_LO_C, constant_parameters["p_small"])
    bp_bc = get_concentration_ratio(pelBacteria_LO_P, pelBacteria_LO_C, constant_parameters["p_small"])
    bn_bc = get_concentration_ratio(pelBacteria_LO_N, pelBacteria_LO_C, constant_parameters["p_small"])
    carnivMesozoo_LO_N_carnivMesozoo_LO_C = get_concentration_ratio(carnivMesozoo_LO_N, carnivMesozoo_LO_C, constant_parameters["p_small"])
    carnivMesozoo_LO_P_carnivMesozoo_LO_C = get_concentration_ratio(carnivMesozoo_LO_P, carnivMesozoo_LO_C, constant_parameters["p_small"])
    omnivMesozoo_LO_N_omnivMesozoo_LO_C = get_concentration_ratio(omnivMesozoo_LO_N, omnivMesozoo_LO_C, constant_parameters["p_small"])
    omnivMesozoo_LO_P_omnivMesozoo_LO_C = get_concentration_ratio(omnivMesozoo_LO_P, omnivMesozoo_LO_C, constant_parameters["p_small"])
    microzoo_LO_N_microzoo_LO_C = get_concentration_ratio(microzoo_LO_N, microzoo_LO_C, constant_parameters["p_small"])
    microzoo_LO_P_microzoo_LO_C = get_concentration_ratio(microzoo_LO_P, microzoo_LO_C, constant_parameters["p_small"])
    heteroFlagellates_LO_N_heteroFlagellates_LO_C = get_concentration_ratio(heteroFlagellates_LO_N, heteroFlagellates_LO_C, constant_parameters["p_small"])
    heteroFlagellates_LO_P_heteroFlagellates_LO_C = get_concentration_ratio(heteroFlagellates_LO_P, heteroFlagellates_LO_C, constant_parameters["p_small"])
    labileDOM_NO_P_labileDOM_NO_C = get_concentration_ratio(labileDOM_NO_P, labileDOM_NO_C, constant_parameters["p_small"])
    particOrganDetritus_NO_P_particOrganDetritus_NO_C = get_concentration_ratio(particOrganDetritus_NO_P, particOrganDetritus_NO_C, constant_parameters["p_small"])
    labileDOM_NO_N_labileDOM_NO_C = get_concentration_ratio(labileDOM_NO_N, labileDOM_NO_C, constant_parameters["p_small"])
    particOrganDetritus_NO_N_particOrganDetritus_NO_C = get_concentration_ratio(particOrganDetritus_NO_N, particOrganDetritus_NO_C, constant_parameters["p_small"])

    #--------------------------------------------------------------------------
    #---------------------- Phytoplankton Equations ---------------------------
    #--------------------------------------------------------------------------
    # P1: Diatoms terms
    (ddiatoms_LO_Cdt_gpp_disInorgCarbon_IO_C, ddiatoms_LO_Cdt_rsp_disInorgCarbon_IO_C, ddiatoms_LO_Cdt_lys_labileDOM_NO_C, ddiatoms_LO_Cdt_lys_particOrganDetritus_NO_C, ddiatoms_LO_Cdt_exu_semilabileDOC_NO_C, ddiatoms_LO_Ndt_upt_nitrate_IO_N, ddiatoms_LO_Ndt_upt_ammonium_IO_N, 
     extra_n1, ddiatoms_LO_Ndt_lys_labileDOM_NO_N, ddiatoms_LO_Ndt_lys_particOrganDetritus_NO_N, ddiatoms_LO_Pdt_upt_phospate_IO_P, ddiatoms_LO_Pdt_upt_labileDOM_NO_P, ddiatoms_LO_Pdt_lys_labileDOM_NO_P, ddiatoms_LO_Pdt_lys_particOrganDetritus_NO_P, 
     ddiatoms_LO_Chldt_syn, ddiatoms_LO_Sidt_upt_silicate_IO_Si, ddiatoms_LO_Sidt_lys_particOrganDetritus_NO_Si) = phyto_eqns(conc, phyto1_prameters, environmental_parameters, constant_parameters, del_z, 1, irradiance, diatoms_LO_C, diatoms_LO_N, diatoms_LO_P, diatoms_LO_Chl, suspended_sediments, temper, time, xEPS)
    
    # (ddiatoms_LO_Cdt_gpp_disInorgCarbon_IO_C, ddiatoms_LO_Cdt_rsp_disInorgCarbon_IO_C, ddiatoms_LO_Cdt_lys_labileDOM_NO_C, ddiatoms_LO_Cdt_lys_particOrganDetritus_NO_C, ddiatoms_LO_Cdt_exu_semilabileDOC_NO_C, ddiatoms_LO_Ndt_upt_nitrate_IO_N, ddiatoms_LO_Ndt_upt_ammonium_IO_N, 
    #  extra_n1, ddiatoms_LO_Ndt_lys_labileDOM_NO_N, ddiatoms_LO_Ndt_lys_particOrganDetritus_NO_N, ddiatoms_LO_Pdt_upt_phospate_IO_P, ddiatoms_LO_Pdt_upt_labileDOM_NO_P, ddiatoms_LO_Pdt_lys_labileDOM_NO_P, ddiatoms_LO_Pdt_lys_particOrganDetritus_NO_P, 
    #  ddiatoms_LO_Chldt_syn, ddiatoms_LO_Sidt_upt_silicate_IO_Si, ddiatoms_LO_Sidt_lys_particOrganDetritus_NO_Si) = (0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.)

    # P2: Flagellates terms
    (dnanoflagellates_LO_Cdt_gpp_disInorgCarbon_IO_C, dnanoflagellates_LO_Cdt_rsp_disInorgCarbon_IO_C, dnanoflagellates_LO_Cdt_lys_labileDOM_NO_C, dnanoflagellates_LO_Cdt_lys_particOrganDetritus_NO_C, dnanoflagellates_LO_Cdt_exu_semilabileDOC_NO_C, dnanoflagellates_LO_Ndt_upt_nitrate_IO_N, dnanoflagellates_LO_Ndt_upt_ammonium_IO_N, 
     extra_n2, dnanoflagellates_LO_Ndt_lys_labileDOM_NO_N, dnanoflagellates_LO_Ndt_lys_particOrganDetritus_NO_N, dnanoflagellates_LO_Pdt_upt_phospate_IO_P, dnanoflagellates_LO_Pdt_upt_labileDOM_NO_P, dnanoflagellates_LO_Pdt_lys_labileDOM_NO_P, dnanoflagellates_LO_Pdt_lys_particOrganDetritus_NO_P, 
     dnanoflagellates_LO_Chldt_syn, dP2sdt_upt_silicate_IO_Si, dP2sdt_lys_particOrganDetritus_NO_Si) = phyto_eqns(conc, phyto2_prameters, environmental_parameters, constant_parameters, del_z, 2, irradiance, nanoflagellates_LO_C, nanoflagellates_LO_N, nanoflagellates_LO_P, nanoflagellates_LO_Chl, suspended_sediments, temper, time, xEPS)

    # P3: PicoPhytoplankton terms
    (dpicophyto_LO_Cdt_gpp_disInorgCarbon_IO_C, dpicophyto_LO_Cdt_rsp_disInorgCarbon_IO_C, dpicophyto_LO_Cdt_lys_labileDOM_NO_C, dpicophyto_LO_Cdt_lys_particOrganDetritus_NO_C, dpicophyto_LO_Cdt_exu_semilabileDOC_NO_C, dpicophyto_LO_Ndt_upt_nitrate_IO_N, dpicophyto_LO_Ndt_upt_ammonium_IO_N, 
     extra_n3, dpicophyto_LO_Ndt_lys_labileDOM_NO_N, dpicophyto_LO_Ndt_lys_particOrganDetritus_NO_N, dpicophyto_LO_Pdt_upt_phospate_IO_P, dpicophyto_LO_Pdt_upt_labileDOM_NO_P, dpicophyto_LO_Pdt_lys_labileDOM_NO_P, dpicophyto_LO_Pdt_lys_particOrganDetritus_NO_P, 
     dpicophyto_LO_Chldt_syn, dP3sdt_upt_silicate_IO_Si, dP3sdt_lys_particOrganDetritus_NO_Si) = phyto_eqns(conc, phyto3_prameters, environmental_parameters, constant_parameters, del_z, 3, irradiance, picophyto_LO_C, picophyto_LO_N, picophyto_LO_P, picophyto_LO_Chl, suspended_sediments, temper, time, xEPS)

    # (dpicophyto_LO_Cdt_gpp_disInorgCarbon_IO_C, dpicophyto_LO_Cdt_rsp_disInorgCarbon_IO_C, dpicophyto_LO_Cdt_lys_labileDOM_NO_C, dpicophyto_LO_Cdt_lys_particOrganDetritus_NO_C, dpicophyto_LO_Cdt_exu_semilabileDOC_NO_C, dpicophyto_LO_Ndt_upt_nitrate_IO_N, dpicophyto_LO_Ndt_upt_ammonium_IO_N, 
    #  extra_n3, dpicophyto_LO_Ndt_lys_labileDOM_NO_N, dpicophyto_LO_Ndt_lys_particOrganDetritus_NO_N, dpicophyto_LO_Pdt_upt_phospate_IO_P, dpicophyto_LO_Pdt_upt_labileDOM_NO_P, dpicophyto_LO_Pdt_lys_labileDOM_NO_P, dpicophyto_LO_Pdt_lys_particOrganDetritus_NO_P, 
    #  dpicophyto_LO_Chldt_syn, dP3sdt_upt_silicate_IO_Si, dP3sdt_lys_particOrganDetritus_NO_Si) = (0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.)

    # P4: Large Phytoplankton terms
    (dlargephyto_LO_Cdt_gpp_disInorgCarbon_IO_C, dlargephyto_LO_Cdt_rsp_disInorgCarbon_IO_C, dlargephyto_LO_Cdt_lys_labileDOM_NO_C, dlargephyto_LO_Cdt_lys_particOrganDetritus_NO_C, dlargephyto_LO_Cdt_exu_semilabileDOC_NO_C, dlargephyto_LO_Ndt_upt_nitrate_IO_N, dlargephyto_LO_Ndt_upt_ammonium_IO_N, 
     extra_n4, dlargephyto_LO_Ndt_lys_labileDOM_NO_N, dlargephyto_LO_Ndt_lys_particOrganDetritus_NO_N, dlargephyto_LO_Pdt_upt_phospate_IO_P, dlargephyto_LO_Pdt_upt_labileDOM_NO_P, dlargephyto_LO_Pdt_lys_labileDOM_NO_P, dlargephyto_LO_Pdt_lys_particOrganDetritus_NO_P, 
     dlargephyto_LO_Chldt_syn, dP4sdt_upt_silicate_IO_Si, dP4sdt_lys_particOrganDetritus_NO_Si) = phyto_eqns(conc, phyto4_prameters, environmental_parameters, constant_parameters, del_z, 4, irradiance, largephyto_LO_C, largephyto_LO_N, largephyto_LO_P, largephyto_LO_Chl, suspended_sediments, temper, time, xEPS)

    # (dlargephyto_LO_Cdt_gpp_disInorgCarbon_IO_C, dlargephyto_LO_Cdt_rsp_disInorgCarbon_IO_C, dlargephyto_LO_Cdt_lys_labileDOM_NO_C, dlargephyto_LO_Cdt_lys_particOrganDetritus_NO_C, dlargephyto_LO_Cdt_exu_semilabileDOC_NO_C, dlargephyto_LO_Ndt_upt_nitrate_IO_N, dlargephyto_LO_Ndt_upt_ammonium_IO_N, 
    #  extra_n4, dlargephyto_LO_Ndt_lys_labileDOM_NO_N, dlargephyto_LO_Ndt_lys_particOrganDetritus_NO_N, dlargephyto_LO_Pdt_upt_phospate_IO_P, dlargephyto_LO_Pdt_upt_labileDOM_NO_P, dlargephyto_LO_Pdt_lys_labileDOM_NO_P, dlargephyto_LO_Pdt_lys_particOrganDetritus_NO_P, 
    #  dlargephyto_LO_Chldt_syn, dP4sdt_upt_silicate_IO_Si, dP4sdt_lys_particOrganDetritus_NO_Si) = (0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.)

    #--------------------------------------------------------------------------
    #------------------------- Bacteria Equations -----------------------------
    #--------------------------------------------------------------------------
    (dBcdt_lys_labileDOM_NO_C, dBcdt_lys_labileDOM_NO_N, dBcdt_lys_labileDOM_NO_P, dBcdt_lys_particOrganDetritus_NO_C, dBcdt_lys_particOrganDetritus_NO_N, dBcdt_lys_particOrganDetritus_NO_P, 
     dBcdt_upt_labileDOM_NO_C, dBcdt_upt_particOrganDetritus_NO_C, dBpdt_upt_rel_phospate_IO_P, dBndt_upt_rel_ammonium_IO_N, dBcdt_upt_semilabileDOC_NO_C, dBcdt_upt_semirefractDOC_NO_C, 
     dBcdt_rel_semilabileDOC_NO_C, dBcdt_rel_semirefractDOC_NO_C, dBcdt_rsp_disInorgCarbon_IO_C, flPTreductEquiv_IO_R, f_B_O, f_B_n, f_B_p) = bacteria_eqns(conc, bacteria_parameters, constant_parameters, environmental_parameters, temper)

    # (dBcdt_lys_labileDOM_NO_C, dBcdt_lys_labileDOM_NO_N, dBcdt_lys_labileDOM_NO_P, dBcdt_lys_particOrganDetritus_NO_C, dBcdt_lys_particOrganDetritus_NO_N, dBcdt_lys_particOrganDetritus_NO_P, 
    #  dBcdt_upt_labileDOM_NO_C, dBcdt_upt_particOrganDetritus_NO_C, dBpdt_upt_rel_phospate_IO_P, dBndt_upt_rel_ammonium_IO_N, dBcdt_upt_semilabileDOC_NO_C, dBcdt_upt_semirefractDOC_NO_C, 
    #  dBcdt_rel_semilabileDOC_NO_C, dBcdt_rel_semirefractDOC_NO_C, dBcdt_rsp_disInorgCarbon_IO_C, flPTreductEquiv_IO_R, f_B_O, f_B_n, f_B_p) = (0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.)
    #--------------------------------------------------------------------------
    #----------------------- Zooplankton Equations ----------------------------
    #--------------------------------------------------------------------------
    
    # Mesozooplankton predation terms
    dcarnivMesozoo_LO_Cdt_prd, domnivMesozoo_LO_Cdt_prd, ic3, in3, ip3, ic4, in4, ip4 = get_mesozoo_predation_terms(conc, mesozoo3_parameters, mesozoo4_parameters, zoo_availability_parameters, environmental_parameters, constant_parameters, temper)

    # Microzooplankton predation terms
    dmicrozoo_LO_Cdt_prd, dheteroFlagellates_LO_Cdt_prd, ic5, in5, ip5, ic6, in6, ip6 = get_microzoo_predation_terms(conc, microzoo5_parameters, microzoo6_parameters, zoo_availability_parameters, environmental_parameters, constant_parameters, temper)

    # Z3: Carnivorous Mesozooplankton terms
    (dcarnivMesozoo_LO_Cdt_rel_labileDOM_NO_C, dcarnivMesozoo_LO_Cdt_rel_particOrganDetritus_NO_C, dcarnivMesozoo_LO_Cdt_rsp_disInorgCarbon_IO_C, dcarnivMesozoo_LO_Ndt_rel_labileDOM_NO_N, dcarnivMesozoo_LO_Ndt_rel_particOrganDetritus_NO_N, dcarnivMesozoo_LO_Pdt_rel_labileDOM_NO_P, dcarnivMesozoo_LO_Pdt_rel_particOrganDetritus_NO_P, 
     dcarnivMesozoo_LO_Pdt_rel_phospate_IO_P, dcarnivMesozoo_LO_Ndt_rel_ammonium_IO_N) = mesozoo_eqns(conc, mesozoo3_parameters, constant_parameters, environmental_parameters, carnivMesozoo_LO_C, carnivMesozoo_LO_N, carnivMesozoo_LO_P, ic3, in3, ip3, temper)

    # (dcarnivMesozoo_LO_Cdt_rel_labileDOM_NO_C, dcarnivMesozoo_LO_Cdt_rel_particOrganDetritus_NO_C, dcarnivMesozoo_LO_Cdt_rsp_disInorgCarbon_IO_C, dcarnivMesozoo_LO_Ndt_rel_labileDOM_NO_N, dcarnivMesozoo_LO_Ndt_rel_particOrganDetritus_NO_N, dcarnivMesozoo_LO_Pdt_rel_labileDOM_NO_P, dcarnivMesozoo_LO_Pdt_rel_particOrganDetritus_NO_P, 
    #  dcarnivMesozoo_LO_Pdt_rel_phospate_IO_P, dcarnivMesozoo_LO_Ndt_rel_ammonium_IO_N) = (0.,0.,0.,0.,0.,0.,0.,0.,0.)
    
    # Z4: Omnivorous Mesozooplankton terms
    (domnivMesozoo_LO_Cdt_rel_labileDOM_NO_C, domnivMesozoo_LO_Cdt_rel_particOrganDetritus_NO_C, domnivMesozoo_LO_Cdt_rsp_disInorgCarbon_IO_C, domnivMesozoo_LO_Ndt_rel_labileDOM_NO_N, domnivMesozoo_LO_Ndt_rel_particOrganDetritus_NO_N, domnivMesozoo_LO_Pdt_rel_labileDOM_NO_P, domnivMesozoo_LO_Pdt_rel_particOrganDetritus_NO_P, 
     domnivMesozoo_LO_Pdt_rel_phospate_IO_P, domnivMesozoo_LO_Ndt_rel_ammonium_IO_N) = mesozoo_eqns(conc, mesozoo4_parameters, constant_parameters, environmental_parameters, omnivMesozoo_LO_C, omnivMesozoo_LO_N, omnivMesozoo_LO_P, ic4, in4, ip4, temper)

    # (domnivMesozoo_LO_Cdt_rel_labileDOM_NO_C, domnivMesozoo_LO_Cdt_rel_particOrganDetritus_NO_C, domnivMesozoo_LO_Cdt_rsp_disInorgCarbon_IO_C, domnivMesozoo_LO_Ndt_rel_labileDOM_NO_N, domnivMesozoo_LO_Ndt_rel_particOrganDetritus_NO_N, domnivMesozoo_LO_Pdt_rel_labileDOM_NO_P, domnivMesozoo_LO_Pdt_rel_particOrganDetritus_NO_P, 
    #  domnivMesozoo_LO_Pdt_rel_phospate_IO_P, domnivMesozoo_LO_Ndt_rel_ammonium_IO_N) = (0.,0.,0.,0.,0.,0.,0.,0.,0.)
    
    # Z5: Microzooplankton terms
    (dmicrozoo_LO_Cdt_rel_labileDOM_NO_C, dmicrozoo_LO_Cdt_rel_particOrganDetritus_NO_C, dmicrozoo_LO_Cdt_rsp_disInorgCarbon_IO_C, dmicrozoo_LO_Ndt_rel_labileDOM_NO_N, dmicrozoo_LO_Ndt_rel_particOrganDetritus_NO_N, dmicrozoo_LO_Pdt_rel_labileDOM_NO_P, dmicrozoo_LO_Pdt_rel_particOrganDetritus_NO_P, 
     dmicrozoo_LO_Pdt_rel_phospate_IO_P, dmicrozoo_LO_Ndt_rel_ammonium_IO_N) = microzoo_eqns(conc, microzoo5_parameters, constant_parameters, environmental_parameters, microzoo_LO_C, microzoo_LO_N, microzoo_LO_P, ic5, in5, ip5, temper)

    # Z6: Heterotrophic Nanoflagellates terms
    (dheteroFlagellates_LO_Cdt_rel_labileDOM_NO_C, dheteroFlagellates_LO_Cdt_rel_particOrganDetritus_NO_C, dheteroFlagellates_LO_Cdt_rsp_disInorgCarbon_IO_C, dheteroFlagellates_LO_Ndt_rel_labileDOM_NO_N, dheteroFlagellates_LO_Ndt_rel_particOrganDetritus_NO_N, dheteroFlagellates_LO_Pdt_rel_labileDOM_NO_P, dheteroFlagellates_LO_Pdt_rel_particOrganDetritus_NO_P,  
     dheteroFlagellates_LO_Pdt_rel_phospate_IO_P, dheteroFlagellates_LO_Ndt_rel_ammonium_IO_N) = microzoo_eqns(conc, microzoo6_parameters, constant_parameters, environmental_parameters, heteroFlagellates_LO_C, heteroFlagellates_LO_N, heteroFlagellates_LO_P, ic6, in6, ip6, temper)

    # (dheteroFlagellates_LO_Cdt_rel_labileDOM_NO_C, dheteroFlagellates_LO_Cdt_rel_particOrganDetritus_NO_C, dheteroFlagellates_LO_Cdt_rsp_disInorgCarbon_IO_C, dheteroFlagellates_LO_Ndt_rel_labileDOM_NO_N, dheteroFlagellates_LO_Ndt_rel_particOrganDetritus_NO_N, dheteroFlagellates_LO_Pdt_rel_labileDOM_NO_P, dheteroFlagellates_LO_Pdt_rel_particOrganDetritus_NO_P,  
    #  dheteroFlagellates_LO_Pdt_rel_phospate_IO_P, dheteroFlagellates_LO_Ndt_rel_ammonium_IO_N) = (0.,0.,0.,0.,0.,0.,0.,0.,0.)

    #--------------------------------------------------------------------------
    #------------------------ Non-living components ---------------------------
    #--------------------------------------------------------------------------    
    (dammonium_IO_Ndt_nit_nitrate_IO_N, dnitrate_IO_Ndt_denit, dreductEquiv_IO_Rdt_reox, dparticOrganDetritus_NO_Sidt_rmn_silicate_IO_Si) = pel_chem_eqns(pel_chem_parameters, environmental_parameters, constant_parameters, temper, conc, flPTreductEquiv_IO_R)

    #---------------------- Oxygen airation by wind ---------------------------
    
    dOdt_wind = calculate_oxygen_reaeration(oxygen_reaeration_parameters, environmental_parameters, constant_parameters, conc, temper, salt, wind)

    #------------------------------- CO_2 Flux --------------------------------
    # do3cdt_air_sea_flux = calculate_co2_flux(co2_flux_parameters, environmental_parameters, constant_parameters, conc, temper, wind, salt)
    do3cdt_air_sea_flux = calculate_co2_flux(co2_flux_parameters, environmental_parameters, constant_parameters, conc, temper, wind, salt, rho, del_z)

    #--------------------------------------------------------------------------
    #----------------------------- Rate Equations -----------------------------
    #--------------------------------------------------------------------------
    
    # Dissolved oxygen [mmol O_2 m^-3 s^-1]
    # ddisOxygen_IO_O_dt = (constant_parameters["omega_c"]*((ddiatoms_LO_Cdt_gpp_disInorgCarbon_IO_C - ddiatoms_LO_Cdt_rsp_disInorgCarbon_IO_C) + (dnanoflagellates_LO_Cdt_gpp_disInorgCarbon_IO_C - dnanoflagellates_LO_Cdt_rsp_disInorgCarbon_IO_C) + (dpicophyto_LO_Cdt_gpp_disInorgCarbon_IO_C - dpicophyto_LO_Cdt_rsp_disInorgCarbon_IO_C) + (dlargephyto_LO_Cdt_gpp_disInorgCarbon_IO_C - dlargephyto_LO_Cdt_rsp_disInorgCarbon_IO_C)) - 
    #            constant_parameters["omega_c"]*f_B_O*dBcdt_rsp_disInorgCarbon_IO_C - 
    #            constant_parameters["omega_c"]*(dcarnivMesozoo_LO_Cdt_rsp_disInorgCarbon_IO_C + domnivMesozoo_LO_Cdt_rsp_disInorgCarbon_IO_C + dmicrozoo_LO_Cdt_rsp_disInorgCarbon_IO_C + dheteroFlagellates_LO_Cdt_rsp_disInorgCarbon_IO_C) - 
    #            constant_parameters["omega_n"]*dammonium_IO_Ndt_nit_nitrate_IO_N -
    #            (1.0/constant_parameters["omega_r"])*dreductEquiv_IO_Rdt_reox) + dOdt_wind

    # ddisOxygen_IO_O_dt = (constant_parameters["omega_c"]*((ddiatoms_LO_Cdt_gpp_disInorgCarbon_IO_C - ddiatoms_LO_Cdt_rsp_disInorgCarbon_IO_C) + (dnanoflagellates_LO_Cdt_gpp_disInorgCarbon_IO_C - dnanoflagellates_LO_Cdt_rsp_disInorgCarbon_IO_C) + (dpicophyto_LO_Cdt_gpp_disInorgCarbon_IO_C - dpicophyto_LO_Cdt_rsp_disInorgCarbon_IO_C) + (dlargephyto_LO_Cdt_gpp_disInorgCarbon_IO_C - dlargephyto_LO_Cdt_rsp_disInorgCarbon_IO_C)) - 
    #            constant_parameters["omega_c"]*f_B_O*dBcdt_rsp_disInorgCarbon_IO_C - 
    #            constant_parameters["omega_c"]*(dcarnivMesozoo_LO_Cdt_rsp_disInorgCarbon_IO_C + domnivMesozoo_LO_Cdt_rsp_disInorgCarbon_IO_C + dmicrozoo_LO_Cdt_rsp_disInorgCarbon_IO_C + dheteroFlagellates_LO_Cdt_rsp_disInorgCarbon_IO_C) - 
    #            constant_parameters["omega_n"]*dammonium_IO_Ndt_nit_nitrate_IO_N -
    #            (1.0/constant_parameters["omega_r"])*dreductEquiv_IO_Rdt_reox)
    
    ddisOxygen_IO_O_dt = (constant_parameters["omega_c"]*((ddiatoms_LO_Cdt_gpp_disInorgCarbon_IO_C - ddiatoms_LO_Cdt_rsp_disInorgCarbon_IO_C) + (dnanoflagellates_LO_Cdt_gpp_disInorgCarbon_IO_C - dnanoflagellates_LO_Cdt_rsp_disInorgCarbon_IO_C) + (dpicophyto_LO_Cdt_gpp_disInorgCarbon_IO_C - dpicophyto_LO_Cdt_rsp_disInorgCarbon_IO_C) + (dlargephyto_LO_Cdt_gpp_disInorgCarbon_IO_C - dlargephyto_LO_Cdt_rsp_disInorgCarbon_IO_C)) - 
               constant_parameters["omega_c"]*f_B_O*dBcdt_rsp_disInorgCarbon_IO_C - 
               constant_parameters["omega_c"]*(dcarnivMesozoo_LO_Cdt_rsp_disInorgCarbon_IO_C + domnivMesozoo_LO_Cdt_rsp_disInorgCarbon_IO_C + dmicrozoo_LO_Cdt_rsp_disInorgCarbon_IO_C + dheteroFlagellates_LO_Cdt_rsp_disInorgCarbon_IO_C) - 
               constant_parameters["omega_n"]*dammonium_IO_Ndt_nit_nitrate_IO_N -
               (1.0/constant_parameters["omega_r"])*dreductEquiv_IO_Rdt_reox)

    # Dissolved inorganic nutrient equations
    dphospate_IO_P_dt = - (ddiatoms_LO_Pdt_upt_phospate_IO_P + dnanoflagellates_LO_Pdt_upt_phospate_IO_P + dpicophyto_LO_Pdt_upt_phospate_IO_P + dlargephyto_LO_Pdt_upt_phospate_IO_P) + (dBpdt_upt_rel_phospate_IO_P*insw_vector(dBpdt_upt_rel_phospate_IO_P)) - ((-1)*f_B_p*dBpdt_upt_rel_phospate_IO_P*insw_vector(-dBpdt_upt_rel_phospate_IO_P)) + (dcarnivMesozoo_LO_Pdt_rel_phospate_IO_P + domnivMesozoo_LO_Pdt_rel_phospate_IO_P + dmicrozoo_LO_Pdt_rel_phospate_IO_P + dheteroFlagellates_LO_Pdt_rel_phospate_IO_P)
    dnitrate_IO_N_dt = - (ddiatoms_LO_Ndt_upt_nitrate_IO_N + dnanoflagellates_LO_Ndt_upt_nitrate_IO_N + dpicophyto_LO_Ndt_upt_nitrate_IO_N + dlargephyto_LO_Ndt_upt_nitrate_IO_N) + dammonium_IO_Ndt_nit_nitrate_IO_N - dnitrate_IO_Ndt_denit
    dammonium_IO_N_dt = - (ddiatoms_LO_Ndt_upt_ammonium_IO_N + dnanoflagellates_LO_Ndt_upt_ammonium_IO_N + dpicophyto_LO_Ndt_upt_ammonium_IO_N + dlargephyto_LO_Ndt_upt_ammonium_IO_N) + (dBndt_upt_rel_ammonium_IO_N*insw_vector(dBndt_upt_rel_ammonium_IO_N)) - ((-1)*f_B_n*dBndt_upt_rel_ammonium_IO_N*insw_vector(-dBndt_upt_rel_ammonium_IO_N)) + (dcarnivMesozoo_LO_Ndt_rel_ammonium_IO_N + domnivMesozoo_LO_Ndt_rel_ammonium_IO_N + dmicrozoo_LO_Ndt_rel_ammonium_IO_N + dheteroFlagellates_LO_Ndt_rel_ammonium_IO_N) - dammonium_IO_Ndt_nit_nitrate_IO_N
    dnitrogenSink_dt = dnitrate_IO_Ndt_denit
    dsilicate_IO_Si_dt = - ddiatoms_LO_Sidt_upt_silicate_IO_Si + dparticOrganDetritus_NO_Sidt_rmn_silicate_IO_Si

    # Reduction equivalents
    dreductEquiv_IO_R_dt = constant_parameters["omega_r"]*constant_parameters["omega_c"]*(1.0 - f_B_O)*dBcdt_rsp_disInorgCarbon_IO_C - constant_parameters["omega_r"]*constant_parameters["omega_dn"]*dnitrate_IO_Ndt_denit*insw_vector(-(disOxygen_IO_O - reductEquiv_IO_R)/constant_parameters["omega_r"]) - dreductEquiv_IO_Rdt_reox

    # Bacterioplankton
    dpelBacteria_LO_C_dt = dBcdt_upt_labileDOM_NO_C + dBcdt_upt_particOrganDetritus_NO_C - dBcdt_rsp_disInorgCarbon_IO_C - dBcdt_lys_labileDOM_NO_C - dBcdt_lys_particOrganDetritus_NO_C - (dmicrozoo_LO_Cdt_prd["b1"] + dheteroFlagellates_LO_Cdt_prd["b1"])
    dpelBacteria_LO_N_dt = - dBcdt_lys_labileDOM_NO_N - dBcdt_lys_particOrganDetritus_NO_N + (labileDOM_NO_N_labileDOM_NO_C*dBcdt_upt_labileDOM_NO_C) + (particOrganDetritus_NO_N_particOrganDetritus_NO_C*dBcdt_upt_particOrganDetritus_NO_C) - (dBndt_upt_rel_ammonium_IO_N*insw_vector(dBndt_upt_rel_ammonium_IO_N)) + ((-1)*f_B_n*dBndt_upt_rel_ammonium_IO_N*insw_vector(-dBndt_upt_rel_ammonium_IO_N)) - (bn_bc*(dmicrozoo_LO_Cdt_prd["b1"] + dheteroFlagellates_LO_Cdt_prd["b1"]))
    dpelBacteria_LO_P_dt = (labileDOM_NO_P_labileDOM_NO_C*dBcdt_upt_labileDOM_NO_C) + (particOrganDetritus_NO_P_particOrganDetritus_NO_C*dBcdt_upt_particOrganDetritus_NO_C) - (dBpdt_upt_rel_phospate_IO_P*insw_vector(dBpdt_upt_rel_phospate_IO_P)) + ((-1)*f_B_p*dBpdt_upt_rel_phospate_IO_P*insw_vector(-dBpdt_upt_rel_phospate_IO_P)) - dBcdt_lys_labileDOM_NO_P - dBcdt_lys_particOrganDetritus_NO_P - (bp_bc*(dmicrozoo_LO_Cdt_prd["b1"] + dheteroFlagellates_LO_Cdt_prd["b1"]))

    # Phytoplankton
    ddiatoms_LO_C_dt = ddiatoms_LO_Cdt_gpp_disInorgCarbon_IO_C - ddiatoms_LO_Cdt_exu_semilabileDOC_NO_C - ddiatoms_LO_Cdt_rsp_disInorgCarbon_IO_C - ddiatoms_LO_Cdt_lys_labileDOM_NO_C - ddiatoms_LO_Cdt_lys_particOrganDetritus_NO_C - dcarnivMesozoo_LO_Cdt_prd["p1"] - domnivMesozoo_LO_Cdt_prd["p1"] - dmicrozoo_LO_Cdt_prd["p1"] - dheteroFlagellates_LO_Cdt_prd["p1"]
    ddiatoms_LO_N_dt = ddiatoms_LO_Ndt_upt_nitrate_IO_N + ddiatoms_LO_Ndt_upt_ammonium_IO_N - extra_n1 - ddiatoms_LO_Ndt_lys_labileDOM_NO_N - ddiatoms_LO_Ndt_lys_particOrganDetritus_NO_N - (diatoms_LO_N_diatoms_LO_C*(dcarnivMesozoo_LO_Cdt_prd["p1"] + domnivMesozoo_LO_Cdt_prd["p1"] + dmicrozoo_LO_Cdt_prd["p1"] + dheteroFlagellates_LO_Cdt_prd["p1"]))
    ddiatoms_LO_P_dt = ddiatoms_LO_Pdt_upt_phospate_IO_P - ddiatoms_LO_Pdt_upt_labileDOM_NO_P - ddiatoms_LO_Pdt_lys_labileDOM_NO_P - ddiatoms_LO_Pdt_lys_particOrganDetritus_NO_P - (diatoms_LO_P_diatoms_LO_C*(dcarnivMesozoo_LO_Cdt_prd["p1"] + domnivMesozoo_LO_Cdt_prd["p1"] + dmicrozoo_LO_Cdt_prd["p1"] + dheteroFlagellates_LO_Cdt_prd["p1"]))
    ddiatoms_LO_Chl_dt = ddiatoms_LO_Chldt_syn - (diatoms_LO_Chl_diatoms_LO_C*(dcarnivMesozoo_LO_Cdt_prd["p1"] + domnivMesozoo_LO_Cdt_prd["p1"] + dmicrozoo_LO_Cdt_prd["p1"] + dheteroFlagellates_LO_Cdt_prd["p1"]))
    ddiatoms_LO_Si_dt = ddiatoms_LO_Sidt_upt_silicate_IO_Si - ddiatoms_LO_Sidt_lys_particOrganDetritus_NO_Si - (diatoms_LO_Si_diatoms_LO_C*(dcarnivMesozoo_LO_Cdt_prd["p1"] + domnivMesozoo_LO_Cdt_prd["p1"] + dmicrozoo_LO_Cdt_prd["p1"] + dheteroFlagellates_LO_Cdt_prd["p1"]))
    
    dnanoflagellates_LO_C_dt = dnanoflagellates_LO_Cdt_gpp_disInorgCarbon_IO_C - dnanoflagellates_LO_Cdt_exu_semilabileDOC_NO_C - dnanoflagellates_LO_Cdt_rsp_disInorgCarbon_IO_C - dnanoflagellates_LO_Cdt_lys_labileDOM_NO_C - dnanoflagellates_LO_Cdt_lys_particOrganDetritus_NO_C - dcarnivMesozoo_LO_Cdt_prd["p2"] - domnivMesozoo_LO_Cdt_prd["p2"] - dmicrozoo_LO_Cdt_prd["p2"] - dheteroFlagellates_LO_Cdt_prd["p2"]
    dnanoflagellates_LO_N_dt = dnanoflagellates_LO_Ndt_upt_nitrate_IO_N + dnanoflagellates_LO_Ndt_upt_ammonium_IO_N - extra_n2 - dnanoflagellates_LO_Ndt_lys_labileDOM_NO_N - dnanoflagellates_LO_Ndt_lys_particOrganDetritus_NO_N - (nanoflagellates_LO_N_nanoflagellates_LO_C*(dcarnivMesozoo_LO_Cdt_prd["p2"] + domnivMesozoo_LO_Cdt_prd["p2"] + dmicrozoo_LO_Cdt_prd["p2"] + dheteroFlagellates_LO_Cdt_prd["p2"]))
    dnanoflagellates_LO_P_dt = dnanoflagellates_LO_Pdt_upt_phospate_IO_P - dnanoflagellates_LO_Pdt_upt_labileDOM_NO_P - dnanoflagellates_LO_Pdt_lys_labileDOM_NO_P - dnanoflagellates_LO_Pdt_lys_particOrganDetritus_NO_P - (nanoflagellates_LO_P_nanoflagellates_LO_C*(dcarnivMesozoo_LO_Cdt_prd["p2"] + domnivMesozoo_LO_Cdt_prd["p2"] + dmicrozoo_LO_Cdt_prd["p2"] + dheteroFlagellates_LO_Cdt_prd["p2"]))
    dnanoflagellates_LO_Chl_dt = dnanoflagellates_LO_Chldt_syn - (nanoflagellates_LO_Chl_nanoflagellates_LO_C*(dcarnivMesozoo_LO_Cdt_prd["p2"] + domnivMesozoo_LO_Cdt_prd["p2"] + dmicrozoo_LO_Cdt_prd["p2"] + dheteroFlagellates_LO_Cdt_prd["p2"]))
    
    dpicophyto_LO_C_dt = dpicophyto_LO_Cdt_gpp_disInorgCarbon_IO_C - dpicophyto_LO_Cdt_exu_semilabileDOC_NO_C - dpicophyto_LO_Cdt_rsp_disInorgCarbon_IO_C - dpicophyto_LO_Cdt_lys_labileDOM_NO_C - dpicophyto_LO_Cdt_lys_particOrganDetritus_NO_C - dcarnivMesozoo_LO_Cdt_prd["p3"] - domnivMesozoo_LO_Cdt_prd["p3"] - dmicrozoo_LO_Cdt_prd["p3"] - dheteroFlagellates_LO_Cdt_prd["p3"]
    dpicophyto_LO_N_dt = dpicophyto_LO_Ndt_upt_nitrate_IO_N + dpicophyto_LO_Ndt_upt_ammonium_IO_N - extra_n3 - dpicophyto_LO_Ndt_lys_labileDOM_NO_N - dpicophyto_LO_Ndt_lys_particOrganDetritus_NO_N - (picophyto_LO_N_picophyto_LO_C*(dcarnivMesozoo_LO_Cdt_prd["p3"] + domnivMesozoo_LO_Cdt_prd["p3"] + dmicrozoo_LO_Cdt_prd["p3"] + dheteroFlagellates_LO_Cdt_prd["p3"]))
    dpicophyto_LO_P_dt = dpicophyto_LO_Pdt_upt_phospate_IO_P - dpicophyto_LO_Pdt_upt_labileDOM_NO_P - dpicophyto_LO_Pdt_lys_labileDOM_NO_P - dpicophyto_LO_Pdt_lys_particOrganDetritus_NO_P - (picophyto_LO_P_picophyto_LO_C*(dcarnivMesozoo_LO_Cdt_prd["p3"] + domnivMesozoo_LO_Cdt_prd["p3"] + dmicrozoo_LO_Cdt_prd["p3"] + dheteroFlagellates_LO_Cdt_prd["p3"]))
    dpicophyto_LO_Chl_dt = dpicophyto_LO_Chldt_syn - (picophyto_LO_Chl_picophyto_LO_C*(dcarnivMesozoo_LO_Cdt_prd["p3"] + domnivMesozoo_LO_Cdt_prd["p3"] + dmicrozoo_LO_Cdt_prd["p3"] + dheteroFlagellates_LO_Cdt_prd["p3"]))
    
    dlargephyto_LO_C_dt = dlargephyto_LO_Cdt_gpp_disInorgCarbon_IO_C - dlargephyto_LO_Cdt_exu_semilabileDOC_NO_C - dlargephyto_LO_Cdt_rsp_disInorgCarbon_IO_C - dlargephyto_LO_Cdt_lys_labileDOM_NO_C - dlargephyto_LO_Cdt_lys_particOrganDetritus_NO_C - dcarnivMesozoo_LO_Cdt_prd["p4"] - domnivMesozoo_LO_Cdt_prd["p4"] - dmicrozoo_LO_Cdt_prd["p4"] - dheteroFlagellates_LO_Cdt_prd["p4"]
    dlargephyto_LO_N_dt = dlargephyto_LO_Ndt_upt_nitrate_IO_N + dlargephyto_LO_Ndt_upt_ammonium_IO_N - extra_n4 - dlargephyto_LO_Ndt_lys_labileDOM_NO_N - dlargephyto_LO_Ndt_lys_particOrganDetritus_NO_N - (largephyto_LO_N_largephyto_LO_C*(dcarnivMesozoo_LO_Cdt_prd["p4"] + domnivMesozoo_LO_Cdt_prd["p4"] + dmicrozoo_LO_Cdt_prd["p4"] + dheteroFlagellates_LO_Cdt_prd["p4"]))
    dlargephyto_LO_P_dt = dlargephyto_LO_Pdt_upt_phospate_IO_P - dlargephyto_LO_Pdt_upt_labileDOM_NO_P - dlargephyto_LO_Pdt_lys_labileDOM_NO_P - dlargephyto_LO_Pdt_lys_particOrganDetritus_NO_P - (largephyto_LO_P_largephyto_LO_C*(dcarnivMesozoo_LO_Cdt_prd["p4"] + domnivMesozoo_LO_Cdt_prd["p4"] + dmicrozoo_LO_Cdt_prd["p4"] + dheteroFlagellates_LO_Cdt_prd["p4"]))
    dlargephyto_LO_Chl_dt = dlargephyto_LO_Chldt_syn - (largephyto_LO_Chl_largephyto_LO_C*(dcarnivMesozoo_LO_Cdt_prd["p4"] + domnivMesozoo_LO_Cdt_prd["p4"] + dmicrozoo_LO_Cdt_prd["p4"] + dheteroFlagellates_LO_Cdt_prd["p4"]))

    # mesozooplankton
    dcarnivMesozoo_LO_C_dt = dcarnivMesozoo_LO_Cdt_prd["p1"] + dcarnivMesozoo_LO_Cdt_prd["p2"] + dcarnivMesozoo_LO_Cdt_prd["p3"] + dcarnivMesozoo_LO_Cdt_prd["p4"] + dcarnivMesozoo_LO_Cdt_prd["z4"] + dcarnivMesozoo_LO_Cdt_prd["z5"] + dcarnivMesozoo_LO_Cdt_prd["z6"] - domnivMesozoo_LO_Cdt_prd["z3"] - dcarnivMesozoo_LO_Cdt_rel_labileDOM_NO_C - dcarnivMesozoo_LO_Cdt_rel_particOrganDetritus_NO_C - dcarnivMesozoo_LO_Cdt_rsp_disInorgCarbon_IO_C
    dcarnivMesozoo_LO_N_dt = diatoms_LO_N_diatoms_LO_C*dcarnivMesozoo_LO_Cdt_prd["p1"] + nanoflagellates_LO_N_nanoflagellates_LO_C*dcarnivMesozoo_LO_Cdt_prd["p2"] + picophyto_LO_N_picophyto_LO_C*dcarnivMesozoo_LO_Cdt_prd["p3"] + largephyto_LO_N_largephyto_LO_C*dcarnivMesozoo_LO_Cdt_prd["p4"] + omnivMesozoo_LO_N_omnivMesozoo_LO_C*dcarnivMesozoo_LO_Cdt_prd["z4"] + microzoo_LO_N_microzoo_LO_C*dcarnivMesozoo_LO_Cdt_prd["z5"] + heteroFlagellates_LO_N_heteroFlagellates_LO_C*dcarnivMesozoo_LO_Cdt_prd["z6"] - carnivMesozoo_LO_N_carnivMesozoo_LO_C*domnivMesozoo_LO_Cdt_prd["z3"] - dcarnivMesozoo_LO_Ndt_rel_labileDOM_NO_N - dcarnivMesozoo_LO_Ndt_rel_particOrganDetritus_NO_N - dcarnivMesozoo_LO_Ndt_rel_ammonium_IO_N
    dcarnivMesozoo_LO_P_dt = diatoms_LO_P_diatoms_LO_C*dcarnivMesozoo_LO_Cdt_prd["p1"] + nanoflagellates_LO_P_nanoflagellates_LO_C*dcarnivMesozoo_LO_Cdt_prd["p2"] + picophyto_LO_P_picophyto_LO_C*dcarnivMesozoo_LO_Cdt_prd["p3"] + largephyto_LO_P_largephyto_LO_C*dcarnivMesozoo_LO_Cdt_prd["p4"] + omnivMesozoo_LO_P_omnivMesozoo_LO_C*dcarnivMesozoo_LO_Cdt_prd["z4"] + microzoo_LO_P_microzoo_LO_C*dcarnivMesozoo_LO_Cdt_prd["z5"] + heteroFlagellates_LO_P_heteroFlagellates_LO_C*dcarnivMesozoo_LO_Cdt_prd["z6"] - carnivMesozoo_LO_P_carnivMesozoo_LO_C*domnivMesozoo_LO_Cdt_prd["z3"] - dcarnivMesozoo_LO_Pdt_rel_labileDOM_NO_P - dcarnivMesozoo_LO_Pdt_rel_particOrganDetritus_NO_P - dcarnivMesozoo_LO_Pdt_rel_phospate_IO_P
    
    domnivMesozoo_LO_C_dt = domnivMesozoo_LO_Cdt_prd["p1"] + domnivMesozoo_LO_Cdt_prd["p2"] + domnivMesozoo_LO_Cdt_prd["p3"] + domnivMesozoo_LO_Cdt_prd["p4"] + domnivMesozoo_LO_Cdt_prd["z3"] + domnivMesozoo_LO_Cdt_prd["z5"] + domnivMesozoo_LO_Cdt_prd["z6"] - dcarnivMesozoo_LO_Cdt_prd["z4"] - domnivMesozoo_LO_Cdt_rel_labileDOM_NO_C - domnivMesozoo_LO_Cdt_rel_particOrganDetritus_NO_C - domnivMesozoo_LO_Cdt_rsp_disInorgCarbon_IO_C
    domnivMesozoo_LO_N_dt = diatoms_LO_N_diatoms_LO_C*domnivMesozoo_LO_Cdt_prd["p1"] + nanoflagellates_LO_N_nanoflagellates_LO_C*domnivMesozoo_LO_Cdt_prd["p2"] + picophyto_LO_N_picophyto_LO_C*domnivMesozoo_LO_Cdt_prd["p3"] + largephyto_LO_N_largephyto_LO_C*domnivMesozoo_LO_Cdt_prd["p4"] + carnivMesozoo_LO_N_carnivMesozoo_LO_C*domnivMesozoo_LO_Cdt_prd["z3"] + microzoo_LO_N_microzoo_LO_C*domnivMesozoo_LO_Cdt_prd["z5"] + heteroFlagellates_LO_N_heteroFlagellates_LO_C*domnivMesozoo_LO_Cdt_prd["z6"] - omnivMesozoo_LO_N_omnivMesozoo_LO_C*dcarnivMesozoo_LO_Cdt_prd["z4"] - domnivMesozoo_LO_Ndt_rel_labileDOM_NO_N - domnivMesozoo_LO_Ndt_rel_particOrganDetritus_NO_N - domnivMesozoo_LO_Ndt_rel_ammonium_IO_N
    domnivMesozoo_LO_P_dt = diatoms_LO_P_diatoms_LO_C*domnivMesozoo_LO_Cdt_prd["p1"] + nanoflagellates_LO_P_nanoflagellates_LO_C*domnivMesozoo_LO_Cdt_prd["p2"] + picophyto_LO_P_picophyto_LO_C*domnivMesozoo_LO_Cdt_prd["p3"] + largephyto_LO_P_largephyto_LO_C*domnivMesozoo_LO_Cdt_prd["p4"] + carnivMesozoo_LO_P_carnivMesozoo_LO_C*domnivMesozoo_LO_Cdt_prd["z3"] + microzoo_LO_P_microzoo_LO_C*domnivMesozoo_LO_Cdt_prd["z5"] + heteroFlagellates_LO_P_heteroFlagellates_LO_C*domnivMesozoo_LO_Cdt_prd["z6"] - omnivMesozoo_LO_P_omnivMesozoo_LO_C*dcarnivMesozoo_LO_Cdt_prd["z4"] - domnivMesozoo_LO_Pdt_rel_labileDOM_NO_P - domnivMesozoo_LO_Pdt_rel_particOrganDetritus_NO_P - domnivMesozoo_LO_Pdt_rel_phospate_IO_P
    
    # microzooplankton
    dmicrozoo_LO_C_dt = dmicrozoo_LO_Cdt_prd["b1"] + dmicrozoo_LO_Cdt_prd["p1"] + dmicrozoo_LO_Cdt_prd["p2"] + dmicrozoo_LO_Cdt_prd["p3"] + dmicrozoo_LO_Cdt_prd["p4"] + dmicrozoo_LO_Cdt_prd["z6"] - dcarnivMesozoo_LO_Cdt_prd["z5"] - domnivMesozoo_LO_Cdt_prd["z5"] - dheteroFlagellates_LO_Cdt_prd["z5"] - dmicrozoo_LO_Cdt_rel_labileDOM_NO_C - dmicrozoo_LO_Cdt_rel_particOrganDetritus_NO_C - dmicrozoo_LO_Cdt_rsp_disInorgCarbon_IO_C
    dmicrozoo_LO_N_dt = bn_bc*dmicrozoo_LO_Cdt_prd["b1"] + diatoms_LO_N_diatoms_LO_C*dmicrozoo_LO_Cdt_prd["p1"] + nanoflagellates_LO_N_nanoflagellates_LO_C*dmicrozoo_LO_Cdt_prd["p2"] + picophyto_LO_N_picophyto_LO_C*dmicrozoo_LO_Cdt_prd["p3"] + largephyto_LO_N_largephyto_LO_C*dmicrozoo_LO_Cdt_prd["p4"] + heteroFlagellates_LO_N_heteroFlagellates_LO_C*dmicrozoo_LO_Cdt_prd["z6"] - microzoo_LO_N_microzoo_LO_C*dcarnivMesozoo_LO_Cdt_prd["z5"] - microzoo_LO_N_microzoo_LO_C*domnivMesozoo_LO_Cdt_prd["z5"] - microzoo_LO_N_microzoo_LO_C*dheteroFlagellates_LO_Cdt_prd["z5"] - dmicrozoo_LO_Ndt_rel_labileDOM_NO_N - dmicrozoo_LO_Ndt_rel_particOrganDetritus_NO_N - dmicrozoo_LO_Ndt_rel_ammonium_IO_N
    dmicrozoo_LO_P_dt = bp_bc*dmicrozoo_LO_Cdt_prd["b1"] + diatoms_LO_P_diatoms_LO_C*dmicrozoo_LO_Cdt_prd["p1"] + nanoflagellates_LO_P_nanoflagellates_LO_C*dmicrozoo_LO_Cdt_prd["p2"] + picophyto_LO_P_picophyto_LO_C*dmicrozoo_LO_Cdt_prd["p3"] + largephyto_LO_P_largephyto_LO_C*dmicrozoo_LO_Cdt_prd["p4"] + heteroFlagellates_LO_P_heteroFlagellates_LO_C*dmicrozoo_LO_Cdt_prd["z6"] - microzoo_LO_P_microzoo_LO_C*dcarnivMesozoo_LO_Cdt_prd["z5"] - microzoo_LO_P_microzoo_LO_C*domnivMesozoo_LO_Cdt_prd["z5"] - microzoo_LO_P_microzoo_LO_C*dheteroFlagellates_LO_Cdt_prd["z5"] - dmicrozoo_LO_Pdt_rel_labileDOM_NO_P - dmicrozoo_LO_Pdt_rel_particOrganDetritus_NO_P - dmicrozoo_LO_Pdt_rel_phospate_IO_P
    
    dheteroFlagellates_LO_C_dt = dheteroFlagellates_LO_Cdt_prd["b1"] + dheteroFlagellates_LO_Cdt_prd["p1"] + dheteroFlagellates_LO_Cdt_prd["p2"] + dheteroFlagellates_LO_Cdt_prd["p3"] + dheteroFlagellates_LO_Cdt_prd["p4"] + dheteroFlagellates_LO_Cdt_prd["z5"] - dcarnivMesozoo_LO_Cdt_prd["z6"] - domnivMesozoo_LO_Cdt_prd["z6"] - dmicrozoo_LO_Cdt_prd["z6"] - dheteroFlagellates_LO_Cdt_rel_labileDOM_NO_C - dheteroFlagellates_LO_Cdt_rel_particOrganDetritus_NO_C - dheteroFlagellates_LO_Cdt_rsp_disInorgCarbon_IO_C
    dheteroFlagellates_LO_N_dt = bn_bc*dheteroFlagellates_LO_Cdt_prd["b1"] + diatoms_LO_N_diatoms_LO_C*dheteroFlagellates_LO_Cdt_prd["p1"] + nanoflagellates_LO_N_nanoflagellates_LO_C*dheteroFlagellates_LO_Cdt_prd["p2"] + picophyto_LO_N_picophyto_LO_C*dheteroFlagellates_LO_Cdt_prd["p3"] + largephyto_LO_N_largephyto_LO_C*dheteroFlagellates_LO_Cdt_prd["p4"] + microzoo_LO_N_microzoo_LO_C*dheteroFlagellates_LO_Cdt_prd["z5"] - heteroFlagellates_LO_N_heteroFlagellates_LO_C*dcarnivMesozoo_LO_Cdt_prd["z6"] - heteroFlagellates_LO_N_heteroFlagellates_LO_C*domnivMesozoo_LO_Cdt_prd["z6"] - heteroFlagellates_LO_N_heteroFlagellates_LO_C*dmicrozoo_LO_Cdt_prd["z6"] - dheteroFlagellates_LO_Ndt_rel_labileDOM_NO_N - dheteroFlagellates_LO_Ndt_rel_particOrganDetritus_NO_N - dheteroFlagellates_LO_Ndt_rel_ammonium_IO_N
    dheteroFlagellates_LO_P_dt = bp_bc*dheteroFlagellates_LO_Cdt_prd["b1"] + diatoms_LO_P_diatoms_LO_C*dheteroFlagellates_LO_Cdt_prd["p1"] + nanoflagellates_LO_P_nanoflagellates_LO_C*dheteroFlagellates_LO_Cdt_prd["p2"] + picophyto_LO_P_picophyto_LO_C*dheteroFlagellates_LO_Cdt_prd["p3"] + largephyto_LO_P_largephyto_LO_C*dheteroFlagellates_LO_Cdt_prd["p4"] + microzoo_LO_P_microzoo_LO_C*dheteroFlagellates_LO_Cdt_prd["z5"] - heteroFlagellates_LO_P_heteroFlagellates_LO_C*dcarnivMesozoo_LO_Cdt_prd["z6"] - heteroFlagellates_LO_P_heteroFlagellates_LO_C*domnivMesozoo_LO_Cdt_prd["z6"] - heteroFlagellates_LO_P_heteroFlagellates_LO_C*dmicrozoo_LO_Cdt_prd["z6"] - dheteroFlagellates_LO_Pdt_rel_labileDOM_NO_P - dheteroFlagellates_LO_Pdt_rel_particOrganDetritus_NO_P - dheteroFlagellates_LO_Pdt_rel_phospate_IO_P

    # DOM
    dlabileDOM_NO_C_dt = (ddiatoms_LO_Cdt_lys_labileDOM_NO_C + dnanoflagellates_LO_Cdt_lys_labileDOM_NO_C + dpicophyto_LO_Cdt_lys_labileDOM_NO_C + dlargephyto_LO_Cdt_lys_labileDOM_NO_C) + dBcdt_lys_labileDOM_NO_C - dBcdt_upt_labileDOM_NO_C + (dmicrozoo_LO_Cdt_rel_labileDOM_NO_C + dheteroFlagellates_LO_Cdt_rel_labileDOM_NO_C)
    dlabileDOM_NO_N_dt = (ddiatoms_LO_Ndt_lys_labileDOM_NO_N + dnanoflagellates_LO_Ndt_lys_labileDOM_NO_N + dpicophyto_LO_Ndt_lys_labileDOM_NO_N + dlargephyto_LO_Ndt_lys_labileDOM_NO_N) + (extra_n1 + extra_n2 + extra_n3 + extra_n4) + dBcdt_lys_labileDOM_NO_N - dBcdt_upt_labileDOM_NO_C*labileDOM_NO_N_labileDOM_NO_C + (dmicrozoo_LO_Ndt_rel_labileDOM_NO_N + dheteroFlagellates_LO_Ndt_rel_labileDOM_NO_N)
    dlabileDOM_NO_P_dt = (ddiatoms_LO_Pdt_lys_labileDOM_NO_P + dnanoflagellates_LO_Pdt_lys_labileDOM_NO_P + dpicophyto_LO_Pdt_lys_labileDOM_NO_P + dlargephyto_LO_Pdt_lys_labileDOM_NO_P) + (ddiatoms_LO_Pdt_upt_labileDOM_NO_P + dnanoflagellates_LO_Pdt_upt_labileDOM_NO_P + dpicophyto_LO_Pdt_upt_labileDOM_NO_P + dlargephyto_LO_Pdt_upt_labileDOM_NO_P) + dBcdt_lys_labileDOM_NO_P - dBcdt_upt_labileDOM_NO_C*labileDOM_NO_P_labileDOM_NO_C + (dmicrozoo_LO_Pdt_rel_labileDOM_NO_P + dheteroFlagellates_LO_Pdt_rel_labileDOM_NO_P)
    dsemilabileDOC_NO_C_dt = (ddiatoms_LO_Cdt_exu_semilabileDOC_NO_C + dnanoflagellates_LO_Cdt_exu_semilabileDOC_NO_C + dpicophyto_LO_Cdt_exu_semilabileDOC_NO_C + dlargephyto_LO_Cdt_exu_semilabileDOC_NO_C) - dBcdt_upt_semilabileDOC_NO_C + dBcdt_rel_semilabileDOC_NO_C
    dsemirefractDOC_NO_C_dt = dBcdt_rel_semirefractDOC_NO_C - dBcdt_upt_semirefractDOC_NO_C

    # POM
    dparticOrganDetritus_NO_C_dt = (ddiatoms_LO_Cdt_lys_particOrganDetritus_NO_C + dnanoflagellates_LO_Cdt_lys_particOrganDetritus_NO_C + dpicophyto_LO_Cdt_lys_particOrganDetritus_NO_C + dlargephyto_LO_Cdt_lys_particOrganDetritus_NO_C) + dBcdt_lys_particOrganDetritus_NO_C - dBcdt_upt_particOrganDetritus_NO_C + (dcarnivMesozoo_LO_Cdt_rel_particOrganDetritus_NO_C + domnivMesozoo_LO_Cdt_rel_particOrganDetritus_NO_C + dmicrozoo_LO_Cdt_rel_particOrganDetritus_NO_C + dheteroFlagellates_LO_Cdt_rel_particOrganDetritus_NO_C)
    dparticOrganDetritus_NO_N_dt = (ddiatoms_LO_Ndt_lys_particOrganDetritus_NO_N + dnanoflagellates_LO_Ndt_lys_particOrganDetritus_NO_N + dpicophyto_LO_Ndt_lys_particOrganDetritus_NO_N + dlargephyto_LO_Ndt_lys_particOrganDetritus_NO_N) + dBcdt_lys_particOrganDetritus_NO_N - dBcdt_upt_particOrganDetritus_NO_C*particOrganDetritus_NO_N_particOrganDetritus_NO_C + (dcarnivMesozoo_LO_Ndt_rel_particOrganDetritus_NO_N + domnivMesozoo_LO_Ndt_rel_particOrganDetritus_NO_N + dmicrozoo_LO_Ndt_rel_particOrganDetritus_NO_N + dheteroFlagellates_LO_Ndt_rel_particOrganDetritus_NO_N)
    dparticOrganDetritus_NO_P_dt = (ddiatoms_LO_Pdt_lys_particOrganDetritus_NO_P + dnanoflagellates_LO_Pdt_lys_particOrganDetritus_NO_P + dpicophyto_LO_Pdt_lys_particOrganDetritus_NO_P + dlargephyto_LO_Pdt_lys_particOrganDetritus_NO_P) + dBcdt_lys_particOrganDetritus_NO_P - dBcdt_upt_particOrganDetritus_NO_C*particOrganDetritus_NO_P_particOrganDetritus_NO_C + (dcarnivMesozoo_LO_Pdt_rel_particOrganDetritus_NO_P + domnivMesozoo_LO_Pdt_rel_particOrganDetritus_NO_P + dmicrozoo_LO_Pdt_rel_particOrganDetritus_NO_P + dheteroFlagellates_LO_Pdt_rel_particOrganDetritus_NO_P)
    dparticOrganDetritus_NO_Si_dt = ddiatoms_LO_Sidt_lys_particOrganDetritus_NO_Si - dparticOrganDetritus_NO_Sidt_rmn_silicate_IO_Si + (diatoms_LO_Si_diatoms_LO_C*(dcarnivMesozoo_LO_Cdt_prd["p1"] + domnivMesozoo_LO_Cdt_prd["p1"] + dmicrozoo_LO_Cdt_prd["p1"] + dheteroFlagellates_LO_Cdt_prd["p1"]))

    # Dissolved inorganic carbon
    # ddisInorgCarbon_IO_C_dt = (-ddiatoms_LO_Cdt_gpp_disInorgCarbon_IO_C + ddiatoms_LO_Cdt_rsp_disInorgCarbon_IO_C) + (-dnanoflagellates_LO_Cdt_gpp_disInorgCarbon_IO_C + dnanoflagellates_LO_Cdt_rsp_disInorgCarbon_IO_C) + (-dpicophyto_LO_Cdt_gpp_disInorgCarbon_IO_C + dpicophyto_LO_Cdt_rsp_disInorgCarbon_IO_C) + (-dlargephyto_LO_Cdt_gpp_disInorgCarbon_IO_C + dlargephyto_LO_Cdt_rsp_disInorgCarbon_IO_C) + dBcdt_rsp_disInorgCarbon_IO_C + dcarnivMesozoo_LO_Cdt_rsp_disInorgCarbon_IO_C + domnivMesozoo_LO_Cdt_rsp_disInorgCarbon_IO_C + dmicrozoo_LO_Cdt_rsp_disInorgCarbon_IO_C + dheteroFlagellates_LO_Cdt_rsp_disInorgCarbon_IO_C + do3cdt_air_sea_flux
    ddisInorgCarbon_IO_C_dt = (-ddiatoms_LO_Cdt_gpp_disInorgCarbon_IO_C + ddiatoms_LO_Cdt_rsp_disInorgCarbon_IO_C) + (-dnanoflagellates_LO_Cdt_gpp_disInorgCarbon_IO_C + dnanoflagellates_LO_Cdt_rsp_disInorgCarbon_IO_C) + (-dpicophyto_LO_Cdt_gpp_disInorgCarbon_IO_C + dpicophyto_LO_Cdt_rsp_disInorgCarbon_IO_C) + (-dlargephyto_LO_Cdt_gpp_disInorgCarbon_IO_C + dlargephyto_LO_Cdt_rsp_disInorgCarbon_IO_C) + dBcdt_rsp_disInorgCarbon_IO_C + dcarnivMesozoo_LO_Cdt_rsp_disInorgCarbon_IO_C + domnivMesozoo_LO_Cdt_rsp_disInorgCarbon_IO_C + dmicrozoo_LO_Cdt_rsp_disInorgCarbon_IO_C + dheteroFlagellates_LO_Cdt_rsp_disInorgCarbon_IO_C

    # Total alkalinity (from Alkalinity.F90)
    if pel_chem_parameters["calc_alkalinity"] and disInorgCarbon_IO_C>0.0:
        dtotalAlkalinity_IO_dt = -dnitrate_IO_N_dt + dammonium_IO_N_dt
    else:
        dtotalAlkalinity_IO_dt = 0.0

    rates = np.array([ddisOxygen_IO_O_dt, dphospate_IO_P_dt, dnitrate_IO_N_dt, dammonium_IO_N_dt, dnitrogenSink_dt, dsilicate_IO_Si_dt, dreductEquiv_IO_R_dt, dpelBacteria_LO_C_dt, dpelBacteria_LO_N_dt, dpelBacteria_LO_P_dt, 
            ddiatoms_LO_C_dt, ddiatoms_LO_N_dt, ddiatoms_LO_P_dt, ddiatoms_LO_Chl_dt, ddiatoms_LO_Si_dt, dnanoflagellates_LO_C_dt, dnanoflagellates_LO_N_dt, dnanoflagellates_LO_P_dt, dnanoflagellates_LO_Chl_dt, 
            dpicophyto_LO_C_dt, dpicophyto_LO_N_dt, dpicophyto_LO_P_dt, dpicophyto_LO_Chl_dt, dlargephyto_LO_C_dt, dlargephyto_LO_N_dt, dlargephyto_LO_P_dt, dlargephyto_LO_Chl_dt, dcarnivMesozoo_LO_C_dt, dcarnivMesozoo_LO_N_dt, dcarnivMesozoo_LO_P_dt,
            domnivMesozoo_LO_C_dt, domnivMesozoo_LO_N_dt, domnivMesozoo_LO_P_dt, dmicrozoo_LO_C_dt, dmicrozoo_LO_N_dt, dmicrozoo_LO_P_dt, dheteroFlagellates_LO_C_dt, dheteroFlagellates_LO_N_dt, dheteroFlagellates_LO_P_dt, dlabileDOM_NO_C_dt, dlabileDOM_NO_N_dt, dlabileDOM_NO_P_dt, 
            dsemilabileDOC_NO_C_dt, dsemirefractDOC_NO_C_dt, dparticOrganDetritus_NO_C_dt, dparticOrganDetritus_NO_N_dt, dparticOrganDetritus_NO_P_dt, dparticOrganDetritus_NO_Si_dt, ddisInorgCarbon_IO_C_dt, dtotalAlkalinity_IO_dt])
    
    return rates, dOdt_wind, do3cdt_air_sea_flux

if __name__ == '__main__':
    # Names of species in the system
    species_names = ['disOxygen_IO_O', 'phospate_IO_P', 'nitrate_IO_N', 'ammonium_IO_N', 'nitrogenSink', 'silicate_IO_Si', 'reductEquiv_IO_R', 'pelBacteria_LO_C', 'pelBacteria_LO_N', 'pelBacteria_LO_P',
                     'diatoms_LO_C', 'diatoms_LO_N', 'diatoms_LO_P', 'diatoms_LO_Chl', 'diatoms_LO_Si', 'nanoflagellates_LO_C', 'nanoflagellates_LO_N', 'nanoflagellates_LO_P', 'nanoflagellates_LO_Chl',
                     'picophyto_LO_C', 'picophyto_LO_N', 'picophyto_LO_P', 'picophyto_LO_Chl', 'largephyto_LO_C', 'largephyto_LO_N', 'largephyto_LO_P', 'largephyto_LO_Chl',
                     'carnivMesozoo_LO_C', 'carnivMesozoo_LO_N', 'carnivMesozoo_LO_P', 'omnivMesozoo_LO_C', 'omnivMesozoo_LO_N', 'omnivMesozoo_LO_P', 'microzoo_LO_C', 'microzoo_LO_N', 'microzoo_LO_P',
                     'heteroFlagellates_LO_C', 'heteroFlagellates_LO_N', 'heteroFlagellates_LO_P', 'labileDOM_NO_C', 'labileDOM_NO_N', 'labileDOM_NO_P', 'semilabileDOC_NO_C', 'semirefractDOC_NO_C', 'particOrganDetritus_NO_C',
                     'particOrganDetritus_NO_N', 'particOrganDetritus_NO_P', 'particOrganDetritus_NO_Si', 'disInorgCarbon_IO_C', 'totalAlkalinity_IO']

    # Initial concentrations
    c0 = [300.0,                    # disOxygen_IO_O
          1.0,                      # phospate_IO_P
          5.0,                      # nitrate_IO_N
          1.0,                      # ammonium_IO_N
          200.0,                    # nitrogenSink
          8.0,                      # silicate_IO_Si
          1.0,                      # reductEquiv_IO_R
          1.0,                      # pelBacteria_LO_C
          1.67e-2,                  # pelBacteria_LO_N
          1.85e-3,                  # pelBacteria_LO_P
          1.0,                      # diatoms_LO_C
          1.26e-2,                  # diatoms_LO_N
          7.86e-4,                  # diatoms_LO_P
          2.50e-2,                  # diatoms_LO_Chl
          1.00e-2,                  # diatoms_LO_Si
          1.0,                      # nanoflagellates_LO_C
          1.26e-2,                  # nanoflagellates_LO_N
          7.86e-4,                  # nanoflagellates_LO_P
          1.50e-2,                  # nanoflagellates_LO_Chl
          1.0,                      # picophyto_LO_C
          1.26e-2,                  # picophyto_LO_N
          7.86e-4,                  # picophyto_LO_P
          2.00e-2,                  # picophyto_LO_Chl
          1.0,                      # largephyto_LO_C
          1.26e-2,                  # largephyto_LO_N
          7.86e-4,                  # largephyto_LO_P
          2.00e-2,                  # largephyto_LO_Chl
          1.0,                      # carnivMesozoo_LO_C
          1.5e-2,                   # carnivMesozoo_LO_N
          1.67e-3,                  # carnivMesozoo_LO_P
          1.0,                      # omnivMesozoo_LO_C
          1.5e-2,                   # omnivMesozoo_LO_N
          1.67e-3,                  # omnivMesozoo_LO_P
          1.0,                      # microzoo_LO_C
          1.67e-2,                  # microzoo_LO_N
          1.85e-3,                  # microzoo_LO_P
          1.0,                      # heteroFlagellates_LO_C
          1.67e-2,                  # heteroFlagellates_LO_N
          1.85e-3,                  # heteroFlagellates_LO_P
          1.0,                      # labileDOM_NO_C
          1.26e-2,                  # labileDOM_NO_N
          7.862e-4,                 # labileDOM_NO_P
          0.1,                      # semilabileDOC_NO_C
          1.0,                      # semirefractDOC_NO_C
          1.0,                      # particOrganDetritus_NO_C
          1.26e-2,                  # particOrganDetritus_NO_N
          7.862e-4,                 # particOrganDetritus_NO_P
          1.45e-2,                  # particOrganDetritus_NO_Si
          27060.0,                  # disInorgCarbon_IO_C
          2660.0                    # totalAlkalinity_IO
          ]

    # Time span for integration
    t_span = [0, 86400*365*10]

    # Assign multiplier = 1, indicating all species are present
    multiplier = np.ones(len(c0))

    # Integrate
    solution_full_model = solve_ivp(bfm50_rate_eqns, t_span, c0, method='RK23')
    
    # -------------------------------- Plots --------------------------------------
    # Plot style
    plt.rc('font', family='serif', size=14)
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)
    plt.rc('axes', labelsize=12, linewidth=2)
    
    species_names_units = ['Oxygen (mmol O$_2$/m$^3$)', 
                           'Phosphate (mmol P/m$^3$)',
                           'Nitrate (mmol N/m$^3$)', 
                           'Ammonium (mmol N/m$^3$)',
                           'Nitrogen sink (mmol N/m$^3$)', 
                           'Silicate (mmol Si/m$^3$)',
                           'Reduction Equivalents (mmol S/m$^3$)',
                           'Pelagic Bacteria (mg C/m$^3$)',
                           'Pelagic Bacteria (mmol N/m$^3$)',
                           'Pelagic Bacteria (mmol P/m$^3$)',
                           'Diatoms (mg C/m$^3$)',
                           'Diatoms (mmol N/m$^3$)',
                           'Diatoms (mmol P/m$^3$)',
                           'Diatoms (mg Chl-a/m$^3$)',
                           'Diatoms (mmol Si/m$^3$)',
                           'Flagellates (mg C/m$^3$)',
                           'Flagellates (mmol N/m$^3$)',
                           'Flagellates (mmol P/m$^3$)',
                           'Flagellates (mg Chl-a/m$^3$)',
                           'PicoPhytoplankton (mg C/m$^3$)',
                           'PicoPhytoplankton (mmol N/m$^3$)',
                           'PicoPhytoplankton (mmol P/m$^3$)',
                           'PicoPhytoplankton (mg Chl-a/m$^3$)',
                           'Large Phytoplankton (mg C/m$^3$)',
                           'Large Phytoplankton (mmol N/m$^3$)',
                           'Large Phytoplankton (mmol P/m$^3$)',
                           'Large Phytoplankton (mg Chl-a/m$^3$)',
                           'Carnivorous Mesozooplankton  (mg C/m$^3$)',
                           'Carnivorous Mesozooplankton  (mmol N/m$^3$)',
                           'Carnivorous Mesozooplankton  (mmol P/m$^3$)',
                           'Omnivorous Mesozooplankton (mg C/m$^3$)',
                           'Omnivorous Mesozooplankton (mmol N/m$^3$)',
                           'Omnivorous Mesozooplankton (mmol P/m$^3$)',
                           'Microzooplankton (mg C/m$^3$)',
                           'Microzooplankton (mmol N/m$^3$)',
                           'Microzooplankton (mmol P/m$^3$)',
                           'Heterotrophic Nanoflagellates (mg C/m$^3$)',
                           'Heterotrophic Nanoflagellates (mmol N/m$^3$)',
                           'Heterotrophic Nanoflagellates (mmol P/m$^3$)',
                           'Labile Dissolved Organic Carbon (mg C/m$^3$)',
                           'Labile Dissolved Organic Nitrogen (mmol N/m$^3$)',
                           'Labile Dissolved Organic Phosphate (mmol P/m$^3$)',
                           'Semi-labile Dissolved Organic Carbon (mg C/m$^3$)',
                           'Semi-refractory Dissolved Organic Carbon (mg C/m$^3$)',
                           'Particulate Organic Carbon (mg C/m$^3$)',
                           'Particulate Organic Nitrate (mmol N/m$^3$)',
                           'Particulate Organic Phosphate (mmol P/m$^3$)',
                           'Particulate Organic Silicate (mmol Si/m$^3$)',
                           'Dissolved Inorganic Carbon (mg C/m$^3$)',
                           'Total Alkalinity (mmol Eq/m$^3$)'
                           ]

    plt.close('all')
    with PdfPages('BFM50_plot.pdf') as pdf:
        for index,species in enumerate(species_names_units):
            fig = plt.figure(index)
            ax = fig.add_axes([0, 0, 1, 1])
            ax.grid(linestyle='--')
            ax.plot(solution_full_model.t/86400/365, solution_full_model.y[index])
            ax.set_ylabel(species)
            ax.set_xlabel('Time (year)')
            ax.set_xlim(0,10)
            if index == 4:
                plt.ylim([199.9,200.1])
            if index == 42:
                plt.ylim([0.09,0.11])
            pdf.savefig(fig, bbox_inches = "tight")
            if index == 19 or index == 39 or index == 49:
                plt.close('all')
    plt.close('all')
