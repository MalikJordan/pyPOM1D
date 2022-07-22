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
from pom.constants import current_path, num_boxes


# Names of species in the system
species_names = ['o2o', 'n1p', 'n3n', 'n4n', 'o4n', 'n5s', 'n6r', 'b1c', 'b1n', 'b1p', 
                 'p1c', 'p1n', 'p1p', 'p1l', 'p1s', 'p2c', 'p2n', 'p2p', 'p2l',
                 'p3c', 'p3n', 'p3p', 'p3l', 'p4c', 'p4n', 'p4p', 'p4l',
                 'z3c', 'z3n', 'z3p', 'z4c', 'z4n', 'z4p', 'z5c', 'z5n', 'z5p',
                 'z6c', 'z6n', 'z6p', 'r1c', 'r1n', 'r1p', 'r2c', 'r3c', 'r6c', 
                 'r6n', 'r6p', 'r6s', 'o3c', 'o3h']


# def bfm50_rate_eqns(count, bfm_phys_vars, time, conc, seasonal_cycle=False):
def bfm50_rate_eqns(bfm_phys_vars, time, d3state, seasonal_cycle=False):
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

        # wind = bfm_phys_vars.wind
        # temper = bfm_phys_vars.temperature[count]
        # salt = bfm_phys_vars.salinity[count]
        # # qs = bfm_phys_vars.irradiance
        # xEPS = bfm_phys_vars.vertical_extinction[count,:]
        # irradiance = bfm_phys_vars.irradiance[count,:]
        # del_z = bfm_phys_vars.depth[count]
        # suspended_sediments = bfm_phys_vars.suspended_matter[count]
        # rho = bfm_phys_vars.density[count]

        wind = bfm_phys_vars.wind
        temper = bfm_phys_vars.temperature
        salt = bfm_phys_vars.salinity
        xEPS = bfm_phys_vars.vertical_extinction
        irradiance = bfm_phys_vars.irradiance
        del_z = bfm_phys_vars.depth
        suspended_sediments = bfm_phys_vars.suspended_matter
        rho = bfm_phys_vars.density


    #--------------------------------------------------------------------------
    # State variables
    o2o = d3state[:,0]              # Dissolved oxygen (mg O_2 m^-3)
    n1p = d3state[:,1]              # Phosphate (mmol P m^-3)
    n3n = d3state[:,2]              # Nitrate (mmol N m^-3)
    n4n = d3state[:,3]              # Ammonium (mmol N m^-3)
    o4n = d3state[:,4]              # Nitrogen sink (mmol N m^-3)
    n5s = d3state[:,5]              # Silicate (mmol Si m^-3)
    n6r = d3state[:,6]              # Reduction equivalents (mmol S m^-3)
    b1c = d3state[:,7]              # Pelagic bacteria carbon (mg C m^-3)
    b1n = d3state[:,8]              # Pelagic bacteria nitrogen (mmol N m^-3)
    b1p = d3state[:,9]              # Pelagic bacteria phosphate (mmol P m^-3)
    p1c = d3state[:,10]             # Diatoms carbon (mg C m^-3)
    p1n = d3state[:,11]             # Diatoms nitrogen (mmol N m^-3)
    p1p = d3state[:,12]             # Diatoms phosphate (mmol P m^-3)
    p1l = d3state[:,13]             # Diatoms chlorophyll (mg Chl-a m^-3)
    p1s = d3state[:,14]             # Diatoms silicate (mmol Si m^-3) 
    p2c = d3state[:,15]             # NanoFlagellates carbon (mg C m^-3)
    p2n = d3state[:,16]             # NanoFlagellates nitrogen (mmol N m^-3)
    p2p = d3state[:,17]             # NanoFlagellates phosphate (mmol P m^-3)
    p2l = d3state[:,18]             # NanoFlagellates chlorophyll (mg Chl-a m^-3)
    p3c = d3state[:,19]             # Picophytoplankton carbon (mg C m^-3)
    p3n = d3state[:,20]             # Picophytoplankton nitrogen (mmol N m^-3)
    p3p = d3state[:,21]             # Picophytoplankton phosphate (mmol P m^-3)
    p3l = d3state[:,22]             # Picophytoplankton chlorophyll (mg Chl-a m^-3)
    p4c = d3state[:,23]             # Large phytoplankton carbon (mg C m^-3)
    p4n = d3state[:,24]             # Large phytoplankton nitrogen (mmol N m^-3)
    p4p = d3state[:,25]             # Large phytoplankton phosphate (mmol P m^-3) 
    p4l = d3state[:,26]             # Large phytoplankton chlorophyll (mg Chl-a m^-3)
    z3c = d3state[:,27]             # Carnivorous mesozooplankton carbon (mg C m^-3)
    z3n = d3state[:,28]             # Carnivorous mesozooplankton nitrogen (mmol N m^-3)
    z3p = d3state[:,29]             # Carnivorous mesozooplankton phosphate (mmol P m^-3)
    z4c = d3state[:,30]             # Omnivorous mesozooplankton carbon (mg C m^-3)
    z4n = d3state[:,31]             # Omnivorous mesozooplankton nitrogen (mmol N m^-3)
    z4p = d3state[:,32]             # Omnivorous mesozooplankton phosphate (mmol P m^-3)
    z5c = d3state[:,33]             # Microzooplankton carbon (mg C m^-3)
    z5n = d3state[:,34]             # Microzooplankton nitrogen (mmol N m^-3)
    z5p = d3state[:,35]             # Microzooplankton phosphate (mmol P m^-3)
    z6c = d3state[:,36]             # Heterotrophic flagellates carbon (mg C m^-3)
    z6n = d3state[:,37]             # Heterotrophic flagellates nitrogen (mmol N m^-3)
    z6p = d3state[:,38]             # Heterotrophic flagellates phosphate (mmol P m^-3)
    r1c = d3state[:,39]             # Labile dissolved organic carbon (mg C m^-3)
    r1n = d3state[:,40]             # Labile dissolved organic nitrogen (mmol N m^-3)
    r1p = d3state[:,41]             # Labile dissolved organic phosphate (mmol P m^-3)
    r2c = d3state[:,42]             # Semi-labile dissolved organic carbon (mg C m^-3)
    r3c = d3state[:,43]             # Semi-refractory Dissolved Organic Carbon (mg C m^-3)
    r6c = d3state[:,44]             # Particulate organic carbon (mg C m^-3)
    r6n = d3state[:,45]             # Particulate organic nitrogen (mmol N m^-3)
    r6p = d3state[:,46]             # Particulate organic phosphate (mmol P m^-3)
    r6s = d3state[:,47]             # Particulate organic silicate (mmol Si m^-3)
    o3c = d3state[:,48]             # Dissolved inorganic carbon(mg C m^-3)
    o3h = d3state[:,49]             # Total alkalinity (mmol Eq m^-3)

    #--------------------------------------------------------------------------
    # concentration ratios
    p1n_p1c = get_concentration_ratio(p1n, p1c, constant_parameters["p_small"])
    p1p_p1c = get_concentration_ratio(p1p, p1c, constant_parameters["p_small"])
    p1l_p1c = get_concentration_ratio(p1l, p1c, constant_parameters["p_small"])
    p1s_p1c = get_concentration_ratio(p1s, p1c, constant_parameters["p_small"])
    p2n_p2c = get_concentration_ratio(p2n, p2c, constant_parameters["p_small"])
    p2p_p2c = get_concentration_ratio(p2p, p2c, constant_parameters["p_small"])
    p2l_p2c = get_concentration_ratio(p2l, p2c, constant_parameters["p_small"])
    p3n_p3c = get_concentration_ratio(p3n, p3c, constant_parameters["p_small"])
    p3p_p3c = get_concentration_ratio(p3p, p3c, constant_parameters["p_small"])
    p3l_p3c = get_concentration_ratio(p3l, p3c, constant_parameters["p_small"])
    p4n_p4c = get_concentration_ratio(p4n, p4c, constant_parameters["p_small"])
    p4p_p4c = get_concentration_ratio(p4p, p4c, constant_parameters["p_small"])
    p4l_p4c = get_concentration_ratio(p4l, p4c, constant_parameters["p_small"])
    bp_bc = get_concentration_ratio(b1p, b1c, constant_parameters["p_small"])
    bn_bc = get_concentration_ratio(b1n, b1c, constant_parameters["p_small"])
    z3n_z3c = get_concentration_ratio(z3n, z3c, constant_parameters["p_small"])
    z3p_z3c = get_concentration_ratio(z3p, z3c, constant_parameters["p_small"])
    z4n_z4c = get_concentration_ratio(z4n, z4c, constant_parameters["p_small"])
    z4p_z4c = get_concentration_ratio(z4p, z4c, constant_parameters["p_small"])
    z5n_z5c = get_concentration_ratio(z5n, z5c, constant_parameters["p_small"])
    z5p_z5c = get_concentration_ratio(z5p, z5c, constant_parameters["p_small"])
    z6n_z6c = get_concentration_ratio(z6n, z6c, constant_parameters["p_small"])
    z6p_z6c = get_concentration_ratio(z6p, z6c, constant_parameters["p_small"])
    r1p_r1c = get_concentration_ratio(r1p, r1c, constant_parameters["p_small"])
    r6p_r6c = get_concentration_ratio(r6p, r6c, constant_parameters["p_small"])
    r1n_r1c = get_concentration_ratio(r1n, r1c, constant_parameters["p_small"])
    r6n_r6c = get_concentration_ratio(r6n, r6c, constant_parameters["p_small"])

    #--------------------------------------------------------------------------
    #---------------------- Phytoplankton Equations ---------------------------
    #--------------------------------------------------------------------------
    # P1: Diatoms terms
    (dp1cdt_gpp_o3c, dp1cdt_rsp_o3c, dp1cdt_lys_r1c, dp1cdt_lys_r6c, dp1cdt_exu_r2c, dp1ndt_upt_n3n, dp1ndt_upt_n4n, 
     extra_n1, dp1ndt_lys_r1n, dp1ndt_lys_r6n, dp1pdt_upt_n1p, dp1pdt_upt_r1p, dp1pdt_lys_r1p, dp1pdt_lys_r6p, 
     dp1ldt_syn, dp1sdt_upt_n5s, dp1sdt_lys_r6s) = phyto_eqns(d3state, phyto1_prameters, environmental_parameters, constant_parameters, del_z, 1, irradiance, p1c, p1n, p1p, p1l, suspended_sediments, temper, time, xEPS)
    # (dp1cdt_gpp_o3c, dp1cdt_rsp_o3c, dp1cdt_lys_r1c, dp1cdt_lys_r6c, dp1cdt_exu_r2c, dp1ndt_upt_n3n, dp1ndt_upt_n4n, 
    #  extra_n1, dp1ndt_lys_r1n, dp1ndt_lys_r6n, dp1pdt_upt_n1p, dp1pdt_upt_r1p, dp1pdt_lys_r1p, dp1pdt_lys_r6p, 
    #  dp1ldt_syn, dp1sdt_upt_n5s, dp1sdt_lys_r6s) = (np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), 
    #                                                 np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), 
    #                                                 np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes))
    # P2: Flagellates terms
    (dp2cdt_gpp_o3c, dp2cdt_rsp_o3c, dp2cdt_lys_r1c, dp2cdt_lys_r6c, dp2cdt_exu_r2c, dp2ndt_upt_n3n, dp2ndt_upt_n4n, 
     extra_n2, dp2ndt_lys_r1n, dp2ndt_lys_r6n, dp2pdt_upt_n1p, dp2pdt_upt_r1p, dp2pdt_lys_r1p, dp2pdt_lys_r6p, 
     dp2ldt_syn, dP2sdt_upt_n5s, dP2sdt_lys_r6s) = phyto_eqns(d3state, phyto2_prameters, environmental_parameters, constant_parameters, del_z, 2, irradiance, p2c, p2n, p2p, p2l, suspended_sediments, temper, time, xEPS)
    np.zeros(num_boxes)
    # P3: PicoPhytoplankton terms
    (dp3cdt_gpp_o3c, dp3cdt_rsp_o3c, dp3cdt_lys_r1c, dp3cdt_lys_r6c, dp3cdt_exu_r2c, dp3ndt_upt_n3n, dp3ndt_upt_n4n, 
     extra_n3, dp3ndt_lys_r1n, dp3ndt_lys_r6n, dp3pdt_upt_n1p, dp3pdt_upt_r1p, dp3pdt_lys_r1p, dp3pdt_lys_r6p, 
     dp3ldt_syn, dP3sdt_upt_n5s, dP3sdt_lys_r6s) = phyto_eqns(d3state, phyto3_prameters, environmental_parameters, constant_parameters, del_z, 3, irradiance, p3c, p3n, p3p, p3l, suspended_sediments, temper, time, xEPS)
    # (dp3cdt_gpp_o3c, dp3cdt_rsp_o3c, dp3cdt_lys_r1c, dp3cdt_lys_r6c, dp3cdt_exu_r2c, dp3ndt_upt_n3n, dp3ndt_upt_n4n, 
    #  extra_n3, dp3ndt_lys_r1n, dp3ndt_lys_r6n, dp3pdt_upt_n1p, dp3pdt_upt_r1p, dp3pdt_lys_r1p, dp3pdt_lys_r6p, 
    #  dp3ldt_syn, dP3sdt_upt_n5s, dP3sdt_lys_r6s) = (np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), 
    #                                                 np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), 
    #                                                 np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes))
    # P4: Large Phytoplankton terms
    (dp4cdt_gpp_o3c, dp4cdt_rsp_o3c, dp4cdt_lys_r1c, dp4cdt_lys_r6c, dp4cdt_exu_r2c, dp4ndt_upt_n3n, dp4ndt_upt_n4n, 
     extra_n4, dp4ndt_lys_r1n, dp4ndt_lys_r6n, dp4pdt_upt_n1p, dp4pdt_upt_r1p, dp4pdt_lys_r1p, dp4pdt_lys_r6p, 
     dp4ldt_syn, dP4sdt_upt_n5s, dP4sdt_lys_r6s) = phyto_eqns(d3state, phyto4_prameters, environmental_parameters, constant_parameters, del_z, 4, irradiance, p4c, p4n, p4p, p4l, suspended_sediments, temper, time, xEPS)
    # (dp4cdt_gpp_o3c, dp4cdt_rsp_o3c, dp4cdt_lys_r1c, dp4cdt_lys_r6c, dp4cdt_exu_r2c, dp4ndt_upt_n3n, dp4ndt_upt_n4n, 
    #  extra_n4, dp4ndt_lys_r1n, dp4ndt_lys_r6n, dp4pdt_upt_n1p, dp4pdt_upt_r1p, dp4pdt_lys_r1p, dp4pdt_lys_r6p, 
    #  dp4ldt_syn, dP4sdt_upt_n5s, dP4sdt_lys_r6s) = (np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), 
    #                                                 np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), 
    #                                                 np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes))
    #--------------------------------------------------------------------------
    #------------------------- Bacteria Equations -----------------------------
    #--------------------------------------------------------------------------
    (dBcdt_lys_r1c, dBcdt_lys_r1n, dBcdt_lys_r1p, dBcdt_lys_r6c, dBcdt_lys_r6n, dBcdt_lys_r6p, 
     dBcdt_upt_r1c, dBcdt_upt_r6c, dBpdt_upt_rel_n1p, dBndt_upt_rel_n4n, dBcdt_upt_r2c, dBcdt_upt_r3c, 
     dBcdt_rel_r2c, dBcdt_rel_r3c, dBcdt_rsp_o3c, flPTn6r, f_B_O, f_B_n, f_B_p) = bacteria_eqns(d3state, bacteria_parameters, constant_parameters, environmental_parameters, temper)
    # (dBcdt_lys_r1c, dBcdt_lys_r1n, dBcdt_lys_r1p, dBcdt_lys_r6c, dBcdt_lys_r6n, dBcdt_lys_r6p, 
    #  dBcdt_upt_r1c, dBcdt_upt_r6c, dBpdt_upt_rel_n1p, dBndt_upt_rel_n4n, dBcdt_upt_r2c, dBcdt_upt_r3c, 
    #  dBcdt_rel_r2c, dBcdt_rel_r3c, dBcdt_rsp_o3c, flPTn6r, f_B_O, f_B_n, f_B_p) = (np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), 
    #                                                                                np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), 
    #                                                                                np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes))
    #--------------------------------------------------------------------------
    #----------------------- Zooplankton Equations ----------------------------
    #--------------------------------------------------------------------------
    
    # Mesozooplankton predation terms
    dz3cdt_prd, dz4cdt_prd, ic3, in3, ip3, ic4, in4, ip4 = get_mesozoo_predation_terms(d3state, mesozoo3_parameters, mesozoo4_parameters, zoo_availability_parameters, environmental_parameters, constant_parameters, temper)

    # Microzooplankton predation terms
    dz5cdt_prd, dz6cdt_prd, ic5, in5, ip5, ic6, in6, ip6 = get_microzoo_predation_terms(d3state, microzoo5_parameters, microzoo6_parameters, zoo_availability_parameters, environmental_parameters, constant_parameters, temper)

    # Z3: Carnivorous Mesozooplankton terms
    (dz3cdt_rel_r1c, dz3cdt_rel_r6c, dz3cdt_rsp_o3c, dz3ndt_rel_r1n, dz3ndt_rel_r6n, dz3pdt_rel_r1p, dz3pdt_rel_r6p, 
     dz3pdt_rel_n1p, dz3ndt_rel_n4n) = mesozoo_eqns(d3state, mesozoo3_parameters, constant_parameters, environmental_parameters, z3c, z3n, z3p, ic3, in3, ip3, temper)
    # (dz3cdt_rel_r1c, dz3cdt_rel_r6c, dz3cdt_rsp_o3c, dz3ndt_rel_r1n, dz3ndt_rel_r6n, dz3pdt_rel_r1p, dz3pdt_rel_r6p, 
    #  dz3pdt_rel_n1p, dz3ndt_rel_n4n) = (np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), 
    #                                     np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes))
    # Z4: Omnivorous Mesozooplankton terms
    (dz4cdt_rel_r1c, dz4cdt_rel_r6c, dz4cdt_rsp_o3c, dz4ndt_rel_r1n, dz4ndt_rel_r6n, dz4pdt_rel_r1p, dz4pdt_rel_r6p, 
     dz4pdt_rel_n1p, dz4ndt_rel_n4n) = mesozoo_eqns(d3state, mesozoo4_parameters, constant_parameters, environmental_parameters, z4c, z4n, z4p, ic4, in4, ip4, temper)
    # (dz4cdt_rel_r1c, dz4cdt_rel_r6c, dz4cdt_rsp_o3c, dz4ndt_rel_r1n, dz4ndt_rel_r6n, dz4pdt_rel_r1p, dz4pdt_rel_r6p, 
    #  dz4pdt_rel_n1p, dz4ndt_rel_n4n) = (np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), 
    #                                     np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes))
    # Z5: Microzooplankton terms
    (dz5cdt_rel_r1c, dz5cdt_rel_r6c, dz5cdt_rsp_o3c, dz5ndt_rel_r1n, dz5ndt_rel_r6n, dz5pdt_rel_r1p, dz5pdt_rel_r6p, 
     dz5pdt_rel_n1p, dz5ndt_rel_n4n) = microzoo_eqns(d3state, microzoo5_parameters, constant_parameters, environmental_parameters, z5c, z5n, z5p, ic5, in5, ip5, temper)

    # Z6: Heterotrophic Nanoflagellates terms
    (dz6cdt_rel_r1c, dz6cdt_rel_r6c, dz6cdt_rsp_o3c, dz6ndt_rel_r1n, dz6ndt_rel_r6n, dz6pdt_rel_r1p, dz6pdt_rel_r6p,  
     dz6pdt_rel_n1p, dz6ndt_rel_n4n) = microzoo_eqns(d3state, microzoo6_parameters, constant_parameters, environmental_parameters, z6c, z6n, z6p, ic6, in6, ip6, temper)
    # (dz6cdt_rel_r1c, dz6cdt_rel_r6c, dz6cdt_rsp_o3c, dz6ndt_rel_r1n, dz6ndt_rel_r6n, dz6pdt_rel_r1p, dz6pdt_rel_r6p,  
    #  dz6pdt_rel_n1p, dz6ndt_rel_n4n) = (np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes), 
    #                                     np.zeros(num_boxes), np.zeros(num_boxes), np.zeros(num_boxes))
    #--------------------------------------------------------------------------
    #------------------------ Non-living components ---------------------------
    #--------------------------------------------------------------------------    
    (dn4ndt_nit_n3n, dn3ndt_denit, dn6rdt_reox, dr6sdt_rmn_n5s, dr6cdt_remin_o3c, dr1cdt_remin_o3c, dr6cdt_remin_o2o, dr1cdt_remin_o2o, 
            dr6pdt_remin_n1p, dr1pdt_remin_n1p, dr6ndt_remin_n4n, dr1ndt_remin_n4n) = pel_chem_eqns(pel_chem_parameters, environmental_parameters, constant_parameters, del_z, temper, d3state, flPTn6r)

    #---------------------- Oxygen airation by wind ---------------------------
    (dOdt_wind, jsurO2o) = calculate_oxygen_reaeration(oxygen_reaeration_parameters, environmental_parameters, constant_parameters, d3state, del_z, temper, salt, wind)

    #------------------------------- CO_2 Flux --------------------------------
    do3cdt_air_sea_flux = np.zeros(num_boxes)
    # for i in range(0,num_boxes):
    conc = d3state[0,:]
    temp = temper[0]
    sal = salt[0]
    dens = rho[0]
    dz = del_z[0]
    do3cdt_air_sea_flux[0] = calculate_co2_flux(co2_flux_parameters, environmental_parameters, constant_parameters, conc, temp, wind, sal, dens, dz)
        # do3cdt_air_sea_flux = calculate_co2_flux(co2_flux_parameters, environmental_parameters, constant_parameters, d3state, temper, wind, salt, rho, del_z)
    
    #--------------------------------------------------------------------------
    #----------------------------- Rate Equations -----------------------------
    #--------------------------------------------------------------------------
    
    # Dissolved oxygen [mmol O_2 m^-3 s^-1]
    do2o_dt = (constant_parameters["omega_c"]*((dp1cdt_gpp_o3c - dp1cdt_rsp_o3c) + (dp2cdt_gpp_o3c - dp2cdt_rsp_o3c) + (dp3cdt_gpp_o3c - dp3cdt_rsp_o3c) + (dp4cdt_gpp_o3c - dp4cdt_rsp_o3c)) - 
               constant_parameters["omega_c"]*f_B_O*dBcdt_rsp_o3c - 
               constant_parameters["omega_c"]*(dz3cdt_rsp_o3c + dz4cdt_rsp_o3c + dz5cdt_rsp_o3c + dz6cdt_rsp_o3c) - 
               constant_parameters["omega_n"]*dn4ndt_nit_n3n -
               (1.0/constant_parameters["omega_r"])*dn6rdt_reox - constant_parameters["omega_c"]*(dr6cdt_remin_o2o + dr1cdt_remin_o2o) +
               jsurO2o)

    # Dissolved inorganic nutrient equations
    dn1p_dt = - (dp1pdt_upt_n1p + dp2pdt_upt_n1p + dp3pdt_upt_n1p + dp4pdt_upt_n1p) + (dBpdt_upt_rel_n1p*insw_vector(dBpdt_upt_rel_n1p)) - ((-1)*f_B_p*dBpdt_upt_rel_n1p*insw_vector(-dBpdt_upt_rel_n1p)) + (dz3pdt_rel_n1p + dz4pdt_rel_n1p + dz5pdt_rel_n1p + dz6pdt_rel_n1p) + (dr6pdt_remin_n1p + dr1pdt_remin_n1p)
    dn3n_dt = - (dp1ndt_upt_n3n + dp2ndt_upt_n3n + dp3ndt_upt_n3n + dp4ndt_upt_n3n) + dn4ndt_nit_n3n - dn3ndt_denit
    dn4n_dt = - (dp1ndt_upt_n4n + dp2ndt_upt_n4n + dp3ndt_upt_n4n + dp4ndt_upt_n4n) + (dBndt_upt_rel_n4n*insw_vector(dBndt_upt_rel_n4n)) - ((-1)*f_B_n*dBndt_upt_rel_n4n*insw_vector(-dBndt_upt_rel_n4n)) + (dz3ndt_rel_n4n + dz4ndt_rel_n4n + dz5ndt_rel_n4n + dz6ndt_rel_n4n) - dn4ndt_nit_n3n + (dr6ndt_remin_n4n + dr1ndt_remin_n4n)
    do4n_dt = dn3ndt_denit
    dn5s_dt = - dp1sdt_upt_n5s + dr6sdt_rmn_n5s

    # Reduction equivalents
    dn6r_dt = constant_parameters["omega_r"]*constant_parameters["omega_c"]*(1.0 - f_B_O)*dBcdt_rsp_o3c - constant_parameters["omega_r"]*constant_parameters["omega_dn"]*dn3ndt_denit*insw_vector(-(o2o - n6r)/constant_parameters["omega_r"]) - dn6rdt_reox
    # dn6r_dt = np.zeros(num_boxes)

    # Bacterioplankton
    db1c_dt = dBcdt_upt_r1c + dBcdt_upt_r6c - dBcdt_rsp_o3c - dBcdt_lys_r1c - dBcdt_lys_r6c - (dz5cdt_prd["b1"] + dz6cdt_prd["b1"])
    db1n_dt = - dBcdt_lys_r1n - dBcdt_lys_r6n + (r1n_r1c*dBcdt_upt_r1c) + (r6n_r6c*dBcdt_upt_r6c) - (dBndt_upt_rel_n4n*insw_vector(dBndt_upt_rel_n4n)) + ((-1)*f_B_n*dBndt_upt_rel_n4n*insw_vector(-dBndt_upt_rel_n4n)) - (bn_bc*(dz5cdt_prd["b1"] + dz6cdt_prd["b1"]))
    db1p_dt = (r1p_r1c*dBcdt_upt_r1c) + (r6p_r6c*dBcdt_upt_r6c) - (dBpdt_upt_rel_n1p*insw_vector(dBpdt_upt_rel_n1p)) + ((-1)*f_B_p*dBpdt_upt_rel_n1p*insw_vector(-dBpdt_upt_rel_n1p)) - dBcdt_lys_r1p - dBcdt_lys_r6p - (bp_bc*(dz5cdt_prd["b1"] + dz6cdt_prd["b1"]))
    # db1c_dt = np.zeros(num_boxes)
    # db1n_dt = np.zeros(num_boxes)
    # db1p_dt = np.zeros(num_boxes)
    
    # Phytoplankton
    dp1c_dt = dp1cdt_gpp_o3c - dp1cdt_exu_r2c - dp1cdt_rsp_o3c - dp1cdt_lys_r1c - dp1cdt_lys_r6c - dz3cdt_prd["p1"] - dz4cdt_prd["p1"] - dz5cdt_prd["p1"] - dz6cdt_prd["p1"]
    dp1n_dt = dp1ndt_upt_n3n + dp1ndt_upt_n4n - extra_n1 - dp1ndt_lys_r1n - dp1ndt_lys_r6n - (p1n_p1c*(dz3cdt_prd["p1"] + dz4cdt_prd["p1"] + dz5cdt_prd["p1"] + dz6cdt_prd["p1"]))
    dp1p_dt = dp1pdt_upt_n1p - dp1pdt_upt_r1p - dp1pdt_lys_r1p - dp1pdt_lys_r6p - (p1p_p1c*(dz3cdt_prd["p1"] + dz4cdt_prd["p1"] + dz5cdt_prd["p1"] + dz6cdt_prd["p1"]))
    dp1l_dt = dp1ldt_syn - (p1l_p1c*(dz3cdt_prd["p1"] + dz4cdt_prd["p1"] + dz5cdt_prd["p1"] + dz6cdt_prd["p1"]))
    dp1s_dt = dp1sdt_upt_n5s - dp1sdt_lys_r6s - (p1s_p1c*(dz3cdt_prd["p1"] + dz4cdt_prd["p1"] + dz5cdt_prd["p1"] + dz6cdt_prd["p1"]))
    # dp1c_dt = np.zeros(num_boxes)
    # dp1n_dt = np.zeros(num_boxes)
    # dp1p_dt = np.zeros(num_boxes)
    # dp1l_dt = np.zeros(num_boxes)
    # dp1s_dt = np.zeros(num_boxes)
    
    dp2c_dt = dp2cdt_gpp_o3c - dp2cdt_exu_r2c - dp2cdt_rsp_o3c - dp2cdt_lys_r1c - dp2cdt_lys_r6c - dz3cdt_prd["p2"] - dz4cdt_prd["p2"] - dz5cdt_prd["p2"] - dz6cdt_prd["p2"]
    dp2n_dt = dp2ndt_upt_n3n + dp2ndt_upt_n4n - extra_n2 - dp2ndt_lys_r1n - dp2ndt_lys_r6n - (p2n_p2c*(dz3cdt_prd["p2"] + dz4cdt_prd["p2"] + dz5cdt_prd["p2"] + dz6cdt_prd["p2"]))
    dp2p_dt = dp2pdt_upt_n1p - dp2pdt_upt_r1p - dp2pdt_lys_r1p - dp2pdt_lys_r6p - (p2p_p2c*(dz3cdt_prd["p2"] + dz4cdt_prd["p2"] + dz5cdt_prd["p2"] + dz6cdt_prd["p2"]))
    dp2l_dt = dp2ldt_syn - (p2l_p2c*(dz3cdt_prd["p2"] + dz4cdt_prd["p2"] + dz5cdt_prd["p2"] + dz6cdt_prd["p2"]))
    
    dp3c_dt = dp3cdt_gpp_o3c - dp3cdt_exu_r2c - dp3cdt_rsp_o3c - dp3cdt_lys_r1c - dp3cdt_lys_r6c - dz3cdt_prd["p3"] - dz4cdt_prd["p3"] - dz5cdt_prd["p3"] - dz6cdt_prd["p3"]
    dp3n_dt = dp3ndt_upt_n3n + dp3ndt_upt_n4n - extra_n3 - dp3ndt_lys_r1n - dp3ndt_lys_r6n - (p3n_p3c*(dz3cdt_prd["p3"] + dz4cdt_prd["p3"] + dz5cdt_prd["p3"] + dz6cdt_prd["p3"]))
    dp3p_dt = dp3pdt_upt_n1p - dp3pdt_upt_r1p - dp3pdt_lys_r1p - dp3pdt_lys_r6p - (p3p_p3c*(dz3cdt_prd["p3"] + dz4cdt_prd["p3"] + dz5cdt_prd["p3"] + dz6cdt_prd["p3"]))
    dp3l_dt = dp3ldt_syn - (p3l_p3c*(dz3cdt_prd["p3"] + dz4cdt_prd["p3"] + dz5cdt_prd["p3"] + dz6cdt_prd["p3"]))
    # dp3c_dt = np.zeros(num_boxes)
    # dp3n_dt = np.zeros(num_boxes)
    # dp3p_dt = np.zeros(num_boxes)
    # dp3l_dt = np.zeros(num_boxes)
    
    dp4c_dt = dp4cdt_gpp_o3c - dp4cdt_exu_r2c - dp4cdt_rsp_o3c - dp4cdt_lys_r1c - dp4cdt_lys_r6c - dz3cdt_prd["p4"] - dz4cdt_prd["p4"] - dz5cdt_prd["p4"] - dz6cdt_prd["p4"]
    dp4n_dt = dp4ndt_upt_n3n + dp4ndt_upt_n4n - extra_n4 - dp4ndt_lys_r1n - dp4ndt_lys_r6n - (p4n_p4c*(dz3cdt_prd["p4"] + dz4cdt_prd["p4"] + dz5cdt_prd["p4"] + dz6cdt_prd["p4"]))
    dp4p_dt = dp4pdt_upt_n1p - dp4pdt_upt_r1p - dp4pdt_lys_r1p - dp4pdt_lys_r6p - (p4p_p4c*(dz3cdt_prd["p4"] + dz4cdt_prd["p4"] + dz5cdt_prd["p4"] + dz6cdt_prd["p4"]))
    dp4l_dt = dp4ldt_syn - (p4l_p4c*(dz3cdt_prd["p4"] + dz4cdt_prd["p4"] + dz5cdt_prd["p4"] + dz6cdt_prd["p4"]))
    # dp4c_dt = np.zeros(num_boxes)
    # dp4n_dt = np.zeros(num_boxes)
    # dp4p_dt = np.zeros(num_boxes)
    # dp4l_dt = np.zeros(num_boxes)
    
    # mesozooplankton
    dz3c_dt = dz3cdt_prd["p1"] + dz3cdt_prd["p2"] + dz3cdt_prd["p3"] + dz3cdt_prd["p4"] + dz3cdt_prd["z4"] + dz3cdt_prd["z5"] + dz3cdt_prd["z6"] - dz4cdt_prd["z3"] - dz3cdt_rel_r1c - dz3cdt_rel_r6c - dz3cdt_rsp_o3c
    dz3n_dt = p1n_p1c*dz3cdt_prd["p1"] + p2n_p2c*dz3cdt_prd["p2"] + p3n_p3c*dz3cdt_prd["p3"] + p4n_p4c*dz3cdt_prd["p4"] + z4n_z4c*dz3cdt_prd["z4"] + z5n_z5c*dz3cdt_prd["z5"] + z6n_z6c*dz3cdt_prd["z6"] - z3n_z3c*dz4cdt_prd["z3"] - dz3ndt_rel_r1n - dz3ndt_rel_r6n - dz3ndt_rel_n4n
    dz3p_dt = p1p_p1c*dz3cdt_prd["p1"] + p2p_p2c*dz3cdt_prd["p2"] + p3p_p3c*dz3cdt_prd["p3"] + p4p_p4c*dz3cdt_prd["p4"] + z4p_z4c*dz3cdt_prd["z4"] + z5p_z5c*dz3cdt_prd["z5"] + z6p_z6c*dz3cdt_prd["z6"] - z3p_z3c*dz4cdt_prd["z3"] - dz3pdt_rel_r1p - dz3pdt_rel_r6p - dz3pdt_rel_n1p
    # dz3c_dt = np.zeros(num_boxes)
    # dz3n_dt = np.zeros(num_boxes)
    # dz3p_dt = np.zeros(num_boxes)
    
    dz4c_dt = dz4cdt_prd["p1"] + dz4cdt_prd["p2"] + dz4cdt_prd["p3"] + dz4cdt_prd["p4"] + dz4cdt_prd["z3"] + dz4cdt_prd["z5"] + dz4cdt_prd["z6"] - dz3cdt_prd["z4"] - dz4cdt_rel_r1c - dz4cdt_rel_r6c - dz4cdt_rsp_o3c
    dz4n_dt = p1n_p1c*dz4cdt_prd["p1"] + p2n_p2c*dz4cdt_prd["p2"] + p3n_p3c*dz4cdt_prd["p3"] + p4n_p4c*dz4cdt_prd["p4"] + z3n_z3c*dz4cdt_prd["z3"] + z5n_z5c*dz4cdt_prd["z5"] + z6n_z6c*dz4cdt_prd["z6"] - z4n_z4c*dz3cdt_prd["z4"] - dz4ndt_rel_r1n - dz4ndt_rel_r6n - dz4ndt_rel_n4n
    dz4p_dt = p1p_p1c*dz4cdt_prd["p1"] + p2p_p2c*dz4cdt_prd["p2"] + p3p_p3c*dz4cdt_prd["p3"] + p4p_p4c*dz4cdt_prd["p4"] + z3p_z3c*dz4cdt_prd["z3"] + z5p_z5c*dz4cdt_prd["z5"] + z6p_z6c*dz4cdt_prd["z6"] - z4p_z4c*dz3cdt_prd["z4"] - dz4pdt_rel_r1p - dz4pdt_rel_r6p - dz4pdt_rel_n1p
    # dz4c_dt = np.zeros(num_boxes)
    # dz4n_dt = np.zeros(num_boxes)
    # dz4p_dt = np.zeros(num_boxes)
    
    # microzooplankton
    dz5c_dt = dz5cdt_prd["b1"] + dz5cdt_prd["p1"] + dz5cdt_prd["p2"] + dz5cdt_prd["p3"] + dz5cdt_prd["p4"] + dz5cdt_prd["z6"] - dz3cdt_prd["z5"] - dz4cdt_prd["z5"] - dz6cdt_prd["z5"] - dz5cdt_rel_r1c - dz5cdt_rel_r6c - dz5cdt_rsp_o3c
    dz5n_dt = bn_bc*dz5cdt_prd["b1"] + p1n_p1c*dz5cdt_prd["p1"] + p2n_p2c*dz5cdt_prd["p2"] + p3n_p3c*dz5cdt_prd["p3"] + p4n_p4c*dz5cdt_prd["p4"] + z6n_z6c*dz5cdt_prd["z6"] - z5n_z5c*dz3cdt_prd["z5"] - z5n_z5c*dz4cdt_prd["z5"] - z5n_z5c*dz6cdt_prd["z5"] - dz5ndt_rel_r1n - dz5ndt_rel_r6n - dz5ndt_rel_n4n
    dz5p_dt = bp_bc*dz5cdt_prd["b1"] + p1p_p1c*dz5cdt_prd["p1"] + p2p_p2c*dz5cdt_prd["p2"] + p3p_p3c*dz5cdt_prd["p3"] + p4p_p4c*dz5cdt_prd["p4"] + z6p_z6c*dz5cdt_prd["z6"] - z5p_z5c*dz3cdt_prd["z5"] - z5p_z5c*dz4cdt_prd["z5"] - z5p_z5c*dz6cdt_prd["z5"] - dz5pdt_rel_r1p - dz5pdt_rel_r6p - dz5pdt_rel_n1p
    
    dz6c_dt = dz6cdt_prd["b1"] + dz6cdt_prd["p1"] + dz6cdt_prd["p2"] + dz6cdt_prd["p3"] + dz6cdt_prd["p4"] + dz6cdt_prd["z5"] - dz3cdt_prd["z6"] - dz4cdt_prd["z6"] - dz5cdt_prd["z6"] - dz6cdt_rel_r1c - dz6cdt_rel_r6c - dz6cdt_rsp_o3c
    dz6n_dt = bn_bc*dz6cdt_prd["b1"] + p1n_p1c*dz6cdt_prd["p1"] + p2n_p2c*dz6cdt_prd["p2"] + p3n_p3c*dz6cdt_prd["p3"] + p4n_p4c*dz6cdt_prd["p4"] + z5n_z5c*dz6cdt_prd["z5"] - z6n_z6c*dz3cdt_prd["z6"] - z6n_z6c*dz4cdt_prd["z6"] - z6n_z6c*dz5cdt_prd["z6"] - dz6ndt_rel_r1n - dz6ndt_rel_r6n - dz6ndt_rel_n4n
    dz6p_dt = bp_bc*dz6cdt_prd["b1"] + p1p_p1c*dz6cdt_prd["p1"] + p2p_p2c*dz6cdt_prd["p2"] + p3p_p3c*dz6cdt_prd["p3"] + p4p_p4c*dz6cdt_prd["p4"] + z5p_z5c*dz6cdt_prd["z5"] - z6p_z6c*dz3cdt_prd["z6"] - z6p_z6c*dz4cdt_prd["z6"] - z6p_z6c*dz5cdt_prd["z6"] - dz6pdt_rel_r1p - dz6pdt_rel_r6p - dz6pdt_rel_n1p
    # dz6c_dt = np.zeros(num_boxes)
    # dz6n_dt = np.zeros(num_boxes)
    # dz6p_dt = np.zeros(num_boxes)
    
    # DOM
    dr1c_dt = (dp1cdt_lys_r1c + dp2cdt_lys_r1c + dp3cdt_lys_r1c + dp4cdt_lys_r1c) + dBcdt_lys_r1c - dBcdt_upt_r1c + (dz5cdt_rel_r1c + dz6cdt_rel_r1c) - dr1cdt_remin_o3c #  - dr1cdt_remin_o2o
    dr1n_dt = (dp1ndt_lys_r1n + dp2ndt_lys_r1n + dp3ndt_lys_r1n + dp4ndt_lys_r1n) + (extra_n1 + extra_n2 + extra_n3 + extra_n4) + dBcdt_lys_r1n - dBcdt_upt_r1c*r1n_r1c + (dz5ndt_rel_r1n + dz6ndt_rel_r1n) - dr1ndt_remin_n4n
    dr1p_dt = (dp1pdt_lys_r1p + dp2pdt_lys_r1p + dp3pdt_lys_r1p + dp4pdt_lys_r1p) + (dp1pdt_upt_r1p + dp2pdt_upt_r1p + dp3pdt_upt_r1p + dp4pdt_upt_r1p) + dBcdt_lys_r1p - dBcdt_upt_r1c*r1p_r1c + (dz5pdt_rel_r1p + dz6pdt_rel_r1p) - dr1pdt_remin_n1p
    dr2c_dt = (dp1cdt_exu_r2c + dp2cdt_exu_r2c + dp3cdt_exu_r2c + dp4cdt_exu_r2c) - dBcdt_upt_r2c + dBcdt_rel_r2c
    dr3c_dt = dBcdt_rel_r3c - dBcdt_upt_r3c
    # dr2c_dt = np.zeros(num_boxes)
    # dr3c_dt = np.zeros(num_boxes)
    
    # POM
    dr6c_dt = (dp1cdt_lys_r6c + dp2cdt_lys_r6c + dp3cdt_lys_r6c + dp4cdt_lys_r6c) + dBcdt_lys_r6c - dBcdt_upt_r6c + (dz3cdt_rel_r6c + dz4cdt_rel_r6c + dz5cdt_rel_r6c + dz6cdt_rel_r6c) - dr6cdt_remin_o3c #  - dr6cdt_remin_o2o
    dr6n_dt = (dp1ndt_lys_r6n + dp2ndt_lys_r6n + dp3ndt_lys_r6n + dp4ndt_lys_r6n) + dBcdt_lys_r6n - dBcdt_upt_r6c*r6n_r6c + (dz3ndt_rel_r6n + dz4ndt_rel_r6n + dz5ndt_rel_r6n + dz6ndt_rel_r6n) - dr6ndt_remin_n4n
    dr6p_dt = (dp1pdt_lys_r6p + dp2pdt_lys_r6p + dp3pdt_lys_r6p + dp4pdt_lys_r6p) + dBcdt_lys_r6p - dBcdt_upt_r6c*r6p_r6c + (dz3pdt_rel_r6p + dz4pdt_rel_r6p + dz5pdt_rel_r6p + dz6pdt_rel_r6p) - dr6pdt_remin_n1p
    dr6s_dt = dp1sdt_lys_r6s - dr6sdt_rmn_n5s + (p1s_p1c*(dz3cdt_prd["p1"] + dz4cdt_prd["p1"] + dz5cdt_prd["p1"] + dz6cdt_prd["p1"]))
    # dr6s_dt = np.zeros(num_boxes)
    
    # Dissolved inorganic carbon
    do3c_dt = (-dp1cdt_gpp_o3c + dp1cdt_rsp_o3c) + (-dp2cdt_gpp_o3c + dp2cdt_rsp_o3c) + (-dp3cdt_gpp_o3c + dp3cdt_rsp_o3c) + (-dp4cdt_gpp_o3c + dp4cdt_rsp_o3c) + dBcdt_rsp_o3c + dz3cdt_rsp_o3c + dz4cdt_rsp_o3c + dz5cdt_rsp_o3c + dz6cdt_rsp_o3c + (dr6cdt_remin_o3c + dr1cdt_remin_o3c) + do3cdt_air_sea_flux

    # Total alkalinity (from Alkalinity.F90)
    if pel_chem_parameters["calc_alkalinity"] and o3c>0.0:
        do3h_dt = -dn3n_dt + dn4n_dt
    else:
        do3h_dt = np.zeros(num_boxes)
        # do3h_dt = 0.0

    rates = np.array([do2o_dt, dn1p_dt, dn3n_dt, dn4n_dt, do4n_dt, dn5s_dt, dn6r_dt, db1c_dt, db1n_dt, db1p_dt, 
            dp1c_dt, dp1n_dt, dp1p_dt, dp1l_dt, dp1s_dt, dp2c_dt, dp2n_dt, dp2p_dt, dp2l_dt, 
            dp3c_dt, dp3n_dt, dp3p_dt, dp3l_dt, dp4c_dt, dp4n_dt, dp4p_dt, dp4l_dt, dz3c_dt, dz3n_dt, dz3p_dt,
            dz4c_dt, dz4n_dt, dz4p_dt, dz5c_dt, dz5n_dt, dz5p_dt, dz6c_dt, dz6n_dt, dz6p_dt, dr1c_dt, dr1n_dt, dr1p_dt, 
            dr2c_dt, dr3c_dt, dr6c_dt, dr6n_dt, dr6p_dt, dr6s_dt, do3c_dt, do3h_dt])/constant_parameters["sec_per_day"]
    
    return rates, dOdt_wind, do3cdt_air_sea_flux

if __name__ == '__main__':
    # Names of species in the system
    species_names = ['o2o', 'n1p', 'n3n', 'n4n', 'o4n', 'n5s', 'n6r', 'b1c', 'b1n', 'b1p',
                     'p1c', 'p1n', 'p1p', 'p1l', 'p1s', 'p2c', 'p2n', 'p2p', 'p2l',
                     'p3c', 'p3n', 'p3p', 'p3l', 'p4c', 'p4n', 'p4p', 'p4l',
                     'z3c', 'z3n', 'z3p', 'z4c', 'z4n', 'z4p', 'z5c', 'z5n', 'z5p',
                     'z6c', 'z6n', 'z6p', 'r1c', 'r1n', 'r1p', 'r2c', 'r3c', 'r6c',
                     'r6n', 'r6p', 'r6s', 'o3c', 'o3h']

    # Initial concentrations
    c0 = [300.0,                    # o2o
          1.0,                      # n1p
          5.0,                      # n3n
          1.0,                      # n4n
          200.0,                    # o4n
          8.0,                      # n5s
          1.0,                      # n6r
          1.0,                      # b1c
          1.67e-2,                  # b1n
          1.85e-3,                  # b1p
          1.0,                      # p1c
          1.26e-2,                  # p1n
          7.86e-4,                  # p1p
          2.50e-2,                  # p1l
          1.00e-2,                  # p1s
          1.0,                      # p2c
          1.26e-2,                  # p2n
          7.86e-4,                  # p2p
          1.50e-2,                  # p2l
          1.0,                      # p3c
          1.26e-2,                  # p3n
          7.86e-4,                  # p3p
          2.00e-2,                  # p3l
          1.0,                      # p4c
          1.26e-2,                  # p4n
          7.86e-4,                  # p4p
          2.00e-2,                  # p4l
          1.0,                      # z3c
          1.5e-2,                   # z3n
          1.67e-3,                  # z3p
          1.0,                      # z4c
          1.5e-2,                   # z4n
          1.67e-3,                  # z4p
          1.0,                      # z5c
          1.67e-2,                  # z5n
          1.85e-3,                  # z5p
          1.0,                      # z6c
          1.67e-2,                  # z6n
          1.85e-3,                  # z6p
          1.0,                      # r1c
          1.26e-2,                  # r1n
          7.862e-4,                 # r1p
          0.1,                      # r2c
          1.0,                      # r3c
          1.0,                      # r6c
          1.26e-2,                  # r6n
          7.862e-4,                 # r6p
          1.45e-2,                  # r6s
          27060.0,                  # o3c
          2660.0                    # o3h
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
