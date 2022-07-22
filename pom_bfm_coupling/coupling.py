import numpy as np
from pom.constants import vertical_layers, water_specific_heat_times_density, seconds_per_day
from inputs import params_POMBFM
from bfm.bfm50.BFM50_rate_eqns import bfm50_rate_eqns
from pom_bfm_coupling.calculations import calculate_vertical_diffusivity
from pom_bfm_coupling.initialize_variables import set_initial_conditions
from pom_bfm_coupling.calculate_average_field import calculate_daily_average_field
from pom.constants import current_path, num_boxes
import json

with open(current_path + "/bfm/bfm50/bfm50_parameters.json") as read_parameters:
    parameters = json.load(read_parameters)
constant_parameters = parameters["constants"]
env_parameters = parameters["environmental_parameters"]
phyto1_parameters = parameters["phyto1_parameters"]
phyto2_parameters = parameters["phyto2_parameters"]
phyto3_parameters = parameters["phyto3_parameters"]
phyto4_parameters = parameters["phyto4_parameters"]


def pom_to_bfm(bfm_phys_vars, vertical_grid, temperature, salinity, inorganic_suspended_matter, shortwave_radiation, vertical_density_profile, wind_stress):

    """
    Description: Passes the physical variables to the BFM

    :return: seawater density, temperature and salinity, suspended sediment load,
             photosynthetically available radiation, gridpoint depth, and wind speed
    """

    # try:
    #     import NOPOINTERS
    #     NOPOINTERS = True
    # except FileNotFoundError:
    #     NOPOINTERS = False
    # if NOPOINTERS:
    #     from modules import ETW, ESW, EIR, ESS, ERHO, EWIND, Depth

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   1D ARRAYS FOR BFM
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    bfm_phys_vars.temperature = temperature.backward[0:num_boxes]
    bfm_phys_vars.salinity = salinity.backward[0:num_boxes]
    bfm_phys_vars.density = (vertical_density_profile[0:num_boxes] * 1.E3) + 1.E3
    bfm_phys_vars.suspended_matter = inorganic_suspended_matter[0:num_boxes]
    bfm_phys_vars.depth = vertical_grid.vertical_spacing[0:num_boxes] * params_POMBFM.h

    bfm_phys_vars.irradiance[0,:] = -1. * shortwave_radiation * water_specific_heat_times_density

    wind = np.sqrt(wind_stress.zonal**2 + wind_stress.meridional**2) * 1.E3
    bfm_phys_vars.wind = np.sqrt(wind/(1.25 * 0.0014))

    return bfm_phys_vars


def pom_bfm_1d(i, vertical_grid, time, diffusion, nutrients, bfm_phys_vars, d3state, d3stateb, d3ave):

    # bfm_variables = set_initial_conditions()
    num_boxes = vertical_layers - 1
    bfm_rates = np.zeros((num_boxes,50))

    dOdt_wind = np.zeros(num_boxes)
    do3cdt_air_sea_flux = np.zeros(num_boxes)
    
    # for count in range(0,num_boxes):
    #     conc = d3state[count,:]
    #     bfm_rates[count,:], dOdt_wind[count], do3cdt_air_sea_flux[count] = bfm50_rate_eqns(count, bfm_phys_vars, time, conc, seasonal_cycle=False)
    
    bfm_rates, dOdt_wind, do3cdt_air_sea_flux = bfm50_rate_eqns(bfm_phys_vars, time, d3state, seasonal_cycle=False)

    d3state, d3stateb = calculate_vertical_diffusivity(vertical_grid, diffusion, nutrients, d3state, d3stateb, bfm_rates, bfm_phys_vars, dOdt_wind, do3cdt_air_sea_flux)
    d3ave.count += 1

    # if i == 0:
    #     TT = time - (params_POMBFM.dti/seconds_per_day)
    # TT = TT + (params_POMBFM.dti/seconds_per_day)

    # if d3ave.count == 0:
    #     d3ave = calculate_daily_average_field(d3ave,d3state,'Initialize')
    if (d3ave.count > 0) and (d3ave.count < seconds_per_day/params_POMBFM.dti):
        d3ave = calculate_daily_average_field(d3ave,d3state,'Accumulate')
    elif d3ave.count == seconds_per_day/params_POMBFM.dti:
        d3ave = calculate_daily_average_field(d3ave,d3state,'Mean')
        d3ave.count = 0
    elif d3ave.day == params_POMBFM.idays:
        d3ave = calculate_daily_average_field(d3ave,d3state,'Reset')

    return d3state, d3stateb, d3ave


def calculate_vertical_extinction(bfm_phys_vars, d3state, group):
    # Calcluations taken from bfm/bfm50/Functions/phyto.py - Updated from 0D to 1D

    # from CalcVerticalExtinction.F90 line 82
    bfm_phys_vars.vertical_extinction[:,group] = env_parameters["p_eps0"] + env_parameters["p_epsESS"]*bfm_phys_vars.suspended_matter + env_parameters["p_epsR6"]*d3state[:,44]
    # from CalcVerticalExtinction.F90 line 101 (ChlAttenFlag=1, ChlDynamicsFlag=2)       
    if group == 0: # P1: Diatoms
        bfm_phys_vars.vertical_extinction[:,group] = bfm_phys_vars.vertical_extinction[:,group] + phyto1_parameters["c_P"] * d3state[:,13]
    elif group == 1: # P2: Flagellates
        bfm_phys_vars.vertical_extinction[:,group] = bfm_phys_vars.vertical_extinction[:,group] + phyto2_parameters["c_P"] * d3state[:,18]
    elif group == 2: # P3: PicoPhytoplankton
        bfm_phys_vars.vertical_extinction[:,group] = bfm_phys_vars.vertical_extinction[:,group] + phyto3_parameters["c_P"] * d3state[:,22]
    elif group == 3: # P4: Large Phytoplankton
        bfm_phys_vars.vertical_extinction[:,group] = bfm_phys_vars.vertical_extinction[:,group] + phyto4_parameters["c_P"] * d3state[:,26]
    
    # # from CalcVerticalExtinction.F90 line 82
    # bfm_phys_vars.vertical_extinction[:,group] = env_parameters["p_eps0"] + env_parameters["p_epsESS"]*bfm_phys_vars.suspended_matter + env_parameters["p_epsR6"]*d3state[44,:]
    # # from CalcVerticalExtinction.F90 line 101 (ChlAttenFlag=1, ChlDynamicsFlag=2)       
    # if group == 0: # P1: Diatoms
    #     bfm_phys_vars.vertical_extinction[:,group] = bfm_phys_vars.vertical_extinction[:,group] + phyto1_parameters["c_P"] * d3state[13,:]
    # elif group == 1: # P2: Flagellates
    #     bfm_phys_vars.vertical_extinction[:,group] = bfm_phys_vars.vertical_extinction[:,group] + phyto2_parameters["c_P"] * d3state[18,:]
    # elif group == 2: # P3: PicoPhytoplankton
    #     bfm_phys_vars.vertical_extinction[:,group] = bfm_phys_vars.vertical_extinction[:,group] + phyto3_parameters["c_P"] * d3state[22,:]
    # elif group == 3: # P4: Large Phytoplankton
    #     bfm_phys_vars.vertical_extinction[:,group] = bfm_phys_vars.vertical_extinction[:,group] + phyto4_parameters["c_P"] * d3state[26,:]
        
    return bfm_phys_vars


def calculate_light_distribution(bfm_phys_vars, group):
    # Calcluations taken from bfm/bfm50/Functions/phyto.py - Updated from 0D to 1D
    # Array calculations from CalcLightDistribution.F90

    bfm_phys_vars.irradiance[0,group] = bfm_phys_vars.irradiance[0,group]*env_parameters["epsilon_PAR"]/constant_parameters["e2w"]
    num_boxes = vertical_layers - 1
    for i in range(1,num_boxes):
        bfm_phys_vars.irradiance[i,group] = bfm_phys_vars.irradiance[i-1,group] * np.exp( -1. * bfm_phys_vars.vertical_extinction[i-1,group] * bfm_phys_vars.depth[i-1])

    return bfm_phys_vars


# def calculate_vertical_extinction(bfm_phys_vars, d3state):
#     # Calcluations taken from bfm/bfm50/Functions/phyto.py - Updated from 0D to 1D
    
#     # from CalcVerticalExtinction.F90 line 82
#     bfm_phys_vars.vertical_extinction = env_parameters["p_eps0"] + env_parameters["p_epsESS"]*bfm_phys_vars.suspended_matter + env_parameters["p_epsR6"]*d3state[:,44]
#     # from CalcVerticalExtinction.F90 line 101 (ChlAttenFlag=1, ChlDynamicsFlag=2)       
#     bfm_phys_vars.vertical_extinction = bfm_phys_vars.vertical_extinction + d3state[:,13] + d3state[:,18]*2 + d3state[:,22]*3 + d3state[:,26]*4

#     return bfm_phys_vars


# def calculate_light_distribution(bfm_phys_vars):
#     # Calcluations taken from bfm/bfm50/Functions/phyto.py - Updated from 0D to 1D
#     # Array calculations from CalcLightDistribution.F90

#     bfm_phys_vars.irradiance[0] = bfm_phys_vars.irradiance[0]*env_parameters["epsilon_PAR"]/constant_parameters["e2w"]
#     num_boxes = vertical_layers - 1
#     for i in range(1,num_boxes):
#         bfm_phys_vars.irradiance[i] = bfm_phys_vars.irradiance[i-1] * np.exp( -1. * bfm_phys_vars.vertical_extinction[i-1] * bfm_phys_vars.depth[i-1])


#     return bfm_phys_vars
