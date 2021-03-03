# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
import numpy as np
from os import path
from cppdefs import *

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# SUBROUTINE: read_pom_input
#
# DESCRIPTION: Opens forcing files reading path specified in pom_input nml.
# (formerly opendat)
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def read_pom_input():

    from main_pombfm1d import current_path
    from pom.modules import KB
    from pom.modules import path_error

    # PATHS TO INPUT DATA FILES
    wind_stress_data_path = current_path + '/inputs/POM_BFM17/monthly_surf_wind_stress_bermuda_killworth2.da'
    surface_salinity_data_path = current_path + '/inputs/POM_BFM17/monthly_surf_salt_bermuda_150m_killworth2.da'
    shortwave_solar_radiation_data_path = current_path + '/inputs/POM_BFM17/monthly_surf_qs_bermuda_killworth2.da'
    inorganic_suspended_matter_data_path = current_path + '/inputs/POM_BFM17/monthly_clima_ISM_150m_bermuda_killworth.da'
    salinity_vertical_profile_data_path = current_path + '/inputs/POM_BFM17/monthly_clima_salt_150m_bermuda_killworth2.da'
    temperature_vertical_profile_data_path = current_path + '/inputs/POM_BFM17/monthly_clima_temp_150m_bermuda_killworth2.da'
    general_circulation_w_velocity_data_path = current_path + '/inputs/POM_BFM17/monthly_clima_w_150m_bermuda_ekman.da'
    intermediate_eddy_w_velocity_1_data_path = current_path + '/inputs/POM_BFM17/bimonthly_random_eddy_w_150m_bermuda_norm1.da'
    intermediate_eddy_w_velocity_2_data_path = current_path + '/inputs/POM_BFM17/bimonthly_random_eddy_w_150m_bermuda_norm2.da'
    salinity_initial_conditions_data_path = current_path + '/inputs/POM_BFM17/init_prof_S_150m_bermuda_killworth2.da'
    temperature_initial_conditions_data_path = current_path + '/inputs/POM_BFM17/init_prof_T_150m_bermuda_killworth2.da'
    heat_flux_loss_data_path = current_path + '/inputs/POM_BFM17/monthly_surf_rad_bermuda_killworth2.da'
    surface_nutrients_data_path = current_path + '/inputs/POM_BFM17/NutrientsARPAOGS.da'
    bottom_nutrients_data_path = current_path + '/inputs/POM_BFM17/monthly_bott_nut_bermuda_150m_killworth.da'

    # LENGTH OF INPUT ARRAYS
    array_length = 13   # MONTHS (D-J-F-M-A-M-J-J-A-S-O-N-D)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   WIND SPEED (u,v)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(wind_stress_data_path):
        path_error(wind_stress_data_path)
    wind_speed_data = np.fromfile(wind_stress_data_path,dtype=float)
    wind_speed_u   = np.zeros(array_length,dtype=float)
    wind_speed_v   = np.zeros(array_length,dtype=float)
    for k in range(0,array_length):
        wind_speed_u[k] = wind_speed_data[2*k + 0]
        wind_speed_v[k] = wind_speed_data[2*k + 1]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SURFACE SALINITY
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(surface_salinity_data_path):
        path_error(surface_salinity_data_path)
    surface_salinity_data = np.fromfile(surface_salinity_data_path,dtype=float)
    surface_salinity   = np.zeros(array_length,dtype=float)
    for k in range(0,array_length):
        surface_salinity[k] = surface_salinity_data[k]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   RADIANCE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(shortwave_solar_radiation_data_path):
        path_error(shortwave_solar_radiation_data_path)
    shortwave_solar_radiation_data = np.fromfile(shortwave_solar_radiation_data_path,dtype=float)
    solar_radiation  = np.zeros(array_length,dtype=float)
    for k in range(0,array_length):
        solar_radiation[k] = shortwave_solar_radiation_data[k]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INORGANIC SUSPENDED MATTER
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(inorganic_suspended_matter_data_path):
        path_error(inorganic_suspended_matter_data_path)
    inorganic_suspended_matter_data = np.fromfile(inorganic_suspended_matter_data_path,dtype=float)
    inorganic_suspended_matter   = np.zeros((array_length,KB),dtype=float)
    for k in range(0,array_length):
        for x in range(0,KB):
            inorganic_suspended_matter[k][x] = inorganic_suspended_matter_data[KB*k + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SALINITY CLIMATOLOGY (DIAGNOSTIC MODE)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(salinity_vertical_profile_data_path):
        path_error(salinity_vertical_profile_data_path)
    salinity_vertical_profile_data = np.fromfile(salinity_vertical_profile_data_path,dtype=float)
    salinity_climatology = np.zeros((array_length,KB),dtype=float)
    for k in range(0,array_length):
        for x in range(0,KB):
            salinity_climatology[k][x] = salinity_vertical_profile_data[KB*k + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   TEMPERATURE CLIMATOLOGY (DIAGNOSTIC MODE)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(temperature_vertical_profile_data_path):
        path_error(temperature_vertical_profile_data_path)
    temperature_vertical_profile_data = np.fromfile(temperature_vertical_profile_data_path,dtype=float)
    temperature_climatology = np.zeros((array_length,KB),dtype=float)
    for k in range(0,array_length):
        for x in range(0,KB):
            temperature_climatology[k][x] = temperature_vertical_profile_data[KB*k + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   GENERAL CIRCULATION W VELOITY CLIMATOLOGY
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(general_circulation_w_velocity_data_path):
        path_error(general_circulation_w_velocity_data_path)
    general_circulation_w_velocity_data = np.fromfile(general_circulation_w_velocity_data_path,dtype=float)
    w_velocity_climatology  = np.zeros((array_length,KB),dtype=float)
    for k in range(0,array_length):
        for x in range(0,KB):
            w_velocity_climatology[k][x] = general_circulation_w_velocity_data[KB*k + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INTERMEDIATE EDDY W VELOCITY 1
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(intermediate_eddy_w_velocity_1_data_path):
        path_error(intermediate_eddy_w_velocity_1_data_path)
    intermediate_eddy_w_velocity_1_data = np.fromfile(intermediate_eddy_w_velocity_1_data_path,dtype=float)
    w_eddy_velocity_1  = np.zeros((array_length,KB),dtype=float)
    for k in range(0,array_length):
        for x in range(0,KB):
            w_eddy_velocity_1[k][x] = intermediate_eddy_w_velocity_1_data[KB*k + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INTERMEDIATE EDDY W VELOCITY 2
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(intermediate_eddy_w_velocity_2_data_path):
        path_error(intermediate_eddy_w_velocity_2_data_path)
    intermediate_eddy_w_velocity_2_data = np.fromfile(intermediate_eddy_w_velocity_2_data_path,dtype=float)
    w_eddy_velocity_2 = np.zeros((array_length,KB),dtype=float)
    for k in range(0,array_length):
        for x in range(0,KB):
            w_eddy_velocity_2[k][x] = intermediate_eddy_w_velocity_2_data[KB*k + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SALINITY INITIAL PROFILE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(salinity_initial_conditions_data_path):
        path_error(salinity_initial_conditions_data_path)
    salinity_initial_conditions_data = np.fromfile(salinity_initial_conditions_data_path,dtype=float)
    salinity_initial_profile = np.zeros(KB,dtype=float)
    for k in range(0,KB):
        salinity_initial_profile[k] = salinity_initial_conditions_data[k]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   TEMPERATURE INITIAL PROFILE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(temperature_initial_conditions_data_path):
        path_error(temperature_initial_conditions_data_path)
    temperature_initial_conditions_data = np.fromfile(temperature_initial_conditions_data_path,dtype=float)
    temperature_initial_profile = np.zeros(KB,dtype=float)
    for k in range(0,KB):
        temperature_initial_profile[k] = temperature_initial_conditions_data[k]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   HEAT FLUX
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(heat_flux_loss_data_path):
        path_error(heat_flux_loss_data_path)
    heat_flux_loss_data = np.fromfile(heat_flux_loss_data_path,dtype=float)
    surface_solar_radiation = np.zeros(array_length,dtype=float)
    surface_heat_flux_loss = np.zeros(array_length,dtype=float)
    coriolis_heat_flux_loss = np.zeros(array_length,dtype=float)
    for k in range(0,array_length):
        surface_solar_radiation[k]  = heat_flux_loss_data[3*k + 0]
        surface_heat_flux_loss[k] = heat_flux_loss_data[3*k + 1]
        coriolis_heat_flux_loss[k]  = heat_flux_loss_data[3*k + 2]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SURFACE NUTRIENTS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(surface_nutrients_data_path):
        path_error(surface_nutrients_data_path)
    surface_nutrients_data  = np.fromfile(surface_nutrients_data_path,dtype=float)
    NO3_s1  = np.zeros(array_length,dtype=float)
    NH4_s1  = np.zeros(array_length,dtype=float)
    PO4_s1  = np.zeros(array_length,dtype=float)
    SIO4_s1 = np.zeros(array_length,dtype=float)
    for k in range(0,array_length):
        NO3_s1[k]  = surface_nutrients_data[4*k + 0]
        NH4_s1[k]  = surface_nutrients_data[4*k + 1]
        PO4_s1[k]  = surface_nutrients_data[4*k + 2]
        SIO4_s1[k] = surface_nutrients_data[4*k + 3]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   BOTTOM NUTRIENTS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(bottom_nutrients_data_path):
        path_error(bottom_nutrients_data_path)
    bottom_nutrients_data = np.fromfile(bottom_nutrients_data_path,dtype=float)
    O2_b1   = np.zeros(array_length,dtype=float)
    NO3_b1  = np.zeros(array_length,dtype=float)
    PO4_b1  = np.zeros(array_length,dtype=float)
    PON_b1  = np.zeros(array_length,dtype=float)
    for k in range(0,array_length):
        O2_b1[k]  = bottom_nutrients_data[4*k + 0]
        NO3_b1[k] = bottom_nutrients_data[4*k + 1]
        PO4_b1[k] = bottom_nutrients_data[4*k + 2]
        PON_b1[k] = bottom_nutrients_data[4*k + 3]

    return wind_speed_u, wind_speed_v, surface_salinity, solar_radiation, inorganic_suspended_matter, \
           salinity_climatology, temperature_climatology, w_velocity_climatology, w_eddy_velocity_1, \
           w_eddy_velocity_2, salinity_initial_profile, temperature_initial_profile, \
           surface_solar_radiation, surface_heat_flux_loss, coriolis_heat_flux_loss, \
           NO3_s1, NH4_s1, PO4_s1, SIO4_s1, O2_b1, NO3_b1, PO4_b1, PON_b1


# wind_speed_u, wind_speed_v, surface_salinity, solar_radiation, inorganic_suspended_matter, \
#     salinity_climatology, temperature_climatology, w_velocity_climatology, w_eddy_velocity_1, \
#     w_eddy_velocity_2, salinity_initial_profile, temperature_initial_profile, \
#     surface_solar_radiation, surface_heat_flux_loss, coriolis_heat_flux_loss, \
#     NO3_s1, NH4_s1, PO4_s1, SIO4_s1, O2_b1, NO3_b1, PO4_b1, PON_b1                              = read_pom_input()




# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# SUBROUTINE: get_TS_IC
#
# DESCRIPTION: This subroutine opens and reads files containing the T&S initial conditions
#              Files are read in direct access mode reading path specified in pom_input nml
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


def get_temperature_and_salinity_initial_coditions():

    from main_pombfm1d import current_path
    from pom.modules import KB
    from pom.modules import path_error

    salinity_initial_conditions_data_path = current_path + '/inputs/POM_BFM17/init_prof_S_150m_bermuda_killworth2.da'
    temperature_initial_conditions_data_path = current_path + '/inputs/POM_BFM17/init_prof_T_150m_bermuda_killworth2.da'

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SALINITY INITIAL PROFILE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(salinity_initial_conditions_data_path):
        path_error(salinity_initial_conditions_data_path)
    salinity_initial_conditions_data = np.fromfile(salinity_initial_conditions_data_path,dtype=float)
    salinity_backwards_time_level = np.zeros(KB,dtype=float)
    for k in range(0, KB):
        salinity_backwards_time_level[k] = salinity_initial_conditions_data[k]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   TEMPERATURE INITIAL PROFILE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(temperature_initial_conditions_data_path):
        path_error(temperature_initial_conditions_data_path)
    temperature_initial_conditions_data = np.fromfile(temperature_initial_conditions_data_path,dtype=float)
    temperature_backwards_time_level = np.zeros(KB,dtype=float)
    for k in range(0, KB):
        temperature_backwards_time_level[k] = temperature_initial_conditions_data[k]

    temperature_current_time_level = np.zeros(KB,dtype=float)
    salinity_current_time_level = np.zeros(KB,dtype=float)

    temperature_current_time_level[:] = temperature_backwards_time_level[:]
    salinity_current_time_level[:] = salinity_backwards_time_level[:]

    return temperature_current_time_level, temperature_backwards_time_level, \
           salinity_current_time_level, salinity_backwards_time_level


# temperature_current_time_level, temperature_backwards_time_level, \
#     salinity_current_time_level, salinity_backwards_time_level      = get_temperature_and_salinity_initial_coditions()

#   EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
