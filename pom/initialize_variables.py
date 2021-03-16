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
    from pom.modules import vertical_layers
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
    wind_speed_zonal   = np.zeros(array_length,dtype=float)
    wind_speed_meridional   = np.zeros(array_length,dtype=float)
    for i in range(0,array_length):
        wind_speed_zonal[i] = wind_speed_data[2*i + 0]
        wind_speed_meridional[i] = wind_speed_data[2*i + 1]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SURFACE SALINITY
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(surface_salinity_data_path):
        path_error(surface_salinity_data_path)
    surface_salinity_data = np.fromfile(surface_salinity_data_path,dtype=float)
    surface_salinity   = np.zeros(array_length,dtype=float)
    for i in range(0,array_length):
        surface_salinity[i] = surface_salinity_data[i]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   RADIANCE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(shortwave_solar_radiation_data_path):
        path_error(shortwave_solar_radiation_data_path)
    shortwave_solar_radiation_data = np.fromfile(shortwave_solar_radiation_data_path,dtype=float)
    solar_radiation  = np.zeros(array_length,dtype=float)
    for i in range(0,array_length):
        solar_radiation[i] = shortwave_solar_radiation_data[i]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INORGANIC SUSPENDED MATTER
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(inorganic_suspended_matter_data_path):
        path_error(inorganic_suspended_matter_data_path)
    inorganic_suspended_matter_data = np.fromfile(inorganic_suspended_matter_data_path,dtype=float)
    inorganic_suspended_matter   = np.zeros((array_length, vertical_layers), dtype=float)
    for i in range(0,array_length):
        for x in range(0, vertical_layers):
            inorganic_suspended_matter[i][x] = inorganic_suspended_matter_data[vertical_layers * i + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SALINITY CLIMATOLOGY (DIAGNOSTIC MODE)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(salinity_vertical_profile_data_path):
        path_error(salinity_vertical_profile_data_path)
    salinity_vertical_profile_data = np.fromfile(salinity_vertical_profile_data_path,dtype=float)
    salinity_climatology = np.zeros((array_length, vertical_layers), dtype=float)
    for i in range(0,array_length):
        for x in range(0, vertical_layers):
            salinity_climatology[i][x] = salinity_vertical_profile_data[vertical_layers * i + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   TEMPERATURE CLIMATOLOGY (DIAGNOSTIC MODE)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(temperature_vertical_profile_data_path):
        path_error(temperature_vertical_profile_data_path)
    temperature_vertical_profile_data = np.fromfile(temperature_vertical_profile_data_path,dtype=float)
    temperature_climatology = np.zeros((array_length, vertical_layers), dtype=float)
    for i in range(0,array_length):
        for x in range(0, vertical_layers):
            temperature_climatology[i][x] = temperature_vertical_profile_data[vertical_layers * i + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   GENERAL CIRCULATION W VELOITY CLIMATOLOGY
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(general_circulation_w_velocity_data_path):
        path_error(general_circulation_w_velocity_data_path)
    general_circulation_w_velocity_data = np.fromfile(general_circulation_w_velocity_data_path,dtype=float)
    w_velocity_climatology  = np.zeros((array_length, vertical_layers), dtype=float)
    for i in range(0,array_length):
        for x in range(0, vertical_layers):
            w_velocity_climatology[i][x] = general_circulation_w_velocity_data[vertical_layers * i + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INTERMEDIATE EDDY W VELOCITY 1
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(intermediate_eddy_w_velocity_1_data_path):
        path_error(intermediate_eddy_w_velocity_1_data_path)
    intermediate_eddy_w_velocity_1_data = np.fromfile(intermediate_eddy_w_velocity_1_data_path,dtype=float)
    w_eddy_velocity_1  = np.zeros((array_length, vertical_layers), dtype=float)
    for i in range(0,array_length):
        for x in range(0, vertical_layers):
            w_eddy_velocity_1[i][x] = intermediate_eddy_w_velocity_1_data[vertical_layers * i + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INTERMEDIATE EDDY W VELOCITY 2
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(intermediate_eddy_w_velocity_2_data_path):
        path_error(intermediate_eddy_w_velocity_2_data_path)
    intermediate_eddy_w_velocity_2_data = np.fromfile(intermediate_eddy_w_velocity_2_data_path,dtype=float)
    w_eddy_velocity_2 = np.zeros((array_length, vertical_layers), dtype=float)
    for i in range(0,array_length):
        for x in range(0, vertical_layers):
            w_eddy_velocity_2[i][x] = intermediate_eddy_w_velocity_2_data[vertical_layers * i + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SALINITY INITIAL PROFILE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(salinity_initial_conditions_data_path):
        path_error(salinity_initial_conditions_data_path)
    salinity_initial_conditions_data = np.fromfile(salinity_initial_conditions_data_path,dtype=float)
    salinity_initial_profile = np.zeros(vertical_layers, dtype=float)
    for i in range(0, vertical_layers):
        salinity_initial_profile[i] = salinity_initial_conditions_data[i]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   TEMPERATURE INITIAL PROFILE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(temperature_initial_conditions_data_path):
        path_error(temperature_initial_conditions_data_path)
    temperature_initial_conditions_data = np.fromfile(temperature_initial_conditions_data_path,dtype=float)
    temperature_initial_profile = np.zeros(vertical_layers, dtype=float)
    for i in range(0, vertical_layers):
        temperature_initial_profile[i] = temperature_initial_conditions_data[i]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   HEAT FLUX
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(heat_flux_loss_data_path):
        path_error(heat_flux_loss_data_path)
    heat_flux_loss_data = np.fromfile(heat_flux_loss_data_path,dtype=float)
    surface_solar_radiation = np.zeros(array_length,dtype=float)
    surface_heat_flux_loss = np.zeros(array_length,dtype=float)
    kinetic_energy_loss = np.zeros(array_length,dtype=float)
    for i in range(0,array_length):
        surface_solar_radiation[i]  = heat_flux_loss_data[3*i + 0]
        surface_heat_flux_loss[i] = heat_flux_loss_data[3*i + 1]
        kinetic_energy_loss[i]  = heat_flux_loss_data[3*i + 2]

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
    for i in range(0,array_length):
        NO3_s1[i]  = surface_nutrients_data[4*i + 0]
        NH4_s1[i]  = surface_nutrients_data[4*i + 1]
        PO4_s1[i]  = surface_nutrients_data[4*i + 2]
        SIO4_s1[i] = surface_nutrients_data[4*i + 3]

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
    for i in range(0,array_length):
        O2_b1[i]  = bottom_nutrients_data[4*i + 0]
        NO3_b1[i] = bottom_nutrients_data[4*i + 1]
        PO4_b1[i] = bottom_nutrients_data[4*i + 2]
        PON_b1[i] = bottom_nutrients_data[4*i + 3]

    return wind_speed_zonal, wind_speed_meridional, surface_salinity, solar_radiation, inorganic_suspended_matter, \
           salinity_climatology, temperature_climatology, w_velocity_climatology, w_eddy_velocity_1, \
           w_eddy_velocity_2, salinity_initial_profile, temperature_initial_profile, \
           surface_solar_radiation, surface_heat_flux_loss, kinetic_energy_loss, \
           NO3_s1, NH4_s1, PO4_s1, SIO4_s1, O2_b1, NO3_b1, PO4_b1, PON_b1


# wind_speed_zonal, wind_speed_meridional, surface_salinity, solar_radiation, inorganic_suspended_matter, \
#            salinity_climatology, temperature_climatology, w_velocity_climatology, w_eddy_velocity_1, \
#            w_eddy_velocity_2, salinity_initial_profile, temperature_initial_profile, \
#            surface_solar_radiation, surface_heat_flux_loss, kinetic_energy_loss, \
#            NO3_s1, NH4_s1, PO4_s1, SIO4_s1, O2_b1, NO3_b1, PO4_b1, PON_b1                            = read_pom_input()




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
    from pom.modules import vertical_layers
    from pom.modules import path_error

    salinity_initial_conditions_data_path = current_path + '/inputs/POM_BFM17/init_prof_S_150m_bermuda_killworth2.da'
    temperature_initial_conditions_data_path = current_path + '/inputs/POM_BFM17/init_prof_T_150m_bermuda_killworth2.da'

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SALINITY INITIAL PROFILE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(salinity_initial_conditions_data_path):
        path_error(salinity_initial_conditions_data_path)
    salinity_initial_conditions_data = np.fromfile(salinity_initial_conditions_data_path,dtype=float)
    salinity_backward_time_level = np.zeros(vertical_layers, dtype=float)
    for i in range(0, vertical_layers):
        salinity_backward_time_level[i] = salinity_initial_conditions_data[i]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   TEMPERATURE INITIAL PROFILE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(temperature_initial_conditions_data_path):
        path_error(temperature_initial_conditions_data_path)
    temperature_initial_conditions_data = np.fromfile(temperature_initial_conditions_data_path,dtype=float)
    temperature_backward_time_level = np.zeros(vertical_layers, dtype=float)
    for i in range(0, vertical_layers):
        temperature_backward_time_level[i] = temperature_initial_conditions_data[i]

    temperature_current_time_level = np.zeros(vertical_layers, dtype=float)
    salinity_current_time_level = np.zeros(vertical_layers, dtype=float)

    temperature_current_time_level[:] = temperature_backward_time_level[:]
    salinity_current_time_level[:] = salinity_backward_time_level[:]

    return temperature_current_time_level, temperature_backward_time_level, \
           salinity_current_time_level, salinity_backward_time_level


# temperature_current_time_level, temperature_backwards_time_level, \
#     salinity_current_time_level, salinity_backwards_time_level      = get_temperature_and_salinity_initial_coditions()

#   EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
