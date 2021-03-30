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

    """
    Description: Opens forcing files reading the paths specified in the pom_input namelist.

    :return: data arrays for wind stress, surface salinity, solar radiation, inorganic
             suspended matter, salinity and temperature vertical profiles, general circulation
             for w velocity, intermediate eddy velocities, salinity and temperature initial
             conditions, heat flux loss, and surface and bottom nutrients
    """

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
    shortwave_radiation = np.zeros(array_length,dtype=float)
    surface_heat_flux_loss = np.zeros(array_length,dtype=float)
    kinetic_energy_loss = np.zeros(array_length,dtype=float)
    for i in range(0,array_length):
        shortwave_radiation[i]  = heat_flux_loss_data[3*i + 0]
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
           shortwave_radiation, surface_heat_flux_loss, kinetic_energy_loss, \
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

    """
    Description: Opens and reads files containing the initial conditions for temperature
                 and salinity. Files are read from the pom_input namelist.

    :return: temperature and salinity at the current and backward time level
    """

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


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# ROUTINE: GetDelta
#
# DESCRIPTION: Get the numeric timestep
#              Transfer the integration time step to the BFM Unit conversion
#              from seconds to days
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def get_numeric_timestep():

    """
    Description: Get the numeric timestep. Transfer the integration timestep to the BFM
                 Unit conversion from seconds to days.

    :return: numeric timestep
    """

    from modules import seconds_per_day, timestep

    numeric_timestep = timestep/seconds_per_day

    return numeric_timestep


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# ROUTINE: Service
#
# DESCRIPTION: This routine passes the physical variables to the BFM
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def pom_to_bfm():

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
    #     from modules import seawater_temperature, seawater_salinity, photosynthetic_radiation, suspended_sediment_load, seawater_density, wind_speed, gridpoint_depth

    from modules import vertical_layers, water_specific_heat_times_density, temperature_backward_time_level, \
        salinity_backward_time_level, density, bottom_depth, vertical_spacing, shortwave_radiation, wind_stress_zonal, \
        wind_stress_meridional, diffusion_coefficient_momentum, diffusion_coefficient_tracers, \
        velocity_zonal_current_time_level, velocity_meridional_current_time_level, kinetic_energy_current_time_level, \
        kinetic_energy_times_length_current_time_level, length_scale

    from modules import inorganic_suspended_matter, w_velocity_interpolation, w_eddy_velocity

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   1D ARRAYS FOR BFM
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    seawater_temperature = np.zeros(vertical_layers,dtype=float)
    seawater_salinity = np.zeros(vertical_layers,dtype=float)
    suspended_sediment_load = np.zeros(vertical_layers,dtype=float)
    seawater_density = np.zeros(vertical_layers,dtype=float)
    gridpoint_depth = np.zeros(vertical_layers,dtype=float)
    photosynthetic_radiation = np.zeros(vertical_layers,dtype=float)

    for i in range(0,vertical_layers):
        seawater_temperature[i] = temperature_backward_time_level[i]
        seawater_salinity[i] = salinity_backward_time_level[i]
        seawater_density[i] = (density[i] * 1.E3) + 1.E3
        suspended_sediment_load[i] = inorganic_suspended_matter[i]
        gridpoint_depth[i] = vertical_spacing[i] * bottom_depth

    photosynthetic_radiation[0] = -1. * shortwave_radiation * water_specific_heat_times_density

    wind_stress = np.sqrt(wind_stress_zonal**2 + wind_stress_meridional**2) * 1.E3
    wind_speed = np.sqrt(wind_stress/(1.25 * 0.0014))

    return seawater_temperature, seawater_salinity, photosynthetic_radiation, \
           suspended_sediment_load, seawater_density, wind_speed


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: get_IC
#
# DESCRIPTION
#
# This subroutine opens and read files containing the initial conditions
# for P2c, Z4c, R1c, R6c, N1p, N3n, N4n, and O2o
# Files are read in direct access mode reading path
# specified in bfm_IC_input nml
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def get_initial_conditions():

    """
    Desription: Opens and reads data files containing the initial conditions for P2c, Z4c,
                R1c, R6c, N1p, N3n, N4n, and O2o

    :return: initial conditions for phytoplankton carbon, zooplankton carbon, particulate
             organic carbon, dissolved organic carbon, phosphate, nitrate, ammonium,
             and oxygen
    """

    from main_pombfm1d import current_path
    from modules import vertical_layers, path_error
    from inputs.namelist_input_data import phyto_input, zoop_input, poc_input, doc_input, \
        phos_input, nit_input, am_input, oxy_input

    phytoplankton_carbon_data_path = current_path + phyto_input
    zooplankton_carbon_data_path = current_path + zoop_input
    particulate_organic_carbon_data_path = current_path + poc_input
    dissolved_organic_carbon_data_path = current_path + doc_input
    phosphate_data_path = current_path + phos_input
    nitrate_data_path = current_path + nit_input
    ammonium_data_path = current_path + am_input
    oxygen_data_path = current_path + oxy_input

    if not path.exists(phytoplankton_carbon_data_path):
        path_error(phytoplankton_carbon_data_path)
    p2c = np.zeros(vertical_layers,dtype=float)
    p2c = np.fromfile(phytoplankton_carbon_data_path,dtype=float)

    if not path.exists(zooplankton_carbon_data_path):
        path_error(zooplankton_carbon_data_path)
    z5c = np.zeros(vertical_layers, dtype=float)
    z5c = np.fromfile(zooplankton_carbon_data_path,dtype=float)

    if not path.exists(particulate_organic_carbon_data_path):
        path_error(particulate_organic_carbon_data_path)
    r6c = np.zeros(vertical_layers,dtype=float)
    r6c = np.fromfile(particulate_organic_carbon_data_path,dtype=float)

    if not path.exists(dissolved_organic_carbon_data_path):
        path_error(dissolved_organic_carbon_data_path)
    r1c = np.zeros(vertical_layers,dtype=float)
    r1c = np.fromfile(dissolved_organic_carbon_data_path,dtype=float)

    if not path.exists(phosphate_data_path):
        path_error(phosphate_data_path)
    n1p = np.zeros(vertical_layers,dtype=float)
    n1p = np.fromfile(phosphate_data_path,dtype=float)

    if not path.exists(nitrate_data_path):
        path_error(nitrate_data_path)
    n3n = np.zeros(vertical_layers,dtype=float)
    n3n = np.fromfile(nitrate_data_path,dtype=float)

    if not path.exists(ammonium_data_path):
        path_error(ammonium_data_path)
    n4n = np.zeros(vertical_layers,dtype=float)
    n4n = np.fromfile(ammonium_data_path,dtype=float)

    if not path.exists(oxygen_data_path):
        path_error(oxygen_data_path)
    o2o = np.zeros(vertical_layers,dtype=float)
    o2o = np.fromfile(oxygen_data_path,dtype=float)

    return p2c, z5c, r6c, r1c, n1p, n3n, n4n, o2o

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: set_initial_conditions
#
# DESCRIPTION
#
# This routine assigns initial conditioons of biochemical variables in POM
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def set_initial_conditions():

    """
    Description: Assigns initial conditions of biochemical variables in POM

    :return:
    """

    from modules import vertical_layers
    from modules import zero, NML_OPEN, NML_READ, NMLUNIT, error_msg_prn
    # from modules import photosynthetic_radiation

    p_nRc = 0.0126
    p_pRc = 0.7862E-03
    p_sRc = 0.0118
    p_iRc = 1./25.

    try:
        import INCLUDE_BEN
        INCLUDE_BEN = True
    except FileNotFoundError:
        INCLUDE_BEN = False
    if INCLUDE_BEN:
        from inputs.namelist_input_data import y1c0, y2c0, y3c0, y4c0, y5c0, h1c0, h2c0, \
            k1p0, k11p0, k21p0, k4n0, k14n0, k24n0, k3n0, k5s0, k6r0, d1m0, d2m0, d6m0, d7m0, d8m0, d9m0, \
            q6c0, q6n0, q6p0, q6s0, q1c0, q11c0, g2o0, g3c0, g13c0, g23c0, g3h0, g13h0, g23h0

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   DEFINITION OF BIOGEOCHEMICAL GLOBAL VARIABLES
    #   IrrOPT  in the equation of Steele and light
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    photosynthetic_radiation = np.zeros(vertical_layers,dtype=float)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   DEFINITION OF GENERAL PELAGIC STATE VARIABLES:
    #   PELAGIC GASES
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    p2c, z5c, r6c, r1c, n1p, n3n, n4n, o2o = get_initial_conditions()

    # NEED TO FINISH TRANSLATING MODULEMEM




#   EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
