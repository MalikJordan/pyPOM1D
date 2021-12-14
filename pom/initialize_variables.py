# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
import numpy as np
from os import path
from cppdefs import *
from inputs import params_POMBFM
from inputs.namelist_input_data import phyto_input, zoop_input, poc_input, doc_input, phos_input, nit_input, am_input, oxy_input
from pom.constants import current_path, seconds_per_day, vertical_layers

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
    if path.exists(wind_stress_data_path):
        wind_speed_data = np.fromfile(wind_stress_data_path)
    wind_speed_zonal   = np.zeros(array_length)
    wind_speed_meridional   = np.zeros(array_length)
    for i in range(0,array_length):
        wind_speed_zonal[i] = wind_speed_data[2*i + 0]
        wind_speed_meridional[i] = wind_speed_data[2*i + 1]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SURFACE SALINITY
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(surface_salinity_data_path):
        surface_salinity = np.fromfile(surface_salinity_data_path)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   RADIANCE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(shortwave_solar_radiation_data_path):
        solar_radiation = np.fromfile(shortwave_solar_radiation_data_path)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INORGANIC SUSPENDED MATTER
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(inorganic_suspended_matter_data_path):
        inorganic_suspended_matter_data = np.fromfile(inorganic_suspended_matter_data_path)
    inorganic_suspended_matter   = np.zeros((vertical_layers,array_length))
    for i in range(0,array_length):
        for x in range(0, vertical_layers):
            inorganic_suspended_matter[x,i] = inorganic_suspended_matter_data[vertical_layers * i + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SALINITY CLIMATOLOGY (DIAGNOSTIC MODE)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(salinity_vertical_profile_data_path):
        salinity_vertical_profile_data = np.fromfile(salinity_vertical_profile_data_path)
    salinity_climatology = np.zeros((vertical_layers,array_length))
    for i in range(0,array_length):
        for x in range(0, vertical_layers):
            salinity_climatology[x,i] = salinity_vertical_profile_data[vertical_layers * i + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   TEMPERATURE CLIMATOLOGY (DIAGNOSTIC MODE)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(temperature_vertical_profile_data_path):
        temperature_vertical_profile_data = np.fromfile(temperature_vertical_profile_data_path)
    temperature_climatology = np.zeros((vertical_layers,array_length))
    for i in range(0,array_length):
        for x in range(0, vertical_layers):
            temperature_climatology[x,i] = temperature_vertical_profile_data[vertical_layers * i + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   GENERAL CIRCULATION W VELOITY CLIMATOLOGY
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(general_circulation_w_velocity_data_path):
        general_circulation_w_velocity_data = np.fromfile(general_circulation_w_velocity_data_path)
    w_velocity_climatology  = np.zeros((vertical_layers,array_length))
    for i in range(0,array_length):
        for x in range(0, vertical_layers):
            w_velocity_climatology[x,i] = general_circulation_w_velocity_data[vertical_layers * i + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INTERMITTANT EDDY W VELOCITY 1
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(intermediate_eddy_w_velocity_1_data_path):
        intermediate_eddy_w_velocity_1_data = np.fromfile(intermediate_eddy_w_velocity_1_data_path)
    w_eddy_velocity_1  = np.zeros((vertical_layers,array_length))
    for i in range(0,array_length):
        for x in range(0, vertical_layers):
            w_eddy_velocity_1[x,i] = intermediate_eddy_w_velocity_1_data[vertical_layers * i + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INTERMITTANT EDDY W VELOCITY 2
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(intermediate_eddy_w_velocity_2_data_path):
        intermediate_eddy_w_velocity_2_data = np.fromfile(intermediate_eddy_w_velocity_2_data_path)
    w_eddy_velocity_2 = np.zeros((vertical_layers,array_length))
    for i in range(0,array_length):
        for x in range(0, vertical_layers):
            w_eddy_velocity_2[x,i] = intermediate_eddy_w_velocity_2_data[vertical_layers * i + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SALINITY INITIAL PROFILE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(salinity_initial_conditions_data_path):
        salinity_backward = np.fromfile(salinity_initial_conditions_data_path)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   TEMPERATURE INITIAL PROFILE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(temperature_initial_conditions_data_path):
        temperature_backward = np.fromfile(temperature_initial_conditions_data_path)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   HEAT FLUX
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(heat_flux_loss_data_path):
        heat_flux_loss_data = np.fromfile(heat_flux_loss_data_path)
    shortwave_radiation = np.zeros(array_length)
    surface_heat_flux = np.zeros(array_length)
    kinetic_energy_loss = np.zeros(array_length)
    for i in range(0,array_length):
        shortwave_radiation[i]  = heat_flux_loss_data[3*i + 0]
        surface_heat_flux[i] = heat_flux_loss_data[3*i + 1]
        kinetic_energy_loss[i]  = heat_flux_loss_data[3*i + 2]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SURFACE NUTRIENTS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(surface_nutrients_data_path):
        surface_nutrients_data  = np.fromfile(surface_nutrients_data_path)
    NO3_s1  = np.zeros(array_length)
    NH4_s1  = np.zeros(array_length)
    PO4_s1  = np.zeros(array_length)
    SIO4_s1 = np.zeros(array_length)
    for i in range(0,array_length):
        NO3_s1[i]  = surface_nutrients_data[4*i + 0]
        NH4_s1[i]  = surface_nutrients_data[4*i + 1]
        PO4_s1[i]  = surface_nutrients_data[4*i + 2]
        SIO4_s1[i] = surface_nutrients_data[4*i + 3]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   BOTTOM NUTRIENTS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(bottom_nutrients_data_path):
        bottom_nutrients_data = np.fromfile(bottom_nutrients_data_path)
    O2_b1   = np.zeros(array_length)
    NO3_b1  = np.zeros(array_length)
    PO4_b1  = np.zeros(array_length)
    PON_b1  = np.zeros(array_length)
    for i in range(0,array_length):
        O2_b1[i]  = bottom_nutrients_data[4*i + 0]
        NO3_b1[i] = bottom_nutrients_data[4*i + 1]
        PO4_b1[i] = bottom_nutrients_data[4*i + 2]
        PON_b1[i] = bottom_nutrients_data[4*i + 3]

    return wind_speed_zonal, wind_speed_meridional, surface_salinity, solar_radiation, inorganic_suspended_matter, \
           salinity_climatology, temperature_climatology, w_velocity_climatology, w_eddy_velocity_1, \
           w_eddy_velocity_2, salinity_backward, temperature_backward, \
           shortwave_radiation, surface_heat_flux, kinetic_energy_loss, \
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

    salinity_initial_conditions_data_path = current_path + '/inputs/POM_BFM17/init_prof_S_150m_bermuda_killworth2.da'
    temperature_initial_conditions_data_path = current_path + '/inputs/POM_BFM17/init_prof_T_150m_bermuda_killworth2.da'

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SALINITY INITIAL PROFILE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(salinity_initial_conditions_data_path):
        salinity_backward = np.fromfile(salinity_initial_conditions_data_path)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   TEMPERATURE INITIAL PROFILE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(temperature_initial_conditions_data_path):
        temperature_backward = np.fromfile(temperature_initial_conditions_data_path)

    temperature = np.zeros(vertical_layers)
    salinity = np.zeros(vertical_layers)

    temperature[:] = temperature_backward[:]
    salinity[:] = salinity_backward[:]

    return temperature, temperature_backward, \
           salinity, salinity_backward


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

    numeric_timestep = params_POMBFM.dti / seconds_per_day

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

# def pom_to_bfm():

#     """
#     Description: Passes the physical variables to the BFM

#     :return: seawater density, temperature and salinity, suspended sediment load,
#              photosynthetically available radiation, gridpoint depth, and wind speed
#     """

#     # try:
#     #     import NOPOINTERS
#     #     NOPOINTERS = True
#     # except FileNotFoundError:
#     #     NOPOINTERS = False
#     # if NOPOINTERS:
#     #     from modules import ETW, ESW, EIR, ESS, ERHO, EWIND, Depth

#     from modules_old import vertical_layers, water_specific_heat_times_density, temperature_backward, \
#         salinity_backward, density_profile, bottom_depth, vertical_spacing, shortwave_radiation, wind_stress_zonal, \
#         wind_stress_meridional, diffusion_coefficient_momentum, diffusion_coefficient_tracers, \
#         velocity_zonal, velocity_meridional, kinetic_energy, \
#         kinetic_energy_times_length, length_scale

#     from modules_old import inorganic_suspended_matter, interpolated_w_velocity, w_eddy_velocity

#     # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#     #   1D ARRAYS FOR BFM
#     # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#     ETW = np.zeros(vertical_layers - 1)
#     ESW = np.zeros(vertical_layers - 1)
#     ESS = np.zeros(vertical_layers - 1)
#     ERHO = np.zeros(vertical_layers - 1)
#     Depth = np.zeros(vertical_layers - 1)
#     EIR = np.zeros(vertical_layers - 1)

#     for i in range(0,vertical_layers - 1):
#         ETW[i] = temperature_backward[i]
#         ESW[i] = salinity_backward[i]
#         ERHO[i] = (density_profile[i] * 1.E3) + 1.E3
#         ESS[i] = inorganic_suspended_matter[i]
#         Depth[i] = vertical_spacing[i] * params_POMBFM.h

#     EIR[0] = -1. * shortwave_radiation * water_specific_heat_times_density

#     wind_stress = np.sqrt(wind_stress_zonal**2 + wind_stress_meridional**2) * 1.E3
#     EWIND = np.sqrt(wind_stress/(1.25 * 0.0014))

#     return ETW, ESW, EIR, \
#            ESS, ERHO, EWIND


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: get_IC
#
# DESCRIPTION
#
# This subroutine opens and read files containing the initial conditions
# for nanoflagellates_LO_C, omnivMesozoo_LO_C, labileDOM_NO_C, particOrganDetritus_NO_C, phospate_IO_P, nitrate_IO_N, ammonium_IO_N, and disOxygen_IO_O
# Files are read in direct access mode reading path
# specified in bfm_IC_input nml
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def get_initial_conditions():

    """
    Desription: Opens and reads data files containing the initial conditions for nanoflagellates_LO_C, omnivMesozoo_LO_C,
                labileDOM_NO_C, particOrganDetritus_NO_C, phospate_IO_P, nitrate_IO_N, ammonium_IO_N, and disOxygen_IO_O

    :return: initial conditions for phytoplankton carbon, zooplankton carbon, particulate
             organic carbon, dissolved organic carbon, phosphate, nitrate, ammonium,
             and oxygen
    """

    phytoplankton_carbon_data_path = current_path + phyto_input
    zooplankton_carbon_data_path = current_path + zoop_input
    particulate_organic_carbon_data_path = current_path + poc_input
    dissolved_organic_carbon_data_path = current_path + doc_input
    phosphate_data_path = current_path + phos_input
    nitrate_data_path = current_path + nit_input
    ammonium_data_path = current_path + am_input
    oxygen_data_path = current_path + oxy_input

    if path.exists(phytoplankton_carbon_data_path):
        nanoflagellates_LO_C = np.fromfile(phytoplankton_carbon_data_path)

    if path.exists(zooplankton_carbon_data_path):
        microzoo_LO_C = np.fromfile(zooplankton_carbon_data_path)

    if path.exists(particulate_organic_carbon_data_path):
        particOrganDetritus_NO_C = np.fromfile(particulate_organic_carbon_data_path)

    if path.exists(dissolved_organic_carbon_data_path):
        labileDOM_NO_C = np.fromfile(dissolved_organic_carbon_data_path)

    if path.exists(phosphate_data_path):
        phospate_IO_P = np.fromfile(phosphate_data_path)

    if path.exists(nitrate_data_path):
        nitrate_IO_N = np.fromfile(nitrate_data_path)

    if path.exists(ammonium_data_path):
        ammonium_IO_N = np.fromfile(ammonium_data_path)

    if path.exists(oxygen_data_path):
        disOxygen_IO_O = np.fromfile(oxygen_data_path)

    return nanoflagellates_LO_C, microzoo_LO_C, particOrganDetritus_NO_C, labileDOM_NO_C, phospate_IO_P, nitrate_IO_N, ammonium_IO_N, disOxygen_IO_O

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

# def set_initial_conditions():
#
#     """
#     Description: Assigns initial conditions of biochemical variables in POM
#
#     :return:
#     """
#
#     from modules_old import vertical_layers
#     from modules_old import zero, NML_OPEN, NML_READ, NMLUNIT, error_msg_prn
#     # from modules import photosynthetic_radiation
#
#     p_nRc = 0.0126
#     p_pRc = 0.7862E-03
#     p_sRc = 0.0118
#     p_iRc = 1./25.
#
#     try:
#         import INCLUDE_BEN
#         INCLUDE_BEN = True
#     except FileNotFoundError:
#         INCLUDE_BEN = False
#     if INCLUDE_BEN:
#         from inputs.namelist_input_data import y1c0, y2c0, y3c0, y4c0, y5c0, h1c0, h2c0, \
#             k1p0, k11p0, k21p0, k4n0, k14n0, k24n0, k3n0, k5s0, k6r0, d1m0, d2m0, d6m0, d7m0, d8m0, d9m0, \
#             q6c0, q6n0, q6p0, q6s0, q1c0, q11c0, g2o0, g3c0, g13c0, g23c0, g3h0, g13h0, g23h0
#
#     # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#     #   DEFINITION OF BIOGEOCHEMICAL GLOBAL VARIABLES
#     #   IrrOPT  in the equation of Steele and light
#     # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#     photosynthetic_radiation = np.zeros(vertical_layers,dtype=float)
#
#     # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#     #   DEFINITION OF GENERAL PELAGIC STATE VARIABLES:
#     #   PELAGIC GASES
#     # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#     nanoflagellates_LO_C, microzoo_LO_C, particOrganDetritus_NO_C, labileDOM_NO_C, phospate_IO_P, nitrate_IO_N, ammonium_IO_N, disOxygen_IO_O = get_initial_conditions()
#
#     # NEED TO FINISH TRANSLATING MODULEMEM
#
#
#
#
# #   EOC
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# #   MODEL  POM - Princeton Ocean Model
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
