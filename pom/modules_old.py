# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
import numpy as np
from cppdefs import *
from pom.initialize_variables import read_pom_input


array_length = 13
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# MODULE: GLOBALMEM
#
# DESCRIPTION: Definiton of the runtime error messages.
#              This module contains global settings:
#              - general constants for controlling prescision,
#              - parameters defining fle streams and error message numbers
#              - the subroutine for printing the error messages (functions at end of file)
#              - and aborting the simulation
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# DOUBLE PRECISION
RLEN = float()
NMLUNIT = 310

# THE UNIT OF THE LOG FILE IS NOT A PARAMETER TO ALLOW PARALLEL WRITING
LOGUNIT = 0
bfm_lwp = True
zero = 0.
one = 1.
pi = 3.14159265359
base_temperature = 20.

#   NEXT PARAMETERS ARE DEFINED TO CONTROL TYPES OF STATE VARIABLES
OFF = -100
SINKSOURCE = -1
NOTRANSPORT = 0
NOOBCSTATES = 1
HORTRANSPORT = 10
ALLTRANSPORT = 20
DONE = 1.0

#   ERROR CODES:
ALLOC = 10
NML_OPEN = 11
NML_READ = 12
DIM_MISMATCH = 13





# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# MODULE: CONSTANTS
#
# DESCRIPTION: Full list of Fortran parameters  (comparable with Sesame constants)
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# THIS INITIALIZATION HAS TO BE DONE IN MODULEPARAM BECAUSE SOME COMPILERS
# DO NOT ALLOW THE INITIALIZATION OF CONSTANTS WITH INTRINSIC FUNCTIONS

MIN_VAL_EXPFUN = zero
ECOLOGY = 1  # BASE TEMPERATURE FOR Q10
TRANSPORT = 2
zero_kelvin = -273.15
molar_gas_constant = 83.131  # GAS CONSTANT: bar mol^-1 deg-1
carbon_molecular_weight = 12.0  # MOLECULAR WEIGHT CARBON
nitrogen_molecular_weight = 14.0  # MOLECULAR WEIGHT NITROGEN
phosphorus_molecular_weight = 31.0  # MOLECULAR WEIGHT PHOSPHORUS
silica_molecular_weight = 28.0  # MOLECULAR WEIGHT SILICA
molecular_weight_conversion_factor = 0.217  # MOLECULAR WEIGHT CONVERSION FACTOR EINSTEIN->W
seconds_per_day = 86400.  # SECONDS IN DAY
days_per_second = zero  # INVERSE OF SECONDS IN DAY
one_per_day = 1.  # RATE WHICH IS USED IN CASES WHERE IMPLICITLY ASSUMED

NO_BENTHOS = 0
BENTHIC_RETURN = 1
BENTHIC_BIO = 2
BENTHIC_FULL = 3
HOURS_PER_DAY = 24.  # HOURS IN DAY
SOLAR_RADIATION = 1368.  # SOLAR RADIATION CONSTANT
NODATA = 0
DAILYDATA = 1
CLOUDDATA = 3
ALLDATA = 4
ONLYDAILY = 5

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   INTEGER CONSTANTS USED IN THE BENTHIC NUTRIENT DYNAMICS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

DEFINE = -1000
DOUBLE_DEFINE = -1100
PARAMETER_DEFINE = -1200
EQUATION = 0
SDERIVATIVE = -2
DERIVATIVE = -1
INTEGRAL = 1
EXPONENTIAL_INTEGRAL = 11
SHIFT = 21
LINEAR_TERM = 1
CONSTANT_TERM = 0
QUADRATIC_TERM = 2
EXPONENTIAL_TERM = -1
ZERO_EXPONENTIAL_TERM = -2
BESSELI_EXP_TERM = -5
BESSELK_EXP_TERM = -6
ADD = 1000
INPUT_TERM = 6001
START_ADD_TERM = 6002
INPUT_ADD_TERM = 6000
RFLUX = 1
MASS = 2
AVERAGE = 3
PARAMETER = 4
STANDARD = 0
SET_CONTINUITY = 100
SET_LAYER_INTEGRAL = 200
SET_LAYER_INTEGRAL_UNTIL = 300
SET_BOUNDARY = 400
SET_DEPTH_INTEGRAL = 500
GET = 9000
LABDA_1 = 1
LABDA_2 = 2
COEFFICIENT = 3
COEFF2PARA = -1000
LAYERS = 4
DIFFUSION = 5
POROSITY = 6
ADSORPTION = 7
# INITIALIZE = 0
FLAG = 1
METHOD = 2
LAYER1 = 1
LAYER2 = 2
LAYER3 = 3
LAYER4 = 4
LAYER5 = 5
LAYER6 = 6
LAYER7 = 7
LAYER8 = 8
FOR_ALL_LAYERS = -1
NUMBER_OF_PROFILES = 12
NCOEFF = 22
NLAYER = 8





# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# MODULE: POM
#
# DESCRIPTION: Definition and allocation of parameters, scalars, and arrays used (mostly)
#              by the physical component of the system.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# SWITCH FOR PROGNOSTIC/DIAGNOSTIC SIMULATION
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   IDIAGN=0 PROGNOSTIC (T&S PROFILES COMPUTED)
#   IDIAGN=1 DIAGNOSTIC (T&S PROFILES PRESCRIBED)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
prognostic_diagnostic_mode_switch = float()

# # OF SURFACE/BOTTOM LAYERS WITH LOG DISTRIBUTION
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   KL1: SURFACE LAYERS
#   KL2: BOTTOM LAYERS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
KL1 = float(); KL2 = float()

# SWITCH FOR COLD/HOT START
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   IHOTST=0 "COLD" START FROM INITIAL CONDITION
#   IHOTST=1 "HOT" START FROM RESTART FILE
#
#   SEE SUBROUTINE "CALCDEPTH"
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
ihotst = float()

# MODEL TIME STEP
dti = float()

# LENGTH OF THE RUN (DAYS)
length_of_run_days = float()

# ITERATIONS NEEDED FOR AN "IDAYS" RUN
iterations_needed = float()

# COUNTER FOR THE TIME MARCHING LOOP
time_loop_counter = float()

# RUNNING TIME
TIME = float()

# TIME AT RESTART
TIME0 = float()

# BOTTOM DEPTH
bottom_depth = float()

# LATITUDE & LONGITUDE
alat = float(); longitude = float()

# CORIOLIS PARAMETER
coriolis_parameter = float()

# BACKGROUND DIFFUSION FOR U, V, Q2, Q2L
background_diffusion_momentum = float()

# BACKGROUND DIFFUSION FOR T,S & BFM TRACERS
background_diffusion_temperature = float(); background_diffusion_salinity = float(); background_diffusion_bfm = float()

# NUTRIENT RELAXATION TIME
NRT_disOxygen_IO_O = float()
NRT_phospate_IO_P = float()
NRT_nitrate_IO_N = float()
NRT_ammonium_IO_N = float()

# PARAMETER FOR THE HASSELIN FILTER
# (TIME STEP MIXING)
SMOTH = float()

# SPECIFIC HEAT TIMES RHO0
water_specific_heat_times_density = 4.187E6

# 1 DAY IN SECONDS (RECIPROCAL)
DAYI = one / seconds_per_day

# VERTICAL LAYERS
vertical_layers = 151

# FLAGS TO CHOOSE T,S AND BFM TRACERS SURFACE B.C. IN PROFTS
NBCT = float(); NBCS = float(); NBCBFM = float()

# FLAG TO CHOOSE JERLOV WATER TYPE IN PROFTS
NTP = float()

# T&S RELAXATION TIME (DAYS) FOR LATERAL ADVECTION
TRT = float(); SRT = float()

# DEPTH (m) AT WHICH LATERAL ADVECTION STARTS
upperH = float()

# RELAXATION TIME (DAYS) FOR SURFACE SALINITY FLUX
SSRT = float()

# VERTICAL COORDINATE SYSTEM
vertical_coordinates = np.empty(vertical_layers, dtype=float);  vertical_coordinates_staggered = np.empty(vertical_layers, dtype=float)
vertical_spacing = np.empty(vertical_layers, dtype=float); vertical_spacing_staggered = np.empty(vertical_layers, dtype=float); DZR = np.empty(vertical_layers, dtype=float)

# TEMPERATURE
temperature_forward = np.empty(vertical_layers, dtype=float); temperature = np.empty(vertical_layers, dtype=float); temperature_backward = np.empty(vertical_layers, dtype=float)

# SALINITY
salinity_forward = np.empty(vertical_layers, dtype=float); salinity = np.empty(vertical_layers, dtype=float); salinity_backward = np.empty(vertical_layers, dtype=float)

# DESITY
density_profile = np.empty(vertical_layers, dtype=float)

# VELOCITY
velocity_zonal_forward = np.empty(vertical_layers, dtype=float); velocity_zonal = np.empty(vertical_layers, dtype=float); UB = np.empty(vertical_layers, dtype=float)
velocity_meridional_forward = np.empty(vertical_layers, dtype=float); velocity_meridional = np.empty(vertical_layers, dtype=float); VB = np.empty(vertical_layers, dtype=float)

# TURBULENT KINETIC ENERGY (T.K.E.X2)
kinetic_energy_forward = np.empty(vertical_layers, dtype=float); kinetic_energy = np.empty(vertical_layers, dtype=float); kinetic_energy_backward = np.empty(vertical_layers, dtype=float)

# LENGTH SCALE
length_scale = np.empty(vertical_layers, dtype=float)

# (T.K.E.X2) TIMES LENGTH SCALE
kinetic_energy_times_length_forward = np.empty(vertical_layers, dtype=float); kinetic_energy_times_length = np.empty(vertical_layers, dtype=float); kinetic_energy_times_length_backward = np.empty(vertical_layers, dtype=float)

# VERTICAL DIFFUSION COEFFICIENTS
diffusion_coefficient_momentum = np.empty(vertical_layers, dtype=float); diffusion_coefficient_tracers = np.empty(vertical_layers, dtype=float); diffusion_coefficient_kinetic_energy = np.empty(vertical_layers, dtype=float)

# SERVICE ARRAYS USED IN POM ROUTINES
GM = np.empty(vertical_layers, dtype=float); GH = np.empty(vertical_layers, dtype=float)
SM = np.empty(vertical_layers, dtype=float); SH = np.empty(vertical_layers, dtype=float)
KN = np.empty(vertical_layers, dtype=float)
PROD = np.empty(vertical_layers, dtype=float); SPROD = np.empty(vertical_layers, dtype=float); BPROD = np.empty(vertical_layers, dtype=float)
A = np.empty(vertical_layers, dtype=float); C = np.empty(vertical_layers, dtype=float)
VH = np.empty(vertical_layers, dtype=float); VHP = np.empty(vertical_layers, dtype=float)
DTEF = np.empty(vertical_layers, dtype=float); D = np.empty(vertical_layers, dtype=float); DT = np.empty(vertical_layers, dtype=float)

# WIND STRESS
wind_stress_zonal = float(); wind_stress_meridional = float()

# BOTTOM STRESS
bottom_stress_zonal = float(); bottom_stress_meridional = float()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   N.B.
#   WHEN THE MODEL IS RUN IN DIAGNOSTIC MODE, THIS IS USED ONLY TO PROVIDE PAR TO THE BIOGEOCHEMICAL COMPONENT.
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
shortwave_radiation = float()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   N.B.
#   THE FOLLOWING SCALARS ARE USED ONLY WHEN THE MODEL IS RUN IN PROGNOSTIC MODE.
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# LOSS TERM OF THE SURFACE HEAT FLUX
# WTSURF = float()

# PRESCRIBED SURFACE T & S
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   TO BE USED (IF DESIRED) AS SURFACE
#   T & S SURFACE BOUNDARY CONDITION
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# TSURF = float(); SSURF = float()

# SURFACE SALINITY FLUX
# WSSURF = float()

# LATERAL ADVECTION FLUX FOR T & S
WTADV = np.empty(vertical_layers, dtype=float); WSADV = np.empty(vertical_layers, dtype=float)





# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# MODULE: SERVICE
#
# DESCRIPTION: List of Fortran parameters.
#              ONE-DIMENSIONAL BFM-POM SYSTEM
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# SURFACE NUTRIENTS
# PO4SURF = float(); NO3SURF = float(); NH4SURF = float();  SIO4SURF    = float()
# PO4BOTT = float(); NO3BOTT = float(); O2BOTT  = float();  PONBOTTgrad = float()

# SUSPENDED INORGANIC MATTER PROFILE
inorganic_suspended_matter = np.empty(vertical_layers - 1, dtype=float)
interpolated_w_velocity = np.empty((array_length, vertical_layers), dtype=float); w_eddy_velocity = np.empty(vertical_layers, dtype=float)

# FREQUENCY OF AVERAGING FOR OUTPUTi (IN DAYS)
savef = float(); nitend = float()
deltat = float()

# THESE ARE THE PATHWAYS FOR THE IC AND FORCING FILES (READ TROUGH NML)
wind_input = ''
surfaceS_input = ''
radiance_input = ''
ism_input = ''
Sal_input = ''
Temp_input = ''
W_input = ''
Weddy_input1 = ''
Weddy_input2 = ''
Sprofile_input = ''
Tprofile_input = ''
heat_input = ''
surfNut_input = ''
bottNut_input = ''
read_restart = ''

# THESE ARE THE PATHWAYS FOR THE BFM17 IC (READ TROUGH NML)
phyto_input = ''
zoop_input = ''
poc_input = ''
doc_input = ''
phos_input = ''
nit_input = ''
am_input = ''
oxy_input = ''





# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# MODULE: FORCING
#
# DESCRIPTION: DEFINITION AND ALLOCATION OF PARAMETERS, SCALAR AND ARRAYS USED
#              TO DEFINE THE PHYSICAL AND  BIOGEOCHEMICAL FORCING
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# MONTHLY SHORTWAVE RADIATION

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   N.B.!
#   ALWAYS NEEDED: WHEN THE MODEL IS RUN IN DIAGNOSTIC MODE PROVIDES ONLY PAR TO BFM.
#                  IN PROGNOSTIC CONTRIBUTES TO THE DEFINITION OF THE TEMPERATURE SURFACE BOUNDARY CONDITION.
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#
# SWRAD1 = float(); SWRAD2 = float()
# # SLUX1, SLUX2 = float()
#
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# #   N.B.!
# #   THE FOLLOWING SCALARS ARE USED ONLY WHEN THE MODEL IS RUN IN PROGNOSTIC MODE.
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # MONTHLY LOSS TERM OF THE SURFACE HEAT FLUX
# WTSURF1 = float(); WTSURF2 = float()
#
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# #   N.B.!
# #   THE FOLLOWING ARE ALWAYS USED.
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # PRESCRIBED T&S PROFILES
# TSTAR = np.empty(KB, dtype=float); SSTAR = np.empty(KB, dtype=float)
#
# # MONTHLY WIND STRESS
# WSU1 = float(); WSU2 = float()
# WSV1 = float(); WSV2 = float()
#
# # MONTHLY SURFACE SALINITY
# SSS1 = float(); SSS2 = float()
#
# # MONTHLY BOTTOM OXYGEN
# O2_b1 = float(); O2_b2 = float()
#
# # MONTHLY SURFACE AND BOTTOM NITRATE
# NO3_s1 = float(); NO3_s2 = float()
# NO3_b1 = float(); NO3_b2 = float()
#
# # MONTHLY SURFACE AND BOTTOM PHOSPHATE
# PO4_s1 = float(); PO4_s2 = float()
# PO4_b1 = float(); PO4_b2 = float()
#
# # MONTHLY SURFACE AMMONIA
# NH4_s1 = float(); NH4_s2 = float()
#
# # MONTHLY BOTTOM PON GRADIENT
# PON_b1 = float(); PON_b2 = float()
#
# # MONTHLY SURFACE SILICATE
# SIO4_s1 = float(); SIO4_s2 = float()
#
# # MONTHLY PROFILES OF INORGANIC SUSPENDED MATTER
# ISM1 = np.empty(KB, dtype=float); ISM2 = np.empty(KB, dtype=float)
#
# # MONTHLY PROFILES OF T & S
# TCLIM1 = np.empty(KB, dtype=float); TCLIM2 = np.empty(KB, dtype=float)
# SCLIM1 = np.empty(KB, dtype=float); SCLIM2 = np.empty(KB, dtype=float)
# WCLIM1 = np.empty(KB, dtype=float); WCLIM2 = np.empty(KB, dtype=float)
# WEDDY1 = np.empty(KB, dtype=float); WEDDY2 = np.empty(KB, dtype=float)
# WEDDY3 = np.empty(KB, dtype=float); WEDDY4 = np.empty(KB, dtype=float)
#
# SLUX1 = float(); QCORR1 = float(); QCORR2 = float()
#
# # INTERPOLATORS AND COUNTERS
# ICOUNTF = float(); IDOUNTF = float()
# IFCHGE = float();  IFDCHGE = float()
# IFDINT = float();  IFINT = float()
#
# RATIOF = float(); RATIOD = float()

from inputs.params_POMBFM import *
def forcing_manager(time_loop_counter):

    # LENGTH OF INPUT ARRAYS
    array_length = 13

    # LOOP COUNTER
    # K = float()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL ARRAYS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # RLENGTH = float()

    # INITIALISATION AND FIRST FORCING READING
    if time_loop_counter == 0:

        # DOUBLE READING OF DATA (NEEDED TO CARRY OUT THE TIME LINEAR INTERPOLATION)
        wind_speed_zonal1, wind_speed_meridional1, surface_salinity1, solar_radiation1, inorganic_suspended_matter1, \
            salinity_climatology1, temperature_climatology1, w_velocity_climatology1, w_eddy_velocity_1, \
            w_eddy_velocity_2, salinity_backward1, temperature_backward1, \
            shortwave_radiation1, surface_heat_flux1, kinetic_energy_loss1, \
            NO3_s1, NH4_s1, PO4_s1, SIO4_s1, O2_b1, NO3_b1, PO4_b1, PON_b1                          = read_pom_input()

        # NOT ALL VARIABLES FROM SECOND CALLING WILL BE USED
        wind_speed_zonal2, wind_speed_meridional2, surface_salinity2, solar_radiation2, inorganic_suspended_matter2, \
            salinity_climatology2, temperature_climatology2, w_velocity_climatology2, w_eddy_velocity_3, \
            w_eddy_velocity_4, salinity_backward2, temperature_backward2, \
            shortwave_radiation2, surface_heat_flux2, kinetic_energy_loss2, \
            NO3_s2, NH4_s2, PO4_s2, SIO4_s2, O2_b2, NO3_b2, PO4_b2, PON_b2                          = read_pom_input()

        # DAY COUNTER
        day_counter = 1

        # MONTH COUNTER
        month_counter = 1

        # TIME STEPS TO COVER ONE DAY
        timesteps_per_day = int(seconds_per_day) / int(dti)

        # TIME STEPS TO COVER ONE MONTH
        timesteps_per_month = 30 * timesteps_per_day

        # DAY INTERPOLATOR

        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # **                                                                   **
        # **  THE DAILY CLIMATOLOGICAL FORCING DATA ARE ASSUMED TO BE          **
        # **  CENTERED AT h 00.00 OF EACH CLIMATOLOGICAL DAY.THEREFORE         **
        # **  THE MONTH INTERPOLATOR(day_interpolator) IS INITIALISED AT THE VALUE       **
        # **  CORRESPONDING TO MIDNIGHT MINUS 1 TIMESTEP.                      **
        # **                                                                   **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **

        day_interpolator = -1

        # MONTH INTERPOLATOR

        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # **                                                                   **
        # **  THE MONTHLY CLIMATOLOGICAL FORCING DATA ARE ASSUMED TO BE        **
        # **  CENTERED AT DAY 15 OF EACH CLIMATOLOGICAL MONTH. THEREFORE       **
        # **  THE MONTH INTERPOLATOR (month_interpolator) IS INITIALISED AT THE VALUE       **
        # **  (timesteps_per_month/2)-1 CORRESPONDING TO DAY 15 MINUS 1 TIMESTEP.           **
        # **                                                                   **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **

        month_interpolator = (timesteps_per_month / 2.) - 1

        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # **                                                                   **
        # **  INITIAL READING OF THE MONTHLY FORCING                           **
        # **                                                                   **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **

        # WIND STRESS CONVERTED TO POM UNITS (N/m2-->m2/s2)
        wind_speed_zonal1[:] = -wind_speed_zonal1[:]  * 1.E-03
        wind_speed_zonal2[:]  = -wind_speed_zonal2[:]  * 1.E-03
        wind_speed_meridional1[:]  = -wind_speed_meridional1[:]  * 1.E-03
        wind_speed_meridional2[:]  = -wind_speed_meridional2[:]  * 1.E-03

        # HEAT FLUX CONVERTED TO POM UNITS(W/m2-->deg.C*m/s)
        shortwave_radiation1[:] = -shortwave_radiation1[:] / water_specific_heat_times_density
        shortwave_radiation2[:] = -shortwave_radiation2[:] / water_specific_heat_times_density
        surface_heat_flux1[:] = -surface_heat_flux1[:] / water_specific_heat_times_density
        surface_heat_flux2[:] = -surface_heat_flux2[:] / water_specific_heat_times_density

        # VERTICAL VELOCITY CONVERTED TO POM UNITS (m/s-->m/s)
        # w_velocity_climatology1[:] = w_velocity_climatology1[:]
        # w_velocity_climatology2[:] = w_velocity_climatology2[:]

        # UPDATE THE DAY COUNTER
        day_counter = day_counter + 1

        # UPDATE THE MONTH COUNTER
        month_counter = month_counter + 1

    # UPDATE INTERPOLATION COUNTERS
    day_interpolator = day_interpolator + 1
    ratio_day = day_interpolator / timesteps_per_day

    month_interpolator = month_interpolator + 1
    ratio_month = month_interpolator / timesteps_per_month

    # INTERPOLATE WIND STRESS
    wind_stress_zonal = wind_speed_zonal1 + ratio_month * (wind_speed_zonal2 - wind_speed_zonal1)
    wind_stress_meridional = wind_speed_meridional1 + ratio_month * (wind_speed_meridional2 - wind_speed_meridional1)

    # INTERPOLATE HEAT FLUX
    if prognostic_diagnostic_mode_switch == 0:
        surface_heat_flux_loss = surface_heat_flux1 + ratio_month * (surface_heat_flux2 - surface_heat_flux1)
        surface_solar_radiation = shortwave_radiation1 + ratio_month * (shortwave_radiation2 - shortwave_radiation1)
    elif prognostic_diagnostic_mode_switch == 1:
        # DAILY
        # surface_solar_radiation = shortwave_radiation1 + ratio_day * (shortwave_radiation2 - shortwave_radiation1)
        # MONTHLY
        surface_solar_radiation = shortwave_radiation1 + ratio_month * (shortwave_radiation2 - shortwave_radiation1)

    # INTERPOLATE T&S PROFILES
    interpolated_temperature = np.zeros((array_length, vertical_layers), dtype=float)
    interpolated_salinity = np.zeros((array_length, vertical_layers), dtype=float)
    interpolated_temperature[:] = temperature_climatology1[:] + \
                                           ratio_month * (temperature_climatology2[:] - temperature_climatology1[:])
    interpolated_salinity[:] = salinity_climatology1[:] + \
                                        ratio_month * (salinity_climatology2[:] - salinity_climatology1[:])
    interpolated_w_velocity[:]  = w_velocity_climatology1[:] + ratio_month * (w_velocity_climatology2[:] - w_velocity_climatology1[:])

    w_eddy_velocity = np.zeros((array_length, vertical_layers), dtype=float)
    if ratio_month <= 0.5:
        w_eddy_velocity[:] = w_eddy_velocity_1[:]
    else:
        w_eddy_velocity[:] = w_eddy_velocity_2[:]

    if prognostic_diagnostic_mode_switch == 0:
        surface_temperature = interpolated_temperature[0]
        surface_salinity = interpolated_salinity[0]
    elif prognostic_diagnostic_mode_switch == 1:
        temperature_forward[:] = interpolated_temperature[:]
        salinity_forward[:] = interpolated_salinity[:]

    # INTERPOLATE SUSPENDED INORGANIC MATTER
    inorganic_suspended_matter = np.zeros((array_length, vertical_layers), dtype=float)
    inorganic_suspended_matter[:] = inorganic_suspended_matter1[:] + ratio_month * (inorganic_suspended_matter2[:] - inorganic_suspended_matter1[:])

    # INTERPOLATE SURFACE NUTRIENTS
    NO3SURF = NO3_s1 + ratio_month * (NO3_s2 - NO3_s1)
    NH4SURF = NH4_s1 + ratio_month * (NH4_s2 - NH4_s1)
    PO4SURF = PO4_s1 + ratio_month * (PO4_s2 - PO4_s1)
    SIO4SURF = SIO4_s1 + ratio_month * (SIO4_s2 - SIO4_s1)

    # INTERPOLATE BOTTOM NUTRIENTS
    O2BOTT = O2_b1 + ratio_month * (O2_b2 - O2_b1)
    NO3BOTT = NO3_b1 + ratio_month * (NO3_b2 - NO3_b1)
    PO4BOTT = PO4_b1 + ratio_month * (PO4_b2 - PO4_b1)
    PONBOTTgrad = PON_b1 + ratio_month * (PON_b2 - PON_b1)

    if month_interpolator == timesteps_per_month:

        # A MONTH HAS GONE...IT IS NECESSARY TO...
        # ....UPDATE MONTH COUNTER....
        month_counter = month_counter + 1
        print('month_counter = ',month_counter)

        # ....RESET INTERPOLATOR....
        month_interpolator = 0

        # ....SHIFT THE MONTHLY DATA....
        wind_speed_zonal1 = wind_speed_zonal2
        wind_speed_meridional1 = wind_speed_meridional2
        shortwave_radiation1 = shortwave_radiation2
        surface_heat_flux1 = surface_heat_flux2
#         SLUX1 = SLUX2
        NO3_s1 = NO3_s2
        NH4_s1 = NH4_s2
        PO4_s1 = PO4_s2
        SIO4_s1 = SIO4_s2
        NO3_b1 = NO3_b2
        O2_b1 = O2_b2
        PO4_b1 = PO4_b2
        PON_b1 = PON_b2
        inorganic_suspended_matter1[:]     = inorganic_suspended_matter2[:]
        temperature_climatology1[:]   = temperature_climatology2[:]
        salinity_climatology1[:]   = salinity_climatology2[:]
        w_velocity_climatology1[:]   = w_velocity_climatology2[:]
        w_eddy_velocity_1[:]   = w_eddy_velocity_3[:]
        w_eddy_velocity_2[:]   = w_eddy_velocity_4[:]

        if month_counter > 13:
            # IF 12 MONTHS HAVE GONE, RESTART THE READING SEQUENCE
            month_counter = 2

            wind_speed_zonal1, wind_speed_meridional1, surface_salinity1, solar_radiation1, inorganic_suspended_matter1, \
                salinity_climatology1, temperature_climatology1, w_velocity_climatology1, w_eddy_velocity_1, \
                w_eddy_velocity_2, salinity_backward1, temperature_backward1, \
                shortwave_radiation1, surface_heat_flux1, kinetic_energy_loss1, \
                NO3_s1, NH4_s1, PO4_s1, SIO4_s1, O2_b1, NO3_b1, PO4_b1, PON_b1                       = read_pom_input()
            wind_speed_zonal1 = -wind_speed_zonal1 * 1.E-03
            wind_speed_meridional1 = -wind_speed_meridional1 * 1.E-03
            shortwave_radiation1 = -shortwave_radiation1 / water_specific_heat_times_density
            surface_heat_flux1 = -surface_heat_flux1 / water_specific_heat_times_density

            wind_speed_zonal2, wind_speed_meridional2, surface_salinity2, solar_radiation2, inorganic_suspended_matter2, \
                salinity_climatology2, temperature_climatology2, w_velocity_climatology2, w_eddy_velocity_3, \
                w_eddy_velocity_4, salinity_backward2, temperature_backward2, \
                shortwave_radiation2, surface_heat_flux2, kinetic_energy_loss2, \
                NO3_s2, NH4_s2, PO4_s2, SIO4_s2, O2_b2, NO3_b2, PO4_b2, PON_b2                       = read_pom_input()

            wind_speed_zonal2 = -wind_speed_zonal2 * 1.E-03
            wind_speed_meridional2 = -wind_speed_meridional2 * 1.E-03
            shortwave_radiation2 = -shortwave_radiation2 / water_specific_heat_times_density
            surface_heat_flux2 = -surface_heat_flux2 / water_specific_heat_times_density

    return




# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODULE: MEM
#
# DESCRIPTION: Definition of Global Shared Memory.
#
#              This module contains all the structural definitions of the BFM
#              and sets up the memory layout.
#              It is automatically generated from the prototype file
#              BFM/proto/ModuleMem.proto by including the information from
#              BFM/General/GlobalDefsBFM.model
#              Do not directly edit this code because changes will be lost at
#              any new compilation.
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   STATE VARIABLES INFO (PEL)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# !    3d name                                                  description           unit
# ! ---------- ------------------------------------------------------------ ---------------
# !        disOxygen_IO_O                                                       Oxygen      mmol O2/m3
# !        phospate_IO_P                                                    Phosphate       mmol P/m3
# !        nitrate_IO_N                                                      Nitrate       mmol N/m3
# !        ammonium_IO_N                                                     Ammonium       mmol N/m3
# !        O4n                                                 NitrogenSink       mmol N/m3
# !        silicate_IO_Si                                                     Silicate      mmol Si/m3
# !        reductEquiv_IO_R                                        Reduction Equivalents     mmol S--/m3
# !        pelBacteria_LO_C                               Aerobic and Anaerobic Bacteria         mg C/m3
# !        pelBacteria_LO_N                               Aerobic and Anaerobic Bacteria       mmol N/m3
# !        pelBacteria_LO_P                               Aerobic and Anaerobic Bacteria       mmol P/m3
# !        diatoms_LO_C                                                      Diatoms         mg C/m3
# !        diatoms_LO_N                                                      Diatoms       mmol N/m3
# !        diatoms_LO_P                                                      Diatoms       mmol P/m3
# !        diatoms_LO_Chl                                                      Diatoms       mg Chl/m3
# !        diatoms_LO_Si                                                      Diatoms      mmol Si/m3
# !        nanoflagellates_LO_C                                                  Flagellates         mg C/m3
# !        nanoflagellates_LO_N                                                  Flagellates       mmol N/m3
# !        nanoflagellates_LO_P                                                  Flagellates       mmol P/m3
# !        nanoflagellates_LO_Chl                                                  Flagellates       mg Chl/m3
# !        picophyto_LO_C                                            PicoPhytoplankton         mg C/m3
# !        picophyto_LO_N                                            PicoPhytoplankton       mmol N/m3
# !        picophyto_LO_P                                            PicoPhytoplankton       mmol P/m3
# !        picophyto_LO_Chl                                            PicoPhytoplankton       mg Chl/m3
# !        largephyto_LO_C                                          Large Phytoplankton         mg C/m3
# !        largephyto_LO_N                                          Large Phytoplankton       mmol N/m3
# !        largephyto_LO_P                                          Large Phytoplankton       mmol P/m3
# !        largephyto_LO_Chl                                          Large Phytoplankton       mg Chl/m3
# !        carnivMesozoo_LO_C                                  Carnivorous Mesozooplankton         mg C/m3
# !        carnivMesozoo_LO_N                                  Carnivorous Mesozooplankton       mmol N/m3
# !        carnivMesozoo_LO_P                                  Carnivorous Mesozooplankton       mmol P/m3
# !        omnivMesozoo_LO_C                                   Omnivorous Mesozooplankton         mg C/m3
# !        omnivMesozoo_LO_N                                   Omnivorous Mesozooplankton       mmol N/m3
# !        omnivMesozoo_LO_P                                   Omnivorous Mesozooplankton       mmol P/m3
# !        microzoo_LO_C                                             Microzooplankton         mg C/m3
# !        microzoo_LO_N                                             Microzooplankton       mmol N/m3
# !        microzoo_LO_P                                             Microzooplankton       mmol P/m3
# !        heteroFlagellates_LO_C                         Heterotrophic Nanoflagellates (HNAN)         mg C/m3
# !        heteroFlagellates_LO_N                         Heterotrophic Nanoflagellates (HNAN)       mmol N/m3
# !        heteroFlagellates_LO_P                         Heterotrophic Nanoflagellates (HNAN)       mmol P/m3
# !        labileDOM_NO_C                              Labile Dissolved Organic Matter         mg C/m3
# !        labileDOM_NO_N                              Labile Dissolved Organic Matter       mmol N/m3
# !        labileDOM_NO_P                              Labile Dissolved Organic Matter       mmol P/m3
# !        semilabileDOC_NO_C                         Semi-labile Dissolved Organic Carbon         mg C/m3
# !        semirefractDOC_NO_C                     Semi-refractory Dissolved Organic Carbon         mg C/m3
# !        particOrganDetritus_NO_C                                   Particulate Organic Matter         mg C/m3
# !        particOrganDetritus_NO_N                                   Particulate Organic Matter       mmol N/m3
# !        particOrganDetritus_NO_P                                   Particulate Organic Matter       mmol P/m3
# !        particOrganDetritus_NO_Si                                   Particulate Organic Matter      mmol Si/m3
# !        disInorgCarbon_IO_C                                   Dissolved Inorganic Carbon         mg C/m3
# !        totalAlkalinity_IO                                   Dissolved Inorganic Carbon      mmol eq/m3

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   STATE VARIABLES INFO (ICE)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   DEFINITION OF ARRAYS WHICH WILL HOLD ALL STATE VARIABLES AND OTHER
#   GLOBAL VARIABLES USED FOR EXCHANGE BETWEEN SUBMODELS AND/OR OUTPUT
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL SYSTEM CONSTANTS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
iiPel   = 0
iiIce   = 700
iiBen   = 1000
iiReset = -1000

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL SYSTEM CONSTANTS (PEL)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
no_d3_box_states   = 50
no_d3_box_diagnoss = 91
no_d2_box_diagnoss = 162
no_d3_box_flux     = 30
#
# try:
#     import INCLUDE_SEAICE
#     INCLUDE_SEAICE = True
# except FileNotFoundError:
#     INCLUDE_SEAICE = False

try:
    INCLUDE_SEAICE
except NameError:
    INCLUDE_SEAICE = False
else:
    INCLUDE_SEAICE = True
if INCLUDE_SEAICE:
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   GLOBAL SYSTEM CONSTANTS (ICE)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    no_d2_box_states_ice   = 0
    no_d2_box_diagnoss_ice = 0
    no_d2_box_flux_ice     = 0

# try:
#     import INCLUDE_BEN
#     INCLUDE_BEN = True
# except FileNotFoundError:
#     INCLUDE_BEN = False
try:
    INCLUDE_BEN
except NameError:
    INCLUDE_BEN = False
else:
    INCLUDE_BEN = True
if INCLUDE_BEN:
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   GLOBAL SYSTEM CONSTANTS (BEN)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    no_d2_box_states_ben   = 0
    no_d2_box_diagnoss_ben = 0
    no_d2_box_flux_ben     = 0

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL DEFINITION OF PELAGIC (D3/D2) STATE VARIABLES
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
ppdisOxygen_IO_O = 1;  ppphospate_IO_P = 2;  ppnitrate_IO_N = 3
ppammonium_IO_N = 4;  ppO4n = 5;  ppsilicate_IO_Si = 6;  ppreductEquiv_IO_R = 7;  pppelBacteria_LO_C = 8;  pppelBacteria_LO_N = 9;  pppelBacteria_LO_P = 10; ppdiatoms_LO_C = 11
ppdiatoms_LO_N = 12; ppdiatoms_LO_P = 13; ppdiatoms_LO_Chl = 14; ppdiatoms_LO_Si = 15; ppnanoflagellates_LO_C = 16; ppnanoflagellates_LO_N = 17; ppnanoflagellates_LO_P = 18
ppnanoflagellates_LO_Chl = 19; ppP2s = 0;  pppicophyto_LO_C = 20; pppicophyto_LO_N = 21; pppicophyto_LO_P = 22; pppicophyto_LO_Chl = 23; ppP3s = 0
pplargephyto_LO_C = 24; pplargephyto_LO_N = 25; pplargephyto_LO_P = 26; pplargephyto_LO_Chl = 27; ppP4s = 0;  ppcarnivMesozoo_LO_C = 28; ppcarnivMesozoo_LO_N = 29
ppcarnivMesozoo_LO_P = 30; ppomnivMesozoo_LO_C = 31; ppomnivMesozoo_LO_N = 32; ppomnivMesozoo_LO_P = 33; ppmicrozoo_LO_C = 34; ppmicrozoo_LO_N = 35; ppmicrozoo_LO_P = 36
ppheteroFlagellates_LO_C = 37; ppheteroFlagellates_LO_N = 38; ppheteroFlagellates_LO_P = 39; pplabileDOM_NO_C = 40; pplabileDOM_NO_N = 41; pplabileDOM_NO_P = 42; ppR1s = 0
ppsemilabileDOC_NO_C = 43; ppR2n = 0;  ppR2p = 0;  ppR2s = 0;  ppsemirefractDOC_NO_C = 44; ppR3n = 0;  ppR3p = 0; ppR3s = 0
ppparticOrganDetritus_NO_C = 45; ppparticOrganDetritus_NO_N = 46; ppparticOrganDetritus_NO_P = 47; ppparticOrganDetritus_NO_Si = 48; ppdisInorgCarbon_IO_C = 49; pptotalAlkalinity_IO = 50

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL DEFINITION OF SEAICE (D2) STATE VARIABLES
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL DEFINITION OF SEABEN (D2) STATE VARIABLES
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   CONSTITUENT PARAMETERS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
iiC = 1; iiN = 2; iiP = 3
iiL = 4; iiS = 5; iiH = 6
iiLastElement = 6

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GROUP PARAMETERS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
iiPelBacteria      = 1; iiB1 = 1
iiPhytoPlankton    = 4; iiP1 = 1; iiP2 = 2; iiP3 = 3; iiP4 = 4
iiMesoZooPlankton  = 2; iiZ3 = 1; iiZ4 = 2
iiMicroZooPlankton = 2; iiZ5 = 1; iiZ6 = 2
iiPelDetritus      = 4; iiR1 = 1; iiR2 = 2; iiR3 = 3; iiR6 = 4
iiInorganic        = 1; iiO3 = 1

# CalcPelBacteria(iiPelBacteria)           = True
# CalcPhytoPlankton(iiPhytoPlankton)       = True
# CalcMesoZooPlankton(iiMesoZooPlankton)   = True
# CalcMicroZooPlankton(iiMicroZooPlankton) = True
# CalcPelDetritus(iiPelDetritus)           = True
# CalcInorganic(iiInorganic)               = True

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GROUP PARAMETERS (ICE)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GROUP PARAMETERS (BEN)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL VARIABLES (PEL)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL VARIABLES (ICE)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL VARIABLES (BEN)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL DEFINITION OF PELAGIC (D3/D2) STATE VARIABLES WHICH CAN BE OUTPUTTED IN NETCDF
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# !    3d name                                                  description           unit
# ! ---------- ------------------------------------------------------------ ---------------
# !        ETW                                                  temperature               C
# !        ESW                                                     Salinity               -
# !       ERHO                                             Seawater Density           kg/m3
# !        EIR                                                   Irradiance         uE/m2/s
# !        ESS                                          Suspended Sediments            g/m3
# !       exud                                                    exudation         mg C/m3
# !      Depth                                              Gridpoint Depth               m
# !     Volume                                             Gridpoint Volume              m3
# !       Area                                               Gridpoint Area              m2
# !        DIC                                   Dissolved Inorganic Carbon         umol/kg
# !        CO2                                                      CO2(aq)         umol/kg
# !       pCO2                                                 Oceanic pCO2            uatm
# !       HCO3                                                  Bicarbonate         umol/kg
# !        CO3                                                    Carbonate         umol/kg
# !        ALK                                                   Alkalinity      umol eq/kg
# !         pH                                                           pH               -
# !      OCalc                                  Saturation state of Calcite               -
# !      OArag                                Saturation state of Aragonite               -
# !        EPR                                               Water Pressure            dbar
# !    totpelc                                        Total Mass in Pelagic             g C
# !    totpeln                                        Total Mass in Pelagic             g N
# !    totpelp                                        Total Mass in Pelagic             g P
# !    totpels                                        Total Mass in Pelagic            g Si
# !      cxoO2                                            Oxygen Saturation      mmol O2/m3
# !     eO2mO2                                   Relative Oxygen saturation               -
# !       Chla                                                Chlorophyll-a       mg Chl/m3
# !    flPTreductEquiv_IO_R                        Pelagic Anaerobic Mineralization Rate    mmol O2/m3/d
# !    flN3O4n                                 Pelagic Denitrification Rate     mmol N/m3/d
# !    flammonium_IO_Nitrate_IO_N                                   Pelagic Nitrification Rate     mmol N/m3/d
# !     sediR2                                  Detritus sedimentation rate             m/d
# !     sediR6                                  Detritus sedimentation rate             m/d
# !       xEPS                                 Total Extinction Coefficient             1/m
# !   ABIO_eps                               Abiotic Extinction Coefficient             1/m
#
# ! qpcPPY(iiP1)                                                      Diatoms     mmol P/mg C
# ! qpcPPY(iiP2)                                                  Flagellates     mmol P/mg C
# ! qpcPPY(iiP3)                                            PicoPhytoplankton     mmol P/mg C
# ! qpcPPY(iiP4)                                          Large Phytoplankton     mmol P/mg C
# ! qncPPY(iiP1)                                                      Diatoms     mmol N/mg C
# ! qncPPY(iiP2)                                                  Flagellates     mmol N/mg C
# ! qncPPY(iiP3)                                            PicoPhytoplankton     mmol N/mg C
# ! qncPPY(iiP4)                                          Large Phytoplankton     mmol N/mg C
# ! qscPPY(iiP1)                                                      Diatoms    mmol Si/mg C
# ! qscPPY(iiP2)                                                  Flagellates    mmol Si/mg C
# ! qscPPY(iiP3)                                            PicoPhytoplankton    mmol Si/mg C
# ! qscPPY(iiP4)                                          Large Phytoplankton    mmol Si/mg C
# ! qlcPPY(iiP1)                                                      Diatoms    mg Chl /mg C
# ! qlcPPY(iiP2)                                                  Flagellates    mg Chl /mg C
# ! qlcPPY(iiP3)                                            PicoPhytoplankton    mg Chl /mg C
# ! qlcPPY(iiP4)                                          Large Phytoplankton    mg Chl /mg C
# ! qpcMEZ(iiZ3)                                  Carnivorous Mesozooplankton     mmol P/mg C
# ! qpcMEZ(iiZ4)                                   Omnivorous Mesozooplankton     mmol P/mg C
# ! qncMEZ(iiZ3)                                  Carnivorous Mesozooplankton     mmol N/mg C
# ! qncMEZ(iiZ4)                                   Omnivorous Mesozooplankton     mmol N/mg C
# ! qpcMIZ(iiZ5)                                             Microzooplankton     mmol P/mg C
# ! qpcMIZ(iiZ6)                         Heterotrophic Nanoflagellates (HNAN)     mmol P/mg C
# ! qncMIZ(iiZ5)                                             Microzooplankton     mmol N/mg C
# ! qncMIZ(iiZ6)                         Heterotrophic Nanoflagellates (HNAN)     mmol N/mg C
# ! qpcOMT(iiR1)                              Labile Dissolved Organic Matter     mmol N/mg C
# ! qpcOMT(iiR2)                         Semi-labile Dissolved Organic Carbon     mmol N/mg C
# ! qpcOMT(iiR3)                     Semi-refractory Dissolved Organic Carbon     mmol N/mg C
# ! qpcOMT(iiR6)                                   Particulate Organic Matter     mmol N/mg C
# ! qncOMT(iiR1)                              Labile Dissolved Organic Matter     mmol P/mg C
# ! qncOMT(iiR2)                         Semi-labile Dissolved Organic Carbon     mmol P/mg C
# ! qncOMT(iiR3)                     Semi-refractory Dissolved Organic Carbon     mmol P/mg C
# ! qncOMT(iiR6)                                   Particulate Organic Matter     mmol P/mg C
# ! qscOMT(iiR1)                              Labile Dissolved Organic Matter    mmol Si/mg C
# ! qscOMT(iiR2)                         Semi-labile Dissolved Organic Carbon    mmol Si/mg C
# ! qscOMT(iiR3)                     Semi-refractory Dissolved Organic Carbon    mmol Si/mg C
# ! qscOMT(iiR6)                                   Particulate Organic Matter    mmol Si/mg C
# ! qpcPBA(iiB1)                               Aerobic and Anaerobic Bacteria     mmol P/mg C
# ! qncPBA(iiB1)                               Aerobic and Anaerobic Bacteria     mmol N/mg C
# ! sediPPY(iiP1)                                                      Diatoms             m/d
# ! sediPPY(iiP2)                                                  Flagellates             m/d
# ! sediPPY(iiP3)                                            PicoPhytoplankton             m/d
# ! sediPPY(iiP4)                                          Large Phytoplankton             m/d
# ! sediMIZ(iiZ5)                                             Microzooplankton             m/d
# ! sediMIZ(iiZ6)                         Heterotrophic Nanoflagellates (HNAN)             m/d
# ! sediMEZ(iiZ3)                                  Carnivorous Mesozooplankton             m/d
# ! sediMEZ(iiZ4)                                   Omnivorous Mesozooplankton             m/d
# ! sunPPY(iiP1)                                                      Diatoms             1/d
# ! sunPPY(iiP2)                                                  Flagellates             1/d
# ! sunPPY(iiP3)                                            PicoPhytoplankton             1/d
# ! sunPPY(iiP4)                                          Large Phytoplankton             1/d
# ! eiPPY(iiP1)                                                      Diatoms               -
# ! eiPPY(iiP2)                                                  Flagellates               -
# ! eiPPY(iiP3)                                            PicoPhytoplankton               -
# ! eiPPY(iiP4)                                          Large Phytoplankton               -
# ! ELiPPY(iiP1)                                                      Diatoms            W/m2
# ! ELiPPY(iiP2)                                                  Flagellates            W/m2
# ! ELiPPY(iiP3)                                            PicoPhytoplankton            W/m2
# ! ELiPPY(iiP4)                                          Large Phytoplankton            W/m2

ppETW     = 1;  ppESW     = 2;  ppERHO    = 3;  ppEIR = 4
ppESS     = 5;  ppexud    = 6;  ppDepth   = 7;  ppVolume  = 8;  ppArea     = 9;  ppDIC   = 10; ppCO2   = 11
pppCO2    = 12; ppHCO3    = 13; ppCO3     = 14; ppALK     = 15; pppH       = 16; ppOCalc = 17; ppOArag = 18
ppEPR     = 19; pptotpelc = 20; pptotpeln = 21; pptotpelp = 22; pptotpels  = 23
ppcxoO2   = 24; ppeO2mO2  = 25; ppChla    = 26; ppflPTreductEquiv_IO_R = 27; ppflN3O4n  = 28
ppflammonium_IO_Nitrate_IO_N = 29; ppsediR2  = 30; ppsediR6  = 31; ppxEPS    = 32; ppABIO_eps = 33

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# !    2d name                                                  description           unit
# ! ---------- ------------------------------------------------------------ ---------------
# !      ETAUB                                                Bottom Stress           N m/s
# !   EPCO2air                             Atmospheric CO2 Partial Pressure            uatm
# ! CO2airflux                                             Sea-air CO2 Flux       mmol/m2/d
# !     Area2d                                           2-D Gridpoint Area              m2
# ! ThereIsLight                                   Switch for day/night cycle               -
# !       SUNQ                                           Daylength in hours               h
# !      EWIND                                                   Wind speed             m/s
# !    totsysc                                                   total mass             g C
# !    totsysn                                                   total mass             g N
# !    totsysp                                                   total mass             g P
# !    totsyss                                                   total mass            g Si
# !       EICE                                             Sea-ice fraction               -

ppETAUB   = 1; ppEPCO2air     = 2;  ppCO2airflux = 3
ppArea2d  = 4; ppThereIsLight = 5;  ppSUNQ       = 6;  ppEWIND = 7; pptotsysc = 8
pptotsysn = 9; pptotsysp      = 10; pptotsyss    = 11; ppEICE  = 12

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL DEFINITION OF SEAICE (D2) STATE VARIABLES WHICH CAN BE OUTPUTTED IN NETCDF
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL DEFINITION OF BENTHIC (D2) STATE VARIABLES WHICH CAN BE OUTPUTTED IN NETCDF
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   BOUNDARY FLUXES
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   OHTER 3D-GLOBAL VARIABLES
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   OHTER 2D-GLOBAL VARIABLES (PEL)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   OHTER 2D-GLOBAL VARIABLES (PEL)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   SHARED GLOBAL FUNCTIONS (PEL) (MUST BE BELOW CONTAINS)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# flux, flux_vector, Source, Source_D3_vector, \
#     fixed_quota_flux_vector
# ppPelBacteria, PelBacteria, ppPhytoPlankton, PhytoPlankton, \
#     ppMesoZooPlankton, MesoZooPlankton, ppMicroZooPlankton, MicroZooPlankton, \
#     ppPelDetritus, PelDetritus, ppInorganic, Inorganic


# try:
#     import INCLUDE_SEAICE
#     INCLUDE_SEAICE = True
# except FileNotFoundError:
#     INCLUDE_SEAICE = False
try:
    INCLUDE_SEAICE
except NameError:
    INCLUDE_SEAICE = False
else:
    INCLUDE_SEAICE = True
if INCLUDE_SEAICE:
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SHARED GLOBAL FUNCTIONS (ICE) (MUST BE BELOW CONTAINS)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Source_D2_vector_ice
    pass

# try:
#     import INCLUDE_BEN
#     INCLUDE_BEN = True
# except FileNotFoundError:
#     INCLUDE_BEN = False
try:
    INCLUDE_BEN
except NameError:
    INCLUDE_BEN = False
else:
    INCLUDE_BEN = True
if INCLUDE_BEN:
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SHARED GLOBAL FUNCTIONS (ICE) (MUST BE BELOW CONTAINS)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Source_D2_vector_ben
    pass

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GROUP PELAGIC (D3) STATE FUNCTIONS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# FUNCTIONS DEFINED IN MODULEMEM
def ppPelBacteria(n,constituent):

    pointers = [pppelBacteria_LO_C, pppelBacteria_LO_N, pppelBacteria_LO_P, 0., 0., 0.]

    if n > 1 or n == 0:
        ppPelBacteria = 0
    elif constituent > 1 or constituent == 0:
        ppPelBacteria = 0
    else:
        ppPelBacteria = pointers[(n-1)*iiLastElement + constituent]

    return ppPelBacteria


def ppPhytoPlankton(n,constituent):

    pointers = [ppdiatoms_LO_C, ppdiatoms_LO_N, ppdiatoms_LO_P, ppdiatoms_LO_Chl, ppdiatoms_LO_Si, 0.,
                ppnanoflagellates_LO_C, ppnanoflagellates_LO_N, ppnanoflagellates_LO_P, ppnanoflagellates_LO_Chl, 0.   , 0.,
                pppicophyto_LO_C, pppicophyto_LO_N, pppicophyto_LO_P, pppicophyto_LO_Chl, 0.   , 0.,
                pplargephyto_LO_C, pplargephyto_LO_N, pplargephyto_LO_P, pplargephyto_LO_Chl, 0.   , 0.]

    if n > 4 or n == 0:
        ppPhytoPlankton = 0
    elif constituent > iiLastElement or constituent == 0:
        ppPhytoPlankton = 0
    else:
        ppPhytoPlankton = pointers[(n-1)*iiLastElement + constituent]

    return ppPhytoPlankton


def ppMesoZooPlankton(n,constituent):

    pointers = [ppcarnivMesozoo_LO_C, ppcarnivMesozoo_LO_N, ppcarnivMesozoo_LO_P, 0., 0., 0.,
                ppomnivMesozoo_LO_C, ppomnivMesozoo_LO_N, ppomnivMesozoo_LO_P, 0., 0., 0.]

    if n > 2 or n == 0:
        ppMesoZooPlankton = 0
    elif constituent > iiLastElement or constituent == 0:
        ppMesoZooPlankton = 0
    else:
        ppMesoZooPlankton = pointers[(n-1)*iiLastElement + constituent]

    return ppMesoZooPlankton


def ppMicroZooPlankton(n,constituent):

    pointers = [ppmicrozoo_LO_C, ppmicrozoo_LO_N, ppmicrozoo_LO_P, 0., 0., 0.,
                ppheteroFlagellates_LO_C, ppheteroFlagellates_LO_N, ppheteroFlagellates_LO_P, 0., 0., 0.]

    if n > 2 or n == 0:
        ppMicroZooPlankton = 0
    elif constituent > iiLastElement or constituent == 0:
        ppMicroZooPlankton = 0
    else:
        ppMicroZooPlankton = pointers[(n-1)*iiLastElement + constituent]

    return ppMicroZooPlankton


def ppPelDetritus(n,constituent):

    pointers = [pplabileDOM_NO_C, pplabileDOM_NO_N, pplabileDOM_NO_P, 0.   , 0.   , 0.,
                ppsemilabileDOC_NO_C, 0.   , 0.   , 0.   , 0.   , 0.,
                ppsemirefractDOC_NO_C, 0.   , 0.   , 0.   , 0.   , 0.,
                ppparticOrganDetritus_NO_C, ppparticOrganDetritus_NO_N, ppparticOrganDetritus_NO_P, 0.   , ppparticOrganDetritus_NO_Si, 0.]

    if n > 4 or n == 0:
        ppPelDetritus = 0
    elif constituent > iiLastElement or constituent == 0:
        ppPelDetritus = 0
    else:
        ppPelDetritus = pointers[(n-1)*iiLastElement + constituent]

    return ppPelDetritus















# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# MODULE: PARAM
#
# DESCRIPTION: List of global model parameters.
#              (global variables that can be changed during the model initialization
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ! Global Switches : turn on/off or choose model components
# ! NAME                          KIND    DESCRIPTION
# ! CalcPelagicFlag               logical Pelagic System
# ! CalcBenthicFlag               numeric Benthic system
# !                                       0 = No Benthic System
# !                                       The following are Not Yet Activated
# !                                       1 = Simple Benthic Return
# !                                       2 = Benthic organisms and intermediate
# !                                           complexity nutrient regeneration
# !                                       3 = Benthic organisms and full nutrient
# !                                           regeneration (early diagenesis)
# ! CalcTransportFlag             logical Compute Transport Term (when coupled
# !                                       with a OGCM)
# ! CalcConservationFlag          logical Mass Conservation Check
# ! CalcPhytoPlankton             logical Pelagic Phytoplankton (vector)
# ! CalcPelBacteria               logical Pelagic Bacteria (vector)
# ! CalcMesoZooPlankton           logical Mesozooplankton (vector)
# ! CalcMicroZooPlankton          logical Microzooplankton (vector)
# ! CalcPelChemistry              logical Pelagic Hydrochemical Processes
# ! AssignPelBenFluxesInBFMFlag   logical Benthic-pelagic fluxes are added to the
# !                                       time integration
# ! AssignAirPelFluxesInBFMFlag   logical Air-sea fluxes are added to the
# !                                       time integration
# ! ChlDynamicsFlag               numeric Choose the dynamics of Chl-a
# !                                       1 = diagnostic, optimal light property
# !                                           in phytoplankton
# !                                           (Ebenhoeh et al 1995, ERSEM-II)
# !                                       2 = state variable, constituent of
# !                                           phytoplankton
# ! LightPeriodFlag               numeric Choose the light averaging period
# !                                       1 = Instantanous irradiance
# !                                       2 = Daily average
# !                                       3 = Daylight average with explicit
# !                                           photoperiod
# ! LightLocationFlag             numeric Choose the parameterization of light
# !                                       location in the discrete grid
# !                                       1 = Light at the top of the cell
# !                                       2 = Light in the middle of the cell
# !                                       3 = Average Light in the cell
# ! check_fixed_quota             numeric Check whether zooplankton have fixed quota
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
CalcPelagicFlag = True
CalcBenthicFlag = 0  # Switch for Benthic system
CalcSeaiceFlag  = True  # Switch for Seaice system

CalcTransportFlag           = False
CalcConservationFlag        = True
CalcPelChemistry            = True
AssignPelBenFluxesInBFMFlag = True
AssignAirPelFluxesInBFMFlag = True

ChlDynamicsFlag   = 2
LightPeriodFlag   = 1
LightLocationFlag = 3
check_fixed_quota = 0

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ! Global Parameters : used throughout the model and not related
# !                     to a specific component
# ! NAME          UNIT          DESCRIPTION
# ! p_small      [-]           Smallest numeric value (the model "zero")
# ! slp0         [mbar]        Reference sea level pressure
# ! p_PAR        [-]           Fraction of Photosynthetically Available Radiation
# ! p_eps0       [1/m]         Background extinction coefficient
# ! p_epsESS     [m2/g]        Specific attenuation coefficient of
# !                            suspended sediments
# ! p_epsChla   [m2/mgChla]    Chla-specific extinction coefficient
# ! p_epsR6      [m2/mgC]      Specific attenuation coefficient of particulate
# !                            detritus
# ! p_pe_labileDOM_NO_C     [-]           Fractional content of C in cytoplasm
# ! p_pe_labileDOM_NO_N     [-]           Fractional content of N in cytoplasm
# ! p_pe_labileDOM_NO_P     [-]           Fractional content of P in cytoplasm
# ! p_qro        [mmolHS-/     Stoichiometric coefficient for
# !               mmolO2]      anaerobic reactions
# ! p_qon_dentri [mmolO2/      Stoichiometric coefficient for
# !               mmolN]       denitrification
# ! p_qon_nitri  [mmolO2/      Stoichiometric coefficient for
# !               mmolN]       nitrification (3/2)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   PELAGIC MODEL PARAMETERS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
p_small      = 1.0E-20
slp0         = 1013.25
p_PAR        = 0.50
p_eps0       = 0.04
p_epsESS     = 0.04E-03
p_epsChla    = 0.03
p_epsR6      = 0.1E-03
p_pe_labileDOM_NO_C     = 0.60
p_pe_labileDOM_NO_N     = 0.72
p_pe_labileDOM_NO_P     = 0.832
p_pe_R1s     = 0.06
p_qro        = 0.50
p_qon_dentri = 1.25
p_qon_nitri  = 1.5

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ! Benthic model parameters
# ! NAME          UNIT          DESCRIPTION
# ! p_sedlevels   [-]           Number of sigma levels for benthic nutrient
# ! p_poro0       [-]           Constant porosity for 0D and 1D runs
# ! p_InitSink    Logical       parameter to Initialize BenthicSInk var.
# ! p_q10diff     [-]           Temperature-dependency porewater diffusion
# ! p_clDxm       [m]           minimal value of D?.m for calculation of the alpha
# ! p_d_tot       [m]           Thickness of modelled benthic sediment layers
# ! p_clD1D2m     [m]           minimum distance between D1m and D2m
# ! p_d_tot_2     [m]           maximal Thickness of D2m
# ! p_sedsigma    [-]           Parameter for sigma level distribution
# ! p_poro        [-]           Sediment porosity
# ! p_p_ae        [-]           Adsorption coefficient
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# 0D-PARAMETERS
p_sedlevels = 20
p_sedsigma  = 2.0
p_d_tot     = 0.30
p_poro0     = 0.4

# 1D-PARAMETERS
# try:
#     import INCLUDE_BEN
#     INCLUDE_BEN = True
# except FileNotFoundError:
#     INCLUDE_BEN = False
try:
    INCLUDE_BEN
except NameError:
    INCLUDE_BEN = False
else:
    INCLUDE_BEN = True
if INCLUDE_BEN:
    p_InitSink = 100.0
    p_q10diff  = 1.49
    p_clDxm    = 0.001
    p_clD1D2m  = 0.01
    p_d_tot_2  = 0.35

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   SEAICE MODEL PARAMETERS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
















# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=






# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# SUBROUTINE: path_error
#
# DESCRIPTION: Print error message if file can't be found.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def path_error(file):
    message = 'Unable to locate:    ' + file
    print(message)





# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# SUBROUTINE: error_msg_prn
#
# DESCRIPTION: Print bfm error messages.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def error_msg_prn(code,infile,what):

    f = open("LOGUNIT","w")
    f.write("*********** RUN TIME ERROR BEGIN ***********")

    if code == ALLOC:
        f.write("Unable to allocate " + what + " in " + infile)
    elif code == NML_OPEN:
        f.write("Unable to open " + what + " in " + infile)
    elif code == NML_READ:
        f.write("Namelist mismatch in " + what + " opened by " + infile)
    elif code == DIM_MISMATCH:
        f.write("Dimension mismatch while reading " + what + " in " + infile)

    f.write("***********  RUN TIME ERROR END  ***********")
    print("BFM error (see logfile)")


#   EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

