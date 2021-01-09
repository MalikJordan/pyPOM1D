# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
import numpy as np
from pom.initialize_variables import read_pom_input

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
bef_lwp = True
ZERO = float(0.0)
ONE = float(1.0)
PI = float(3.14159265359)
BASETEMP = float(20.0)

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

MIN_VAL_EXPFUN = ZERO
ECOLOGY = 1  # BASE TEMPERATURE FOR Q10
TRANSPORT = 2
ZERO_KELVIN = -273.15
Rgas = 83.131  # GAS CONSTANT: bar mol^-1 deg-1
MW_C = 12.0  # MOLECULAR WEIGHT CARBON
MW_N = 14.0  # MOLECULAR WEIGHT NITROGEN
MW_P = 31.0  # MOLECULAR WEIGHT PHOSPHORUS
MW_SI = 28.0  # MOLECULAR WEIGHT SILICA
E2W = 0.217  # MOLECULAR WEIGHT CONVERSION FACTOR EINSTEIN->W
SEC_PER_DAY = 86400.  # SECONDS IN DAY
DAY_PER_SEC = ZERO  # INVERSE OF SECONDS IN DAY
ONE_PER_DAY = 1.  # RATE WHICH IS USED IN CASES WHERE IMPLICITLY ASSUMED

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
IDIAGN = float()

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
IHOTST = float()

# MODEL TIME STEP
DTI = float()

# LENGTH OF THE RUN (DAYS)
IDAYS = float()

# ITERATIONS NEEDED FOR AN "IDAYS" RUN
IEND = float()

# COUNTER FOR THE TIME MARCHING LOOP
intt = float()

# RUNNING TIME
TIME = float()

# TIME AT RESTART
TIME0 = float()

# BOTTOM DEPTH
H = float()

# LATITUDE & LONGITUDE
ALAT = float(); ALON = float()

# CORIOLIS PARAMETER
COR = float()

# BACKGROUND DIFFUSION FOR U, V, Q2, Q2L
UMOL = float()

# BACKGROUND DIFFUSION FOR T,S & BFM TRACERS
UMOLT = float(); UMOLS = float(); UMOLBFM = float()

# NUTRIENT RELAXATION TIME
NRT_o2o = float()
NRT_n1p = float()
NRT_n3n = float()
NRT_n4n = float()

# PARAMETER FOR THE HASSELIN FILTER
# (TIME STEP MIXING)
SMOTH = float()

# SPECIFIC HEAT TIMES RHO0
RCP = 4.187E6

# 1 DAY IN SECONDS (RECIPROCAL)
DAYI = ONE/SEC_PER_DAY

# VERTICAL LAYERS
KB = 151

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
Z = np.empty(KB,dtype=float);  ZZ = np.empty(KB,dtype=float)
DZ = np.empty(KB,dtype=float); DZZ = np.empty(KB,dtype=float); DZR = np.empty(KB,dtype=float)

# TEMPERATURE
TF = np.empty(KB,dtype=float); T = np.empty(KB,dtype=float); TB = np.empty(KB,dtype=float)

# SALINITY
SF = np.empty(KB,dtype=float); S = np.empty(KB,dtype=float); SB = np.empty(KB,dtype=float)

# DESITY
RHO = np.empty(KB,dtype=float)

# VELOCITY
UF = np.empty(KB,dtype=float); U = np.empty(KB,dtype=float); UB = np.empty(KB,dtype=float)
VF = np.empty(KB,dtype=float); V = np.empty(KB,dtype=float); VB = np.empty(KB,dtype=float)

# TURBULENT KINETIC ENERGY (T.K.E.X2)
Q2F = np.empty(KB,dtype=float); Q2 = np.empty(KB,dtype=float); Q2B = np.empty(KB,dtype=float)

# LENGTH SCALE
L = np.empty(KB,dtype=float)

# (T.K.E.X2) TIMES LENGTH SCALE
Q2LF = np.empty(KB,dtype=float); Q2L = np.empty(KB,dtype=float); Q2LB = np.empty(KB,dtype=float)

# VERTICAL DIFFUSION COEFFICIENTS
KM = np.empty(KB,dtype=float); KH = np.empty(KB,dtype=float); KQ = np.empty(KB,dtype=float)

# SERVICE ARRAYS USED IN POM ROUTINES
GM = np.empty(KB,dtype=float); GH = np.empty(KB,dtype=float)
SM = np.empty(KB,dtype=float); SH = np.empty(KB,dtype=float)
KN = np.empty(KB,dtype=float)
PROD = np.empty(KB,dtype=float); SPROD = np.empty(KB,dtype=float); BPROD = np.empty(KB,dtype=float)
A = np.empty(KB,dtype=float); C = np.empty(KB,dtype=float)
VH = np.empty(KB,dtype=float); VHP = np.empty(KB,dtype=float)
DTEF = np.empty(KB,dtype=float); D = np.empty(KB,dtype=float); DT = np.empty(KB,dtype=float)

# WIND STRESS
WUSURF = float(); WVSURF = float()

# BOTTOM STRESS
WUBOT = float(); WVBOT = float()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   N.B.
#   WHEN THE MODEL IS RUN IN DIAGNOSTIC MODE, THIS IS USED ONLY TO PROVIDE PAR TO THE BIOGEOCHEMICAL COMPONENT.
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# SWRAD = float()

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
WTADV = np.empty(KB,dtype=float); WSADV = np.empty(KB,dtype=float)





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
ISM = np.empty(KB-1, dtype=float)
WGEN = np.empty(KB, dtype=float); WEDDY = np.empty(KB, dtype=float)

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


def FORCING_MANAGER():

    # LENGTH OF INPUT ARRAYS
    array_length = 13

    # LOOP COUNTER
    K = float()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL ARRAYS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    RLENGTH = float()

    # INITIALISATION AND FIRST FORCING READING
    if intt == int(ONE):

        # DOUBLE READING OF DATA (NEEDED TO CARRY OUT THE TIME LINEAR INTERPOLATION)
        WSU1, WSV1, SSS1, SLUX1, ISM1, SCLIM1, TCLIM1, WCLIM1, WEDDY1, WEDDY2, SB, TB, \
            SWRAD1, WTSURF1, QCORR1, NO3_s1, NH4_s1, PO4_s1, SIO4_s1, O2_b1, NO3_b1, PO4_b1, PON_b1 = read_pom_input()

        # NOT ALL VARIABLES FROM SECOND CALLING WILL BE USED
        WSU2, WSV2, SSS2, SLUX2, ISM2, SCLIM2, TCLIM2, WCLIM2, WEDDY3, WEDDY4, SB2, TB2, \
            SWRAD2, WTSURF2, QCORR2, NO3_s2, NH4_s2, PO4_s2, SIO4_s2, O2_b2, NO3_b2, PO4_b2, PON_b2 = read_pom_input()

        # DAY COUNTER
        IDOUNTF = 1

        # MONTH COUNTER
        ICOUNTF = 1

        # TIME STEPS TO COVER ONE DAY
        IFDCHGE = int(SEC_PER_DAY)/int(DTI)

        # TIME STEPS TO COVER ONE MONTH
        IFCHGE = float(30)*IFDCHGE

        # DAY INTERPOLATOR

        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # **                                                                   **
        # **  THE DAILY CLIMATOLOGICAL FORCING DATA ARE ASSUMED TO BE          **
        # **  CENTERED AT h 00.00 OF EACH CLIMATOLOGICAL DAY.THEREFORE         **
        # **  THE MONTH INTERPOLATOR(IFDINT) IS INITIALISED AT THE VALUE       **
        # **  CORRESPONDING TO MIDNIGHT MINUS 1 TIMESTEP.                      **
        # **                                                                   **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **

        IFDINT = -int(ONE)

        # MONTH INTERPOLATOR

        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # **                                                                   **
        # **  THE MONTHLY CLIMATOLOGICAL FORCING DATA ARE ASSUMED TO BE        **
        # **  CENTERED AT DAY 15 OF EACH CLIMATOLOGICAL MONTH. THEREFORE       **
        # **  THE MONTH INTERPOLATOR (IFINT) IS INITIALISED AT THE VALUE       **
        # **  (IFCHGE/2)-1 CORRESPONDING TO DAY 15 MINUS 1 TIMESTEP.           **
        # **                                                                   **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **

        IFINT = (IFCHGE / float(2)) - int(ONE)

        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # **                                                                   **
        # **  INITIAL READING OF THE MONTHLY FORCING                           **
        # **                                                                   **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **

        # WIND STRESS CONVERTED TO POM UNITS (N/m2-->m2/s2)
        WSU1 = -WSU1 * float(1.E-03)
        WSU2 = -WSU2 * float(1.E-03)
        WSV1 = -WSV1 * float(1.E-03)
        WSV2 = -WSV2 * float(1.E-03)

        # HEAT FLUX CONVERTED TO POM UNITS(W/m2-->deg.C*m/s)
        SWRAD1 = -SWRAD1 / RCP
        SWRAD2 = -SWRAD2 / RCP
        WTSURF1 = -WTSURF1 / RCP
        WTSURF2 = -WTSURF2 / RCP

        # VERTICAL VELOCITY CONVERTED TO POM UNITS (m/s-->m/s)
        WCLIM1 = WCLIM1
        WCLIM2 = WCLIM2

        # UPDATE THE DAY COUNTER
        IDOUNTF = IDOUNTF + int(ONE)

        # UPDATE THE MONTH COUNTER
        ICOUNTF = ICOUNTF + int(ONE)

    # UPDATE INTERPOLATION COUNTERS
    IFDINT = IFDINT + int(ONE)
    RATIOD = float(IFDINT) / float(IFDCHGE)

    IFINT = IFINT + int(ONE)
    RATIOF = float(IFINT) / float(IFCHGE)

    # INTERPOLATE WIND STRESS
    WUSURF = WSU1 + RATIOF * (WSU2 - WSU1)
    WVSURF = WSV1 + RATIOF * (WSV2 - WSV1)

    # INTERPOLATE HEAT FLUX
    if IDIAGN == int(ZERO):
        WTSURF = WTSURF1 + RATIOF * (WTSURF2 - WTSURF1)
        SWRAD = SWRAD1 + RATIOF * (SWRAD2 - SWRAD1)
    elif IDIAGN == int(ONE):
        # DAILY
        # SWRAD = SWRAD1 + RATIOD * (SWRAD2 - SWRAD1)
        # MONTHLY
        SWRAD = SWRAD1 + RATIOF * (SWRAD2 - SWRAD1)

    # INTERPOLATE T&S PROFILES
    TSTAR = np.zeros((array_length,KB))
    SSTAR = np.zeros((array_length,KB))
    TSTAR[:] = TCLIM1[:] + RATIOF * (TCLIM2[:] - TCLIM1[:])
    SSTAR[:] = SCLIM1[:] + RATIOF * (SCLIM2[:] - SCLIM1[:])
    WGEN[:]  = WCLIM1[:] + RATIOF * (WCLIM2[:] - WCLIM1[:])
    if RATIOF <= 0.5:
        WEDDY[:] = WEDDY1[:]
    else:
        WEDDY[:] = WEDDY2[:]

    if IDIAGN == int(ZERO):
        TSURF = TSTAR[0]
        SSURF = SSTAR[0]
    elif IDIAGN == int(ONE):
        TF[:] = TSTAR[:]
        SF[:] = SSTAR[:]

    # INTERPOLATE SUSPENDED INORGANIC MATTER
    ISM[:] = ISM1[:] + RATIOF * (ISM2[:] - ISM1[:])

    # INTERPOLATE SURFACE NUTRIENTS
    NO3SURF = NO3_s1 + RATIOF * (NO3_s2 - NO3_s1)
    NH4SURF = NH4_s1 + RATIOF * (NH4_s2 - NH4_s1)
    PO4SURF = PO4_s1 + RATIOF * (PO4_s2 - PO4_s1)
    SIO4SURF = SIO4_s1 + RATIOF * (SIO4_s2 - SIO4_s1)

    # INTERPOLATE BOTTOM NUTRIENTS
    O2BOTT = O2_b1 + RATIOF * (O2_b2 - O2_b1)
    NO3BOTT = NO3_b1 + RATIOF * (NO3_b2 - NO3_b1)
    PO4BOTT = PO4_b1 + RATIOF * (PO4_b2 - PO4_b1)
    PONBOTTgrad = PON_b1 + RATIOF * (PON_b2 - PON_b1)

    if IFINT == IFCHGE:

        # A MONTH HAS GONE...IT IS NECESSARY TO...
        # ....UPDATE MONTH COUNTER....
        ICOUNTF = ICOUNTF + 1
        print('ICOUNTF = ',ICOUNTF)

        # ....RESET INTERPOLATOR....
        IFINT = int(ZERO)

        # ....SHIFT THE MONTHLY DATA....
        WSU1 = WSU2
        WSV1 = WSV2
        SWRAD1 = SWRAD2
        WTSURF1 = WTSURF2
#         SLUX1 = SLUX2
        NO3_s1 = NO3_s2
        NH4_s1 = NH4_s2
        PO4_s1 = PO4_s2
        SIO4_s1 = SIO4_s2
        NO3_b1 = NO3_b2
        O2_b1 = O2_b2
        PO4_b1 = PO4_b2
        PON_b1 = PON_b2
        ISM1[:]     = ISM2[:]
        TCLIM1[:]   = TCLIM2[:]
        SCLIM1[:]   = SCLIM2[:]
        WCLIM1[:]   = WCLIM2[:]
        WEDDY1[:]   = WEDDY3[:]
        WEDDY2[:]   = WEDDY4[:]

        if ICOUNTF > 13:
            # IF 12 MONTHS HAVE GONE, RESTART THE READING SEQUENCE
            ICOUNTF = 2

            WSU1, WSV1, SSS1, SLUX1, ISM1, SCLIM1, TCLIM1, WCLIM1, WEDDY1, WEDDY2, SB, TB, \
                SWRAD1, WTSURF1, QCORR1, NO3_s1, NH4_s1, PO4_s1, SIO4_s1, O2_b1, NO3_b1, PO4_b1, PON_b1 = read_pom_input()

            WSU1 = -WSU1 * float(1.E-03)
            WSV1 = -WSV1 * float(1.E-03)
            SWRAD1 = -SWRAD1 / RCP
            WTSURF1 = -WTSURF1 / RCP

            WSU2, WSV2, SSS2, SLUX2, ISM2, SCLIM2, TCLIM2, WCLIM2, WEDDY3, WEDDY4, SB2, TB2, \
                SWRAD2, WTSURF2, QCORR2, NO3_s2, NH4_s2, PO4_s2, SIO4_s2, O2_b2, NO3_b2, PO4_b2, PON_b2 = read_pom_input()

            WSU2 = -WSU2 * float(1.E-03)
            WSV2 = -WSV2 * float(1.E-03)
            SWRAD2 = -SWRAD2 / RCP
            WTSURF2 = -WTSURF2 / RCP

    return




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

