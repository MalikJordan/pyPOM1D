# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ONE-DIMENSIONAL BFM-POM SYSTEM
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # MODULE POM
#
# DESCRIPTION
#
#   Definition and allocation of parameters, scalars, and arrays used (mostly)
#   by the physical component of the system.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from decimal import *
from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleGlobalMem import RLEN, ONE
from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleConstants import SEC_PER_DAY
import numpy as np

# SET INTEGER PRECISION
getcontext().prec = 12

# SWITCH FOR PROGNOSTIC/DIAGNOSTIC SIMULATION
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   IDIAGN=0 PROGNOSTIC (T&S PROFILES COMPUTED)
#   IDIAGN=1 DIAGNOSTIC (T&S PROFILES PRESCRIBED)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
IDIAGN = Decimal()

# # OF SURFACE/BOTTOM LAYERS WITH LOG DISTRIBUTION
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   KL1: SURFACE LAYERS
#   KL2: BOTTOM LAYERS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
KL1, KL2 = Decimal()

# SWITCH FOR COLD/HOT START
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   IHOTST=0 "COLD" START FROM INITIAL CONDITION
#   IHOTST=1 "HOT" START FROM RESTART FILE
#
#   SEE SUBROUTINE "CALCDEPTH"
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
IHOTST = Decimal()

# MODEL TIME STEP
DTI = Decimal()

# LENGTH OF THE RUN (DAYS)
IDAYS = Decimal()

# ITERATIONS NEEDED FOR AN "IDAYS" RUN
IEND = Decimal()

# COUNTER FOR THE TIME MARCHING LOOP
intt = Decimal

# RUNNING TIME
TIME = Decimal()

# TIME AT RESTART
TIME0 = Decimal()

# BOTTOM DEPTH
H = Decimal()

# LATITUDE & LONGITUDE
ALAT, ALON = Decimal()

# CORIOLIS PARAMETER
COR = Decimal()

# BACKGROUND DIFFUSION FOR U, V, Q2, Q2L
UMOL = Decimal()

# BACKGROUND DIFFUSION FOR T,S & BFM TRACERS
UMOLT, UMOLS, UMOLBFM = Decimal()

# NUTRIENT RELAXATION TIME
NRT_o2o = Decimal()
NRT_n1p = Decimal()
NRT_n3n = Decimal()
NRT_n4n = Decimal()

# PARAMETER FOR THE HASSELIN FILTER
# (TIME STEP MIXING)
SMOTH = Decimal()

# SPECIFIC HEAT TIMES RHO0
RCP = 4.187E6

# 1 DAY IN SECONDS (RECIPROCAL)
DAYI = ONE/SEC_PER_DAY

# VERTICAL LAYERS
KB = 151

# FLAGS TO CHOOSE T,S AND BFM TRACERS SURFACE B.C. IN PROFTS
NBCT, NBCS, NBCBFM = Decimal()

# FLAG TO CHOOSE JERLOV WATER TYPE IN PROFTS
NTP = Decimal()

# T&S RELAXATION TIME (DAYS) FOR LATERAL ADVECTION
TRT, SRT = Decimal()

# DEPTH (m) AT WHICH LATERAL ADVECTION STARTS
upperH = Decimal()

# RELAXATION TIME (DAYS) FOR SURFACE SALINITY FLUX
SSRT = Decimal()

# VERTICAL COORDINATE SYSTEM
Z, ZZ, DZ, DZZ, DZR = np.empty(KB, dtype=float)

# TEMPERATURE
TF, T, TB = np.empty(KB, dtype=float)

# SALINITY
SF, S, SB = np.empty(KB, dtype=float)

# DESITY
RHO = np.empty(KB, dtype=float)

# VELOCITY
UF, U, UB, VF, V, VB = np.empty(KB, dtype=float)

# TURBULENT KINETIC ENERGY (T.K.E.X2)
Q2F, Q2, Q2B = np.empty(KB, dtype=float)

# LENGTH SCALE
L = np.empty(KB, dtype=float)

# (T.K.E.X2) TIMES LENGTH SCALE
Q2LF, Q2L, Q2LB = np.empty(KB, dtype=float)

# VERTICAL DIFFUSION COEFFICIENTS
KM, KH, KQ = np.empty(KB, dtype=float)

# SERVICE ARRAYS USED IN POM ROUTINES
GM, GH, SM, SH, KN, SPROD, BPROD, A, C, VH, VHP, \
    PROD, DTEF, D, DT = np.empty(KB, dtype=float)

# WIND STRESS
WUSURF, WVSURF = Decimal()

# BOTTOM STRESS
WUBOT, WVBOT = Decimal()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   N.B.
#   WHEN THE MODEL IS RUN IN DIAGNOSTIC MODE, THIS IS USED ONLY TO PROVIDE PAR TO THE BIOGEOCHEMICAL COMPONENT.
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
SWRAD = Decimal()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   N.B.
#   THE FOLLOWING SCALARS ARE USED ONLY WHEN THE MODEL IS RUN IN PROGNOSTIC MODE.
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# LOSS TERM OF THE SURFACE HEAT FLUX
WTSURF = Decimal()

# PRESCRIBED SURFACE T & S
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   TO BE USED (IF DESIRED) AS SURFACE
#   T & S SURFACE BOUNDARY CONDITION
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
TSURF, SSURF = Decimal()

# SURFACE SALINITY FLUX
WSSURF = Decimal()

# LATERAL ADVECTION FLUX FOR T & S
WTADV, WSADV = np.empty(KB, dtype=float)


# EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
