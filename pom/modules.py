# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
import numpy as np

ONE = 1.
SEC_PER_DAY = 86400.
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
SWRAD = float()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   N.B.
#   THE FOLLOWING SCALARS ARE USED ONLY WHEN THE MODEL IS RUN IN PROGNOSTIC MODE.
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# LOSS TERM OF THE SURFACE HEAT FLUX
WTSURF = float()

# PRESCRIBED SURFACE T & S
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   TO BE USED (IF DESIRED) AS SURFACE
#   T & S SURFACE BOUNDARY CONDITION
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
TSURF = float(); SSURF = float()

# SURFACE SALINITY FLUX
WSSURF = float()

# LATERAL ADVECTION FLUX FOR T & S
WTADV = np.empty(KB,dtype=float); WSADV = np.empty(KB,dtype=float)


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

#   EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

