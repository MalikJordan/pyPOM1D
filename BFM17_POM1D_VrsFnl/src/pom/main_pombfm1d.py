# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# This is the "main" program of the coupled numerical model originating by
# the direct on-line coupling of the 1D Version of the Princeton Ocean model
# "POM" and the Biological Flux Model "BFM".
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from decimal import *
import numpy as np
import f90nml

from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleGlobalMem import RLEN, ZERO, PI, ONE
from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleConstants import SEC_PER_DAY
from BFM17_POM1D_VrsFnl.src.pom.coupling.get_TS_IC import get_TS_IC
from BFM17_POM1D_VrsFnl.src.pom.coupling.get_rst import get_rst
from BFM17_POM1D_VrsFnl.src.pom.coupling.pom_ini_bfm_1d import pom_ini_bfm_1D
from BFM17_POM1D_VrsFnl.src.pom.coupling.pom_bfm_1d import pom_bfm_1D
from BFM17_POM1D_VrsFnl.src.pom.coupling.restart_BFM_inPOM import restart_BFM_inPOM
from BFM17_POM1D_VrsFnl.src.pom.coupling.ModuleForcing import FORCING_MANAGER
from BFM17_POM1D_VrsFnl.src.pom.phys.POMModule import IDIAGN, \
    KL1, KL2, \
    IHOTST, \
    DTI, \
    IDAYS, \
    IEND, \
    intt, \
    TIME, \
    TIME0, \
    H, \
    ALAT, \
    COR, \
    UMOL, \
    UMOLT, UMOLS, UMOLBFM, \
    SMOTH, \
    RCP, \
    DAYI, \
    KB, \
    NBCT, NBCS, NBCBFM, \
    NTP, \
    TRT, SRT, \
    upperH, \
    SSRT, \
    Z, ZZ, DZ, DZZ, DZR, \
    TF, T, TB, \
    SF, S, SB, \
    RHO, \
    UF, U, UB, VF, V, VB, \
    Q2F, Q2, Q2B, \
    L, \
    Q2LF, Q2L, Q2LB, \
    KM, KH, KQ, \
    GM, GH, SM, SH, KN, SPROD, BPROD, A, C, VH, VHP, PROD, DTEF, \
    D, DT, \
    WUSURF, WVSURF, \
    WUBOT, WVBOT, \
    SWRAD, \
    WTSURF, \
    TSURF, SSURF, \
    WSSURF, \
    WTADV, WSADV, \
    NRT_o2o, NRT_n1p, NRT_n3n, NRT_n4n
from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleService import savef

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   EXTERNAL SUBROUTINES
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from BFM17_POM1D_VrsFnl.src.pom.phys.dens import DENS
from BFM17_POM1D_VrsFnl.src.pom.phys.profq1d import PROFQ
from BFM17_POM1D_VrsFnl.src.pom.phys.profTS import PROFTS
from BFM17_POM1D_VrsFnl.src.pom.phys.profu import PROFU
from BFM17_POM1D_VrsFnl.src.pom.phys.profv import PROFV
from BFM17_POM1D_VrsFnl.src.pom.phys.calcdepth import CALCDEPTH


getcontext().prec = 12  # 12-digit precision (ilong)

# NAMELIST READING UNIT
namlst = 10

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   LOCAL SCALARS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# LOOP COUNTER
K = int()

# TWICE THE TIME STEP
DT2 = Decimal()

# EARTH ROTATION ANGULAR VELOCITY
OMEGA = Decimal(7.29E-05)

# INTERPOLATED CLIMATOLOGICAL T AND S PROFILES
TSTAR, SSTAR = np.empty(KB, dtype=float)

# FOCING COUNTERS
ICOUNTF, IDOUNTF, IFCHGE, IFDCHGE, \
IFDINT, IFINT, IFNCHGE, IFNINT = Decimal()

# MONTHLY CLIMATOLOGICAL T AND S PROFILES
SCLIM1, SCLIM2, TCLIM1, TCLIM2 = np.empty(KB, dtype=float)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   LOCAL ARRAYS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
WSU1, WSV1, SSS1, SSS2, SLUX1, SLUX2, SWRAD1, WTSURF1, QCORR1, \
    QCORR2, NO3_s1, NO3_s2, NH4_s1, NH4_s2, PO4_s1, PO4_s2, SIO4_s1, SIO4_s2, \
    O2_b1, O2_b2, NO3_b1, NO3_b2, PO4_b1, PO4_b2, WCLIM1, WCLIM2 = Decimal()

params_POMBFM = f90nml.read('params_POMBFM.nml')

# GENERAL INITIALIZATION
Z[:] = ZERO
ZZ[:] = ZERO
DZ[:] = ZERO
DZZ[:] = ZERO
TF[:] = ZERO
T[:] = ZERO
TB[:] = ZERO
SF[:] = ZERO
S[:] = ZERO
SB[:] = ZERO
RHO[:] = ZERO
UF[:] = ZERO
U[:] = ZERO
UB[:] = ZERO
VF[:] = ZERO
V[:] = ZERO
VB[:] = ZERO
Q2F[:] = Decimal(1.E-07)
Q2[:] = Q2F[:]
Q2B[:] = Q2[:]
Q2LF[:] = Decimal(1.E-7)
Q2L[:] = Q2LF[:]
Q2LB[:] = Q2L[:]
L[:] = Decimal(1.0)
L[0] = ZERO
L[KB - 1] = ZERO
KM[:] = ZERO
KH[:] = ZERO
KQ[:] = ZERO
GM[:] = ZERO
GH[:] = ZERO
SM[:] = ZERO
SH[:] = ZERO
KN[:] = ZERO
SPROD[:] = ZERO
BPROD[:] = ZERO
PROD[:] = ZERO
VH[:] = ZERO
VHP[:] = ZERO
DTEF[:] = ZERO
D[:] = ZERO
DT[:] = ZERO
A[:] = ZERO
C[:] = ZERO
WTADV[:] = ZERO
WSADV[:] = ZERO
WUSURF = ZERO
WVSURF = ZERO
WUBOT = ZERO
WVBOT = ZERO
WTSURF = ZERO
SWRAD = ZERO
WSSURF = ZERO
TSURF = ZERO
SSURF = ZERO

# DEFINE VERTICAL COORDINATE
CALCDEPTH(Z, ZZ, DZ, DZZ, KB, KL1, KL2)

DZ[KB - 1] = Decimal(1.E-06)
DZR[:] = ONE / DZ[:]
DZ[KB - 1] = ZERO

# CORIOLIS PARAMETER
COR = Decimal(2.0) * OMEGA * np.sin(ALAT * Decimal(2.0) * PI / Decimal(360))

# TWICE THE TIME STEP
DT2 = Decimal(2.0) * DTI

# ITERATIONS NEEDED TO CARRY OUT AN "IDAYS" SIMULATION
IEND = IDAYS * int(SEC_PER_DAY) / int(DTI)

# READ  T&S INITIAL CONDITIONS (IHOTST=0) OR RESTART FILE (IHOTST=1)
if IHOTST == 0:
    TIME0 = ZERO
    get_TS_IC()
    # DEFINE INITIAL DENSITY FIELD
    DENS(T, S, ZZ, H, RHO, KB)

elif IHOTST == 1:
    # READ RESTART
    get_rst()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   INITIALIZATION OF BFM
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
if POM_only is not None:
    pom_ini_bfm_1D()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   BEGIN THE TIME MARCH
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
print('ICOUNT before time march loop = ',ICOUNTF)

# BFM SET-UP, INITIALISATION AND RESTART READING (IF REQUIRED)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   BEGIN THE TIME MARCH
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
for intt in range(0,IEND):

    # COMPUTE TIME IN DAYS
    TIME = TIME0 + (DTI * float(intt) * DAYI)

    # TURBULENCE CLOSURE
    Q2F[:]  = Q2B[:]
    Q2LF[:] = Q2LB[:]

    PROFQ(DT2)

    # DEFINE ALL FORCINGS
    FORCING_MANAGER()

# T&S COMPUTATION
if IDIAGN == int(ZERO):
    # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
    # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
    # **                                                 **
    # **               PROGNOSTIC MODE:                  **
    # **        T & S FULLY COMPUTED BY MODEL            **
    # **                                                 **
    # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
    # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **

    # COMPUTE LATERAL ADVECTION TERM FOR T&S
    WTADV[:] = ZERO
    WSADV[:] = ZERO

    if TRT != ZERO:
        for K in range(0,KB):
            if -ZZ[K]*H >= upperH:
                WTADV[K] = (TSTAR[K]-T[K])/(TRT*SEC_PER_DAY)

    if SRT != ZERO:
        for K in range(0,KB):
            if -ZZ[K]*H >= upperH:
                WSADV[K] = (SSTAR[K] - S[K]) / (SRT * SEC_PER_DAY)

    # COMPUTE SURFACE SALINITY FLUX
    WSSURF = -(SSURF-S[0])*SSRT/SEC_PER_DAY

    # COMPUTE TEMPERATURE
    TF[:] = TB[:] + (WTADV[:]*DT2)

    PROFTS(TF,WTSURF,0,SWRAD,TSURF,NBCT,DT2,NTP,UMOLT)

    # COMPUTE SALINITY
    SF[:] = SB[:] + (WSADV[:]*DT2)
    PROFTS(SF,WSSURF,0,ZERO,SSURF,NBCS,DT2,NTP,UMOLS)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   MIXING THE TIMESTEP (ASSELIN)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    T[:]  = T[:] + Decimal(0.5) * SMOTH * (TF[:] + TB[:] - Decimal(2.0) * T[:])
    S[:]  = S[:] + Decimal(0.5) * SMOTH * (SF[:] + SB[:] - Decimal(2.0) * S[:])

# COMPUTE VELOCITY
UF[:] = UB[:] + DT2*COR*V[:]
PROFU(DT2)
VF[:] = VB[:] + DT2*COR*U[:]
PROFV(DT2)

# MIX TIME STEP (ASSELIN FILTER)
Q2[:]   = Q2[:]  + Decimal(0.5) * SMOTH * (Q2F[:]  + Q2B[:]  - Decimal(2.0) * Q2[:])
Q2L[:]  = Q2L[:] + Decimal(0.5) * SMOTH * (Q2LF[:] + Q2LB[:] - Decimal(2.0) * Q2L[:])

U[:]  = U[:] + Decimal(0.5) * SMOTH * (UF[:] + UB[:] - Decimal(2.0) * U[:])
V[:]  = V[:] + Decimal(0.5) * SMOTH * (VF[:] + VB[:] - Decimal(2.0) * V[:])

# RESTORE TIME SEQUENCE
Q2B[:]  = Q2[:]
Q2[:]   = Q2F[:]
Q2LB[:] = Q2L[:]
Q2L[:]  = Q2LF[:]

UB[:] = U[:]
U[:]  = UF[:]
VB[:] = V[:]
V[:]  = VF[:]

TB[:] = T[:]
T[:]  = TF[:]
SB[:] = S[:]
S[:]  = SF[:]

# UPDATE DENSITY
DENS(T,S,ZZ,H,RHO,KB)

# PHYSICS HAS BEEN COMPUTED.....GO BFM!!!!

if POM_only is not None:
    pom_bfm_1D()

# WRITING OF POM RESTART
print(TIME,
      U, UB, V, VB,
      T, TB, S, SB,
      Q2, Q2B, Q2L, Q2LB,
      KH, KM, KQ,
      L,
      WUBOT, WVBOT,
      RHO)

# BFM RESTART
restart_BFM_inPOM()

print('MAIN DONE')



















