# WARNING THIS IS A TEST VERSION
# (MONTHLY FREQUENCY)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  BFM - Biogeochemical Flux Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleGlobalMem import RLEN
from BFM17_POM1D_VrsFnl.src.pom.phys.POMModule import KB
from decimal import *
import numpy as np

from BFM17_POM1D_VrsFnl.src.pom.phys.opendat import opendat

getcontext().prec = 12  # 12-digit precision (ilong)

# MONTHLY SHORTWAVE RADIATION

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   N.B.!
#   ALWAYS NEEDED: WHEN THE MODEL IS RUN IN DIAGNOSTIC MODE PROVIDES ONLY PAR TO BFM.
#                  IN PROGNOSTIC CONTRIBUTES TO THE DEFINITION OF THE TEMPERATURE SURFACE BOUNDARY CONDITION.
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


SWRAD1, SWRAD2 = Decimal()
# SLUX1, SLUX2 = Decimal()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   N.B.!
#   THE FOLLOWING SCALARS ARE USED ONLY WHEN THE MODEL IS RUN IN PROGNOSTIC MODE.
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# MONTHLY LOSS TERM OF THE SURFACE HEAT FLUX
WTSURF1, WTSURF2 = Decimal()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   N.B.!
#   THE FOLLOWING ARE ALWAYS USED.
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# PRESCRIBED T&S PROFILES
TSTAR, SSTAR = np.empty(KB, dtype=float)

# MONTHLY WIND STRESS
WSU1, WSU2, WSV1, WSV2 = Decimal()

# MONTHLY SURFACE SALINITY
SSS1, SSS2 = Decimal()

# MONTHLY BOTTOM OXYGEN
O2_b1, O2_b2 = Decimal()

# MONTHLY SURFACE AND BOTTOM NITRATE
NO3_s1, NO3_s2 = Decimal()
NO3_b1, NO3_b2 = Decimal()

# MONTHLY SURFACE AND BOTTOM PHOSPHATE
PO4_s1, PO4_s2 = Decimal()
PO4_b1, PO4_b2 = Decimal()

# MONTHLY SURFACE AMMONIA
NH4_s1, NH4_s2 = Decimal()

# MONTHLY BOTTOM PON GRADIENT
PON_b1, PON_b2 = Decimal()

# MONTHLY SURFACE SILICATE
SIO4_s1, SIO4_s2 = Decimal()

# MONTHLY PROFILES OF INORGANIC SUSPENDED MATTER
ISM1, ISM2 = np.empty(KB, dtype=float)

# MONTHLY PROFILES OF T & S
TCLIM1, TCLIM2 = np.empty(KB, dtype=float)
SCLIM1, SCLIM2 = np.empty(KB, dtype=float)
WCLIM1, WCLIM2 = np.empty(KB, dtype=float)
WEDDY1, WEDDY2, WEDDY3, WEDDY4 = np.empty(KB, dtype=float)

SLUX1, QCORR1, QCORR2 = Decimal()

# INTERPOLATORS AND COUNTERS
ICOUNTF, IDOUNTF, \
IFCHGE, IFDCHGE, \
IFDINT, IFINT = Decimal()

RATIOF, RATIOD = Decimal()


def FORCING_MANAGER():
    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleGlobalMem import RLEN, ZERO, PI, ONE, NML_OPEN, NML_READ, \
        error_msg_prn
    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleConstants import SEC_PER_DAY
    from BFM17_POM1D_VrsFnl.src.pom.phys.POMModule import IDIAGN, DTI, intt, RCP, KB, TF, SF, WUSURF, WVSURF, SWRAD, \
        WTSURF, WSSURF, TSURF, SSURF, TB, SB
    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleService import ISM, savef, PO4SURF, NO3SURF, NH4SURF, SIO4SURF, \
        PO4BOTT, NO3BOTT, O2BOTT, PONBOTTgrad, WGEN, WEDDY, wind_input, radiance_input, ism_input, Sprofile_input, \
        W_input, \
        Tprofile_input, surfNut_input, bottNut_input, Sal_input, Temp_input, heat_input, surfaceS_input
    import numpy as np

    getcontext().prec = 12  # 12-digit precision (ilong)

    # LOOP COUNTER
    K = Decimal()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL ARRAYS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    RLENGTH = Decimal()

    # INITIALISATION AND FIRST FORCING READING
    if intt == int(ONE):
        opendat(WSU1, WSV1, SSS1, SLUX1, ISM1, SCLIM1, TCLIM1, WCLIM1,
                WEDDY1, WEDDY2, SB, TB, SWRAD1, WTSURF1, QCORR1, NO3_s1, NH4_s1,
                PO4_s1, SIO4_s1, O2_b1, NO3_b1, PO4_b1, PON_b1)

    # DAY COUNTER
    IDOUNTF = 1

    # MONTH COUNTER
    ICOUNTF = 1

    # TIME STEPS TO COVER ONE DAY
    IFDCHGE = int(SEC_PER_DAY)/int(DTI)

    # TIME STEPS TO COVER ONE MONTH
    IFCHGE = Decimal(30)*IFDCHGE

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

    IFINT = (IFCHGE / Decimal(2)) - int(ONE)

    # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
    # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
    # **                                                                   **
    # **  INITIAL READING OF THE MONTHLY FORCING                           **
    # **                                                                   **
    # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
    # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **

    # WIND STRESS
    









