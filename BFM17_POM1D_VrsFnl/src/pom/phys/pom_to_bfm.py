# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: Profq
#
# DESCRIPTION
#
#   This subroutine passes the physical variables to the BFM
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
import numpy as np
from decimal import *


def pom_to_bfm():

    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleGlobalMem import RLEN
    from BFM17_POM1D_VrsFnl.src.pom.phys.POMModule import KB, RCP, TB, SB, RHO, H, DZ, SWRAD,\
        WUSURF, WVSURF, KH, KM, U, V, Q2, Q2L, RHO, L
    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleService import ISM, WGEN, WEDDY
    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleMem import ETW, ESW, EIR, ESS, ERHO, ewind, depth

    getcontext().prec = 12  # 12-digit precision (ilong)

    tauw = Decimal()

    for K in range(0,KB - 1):
        ETW[K] = TB[K]
        ESW[K] = SB[K]
        ERHO[K] = (RHO[K] * 1.E3) + 1.E3
        ESS[K] = ISM[K]
        depth[K] = DZ[K] * H

    EIR[0] = (-1.0) * SWRAD * RCP

    tauw = np.sqrt(WUSURF ** 2 + WVSURF ** 2) * 1.E03
    ewind = np.sqrt(tauw/(1.25 * 0.0014))

    return

# EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
