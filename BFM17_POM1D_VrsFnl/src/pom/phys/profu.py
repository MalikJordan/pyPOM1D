# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: PROFU
#
# DESCRIPTION
#
#   This subroutine solves for the equation
#   DT2*(KM*U')' - U= -UB
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
import numpy as np
from decimal import *


def PROFU(DT2):

    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleGlobalMem import RLEN
    from BFM17_POM1D_VrsFnl.src.pom.phys.POMModule import H, KB, A, C, KM, DZ, DZZ, VH, VHP, \
        WUSURF, UF, UB, VB, WUBOT, UMOL, Z, ZZ

    getcontext().prec = 12  # 12-digit precision (ilong)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SCALAR ARGUMENTS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # DT2 = Decimal()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL SCALARS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    CBC, DH, UMOL1 = Decimal()
    K, KI = Decimal()

    UMOL1 = 0.0007

    DH = H  # 85
    for K in range(1, KB - 1):
        A[K - 1] = -DT2 * (KM[K] + UMOL) / (DZ[K - 1] * DZZ[K - 1] * DH * DH)
        C[K] = -DT2 * (KM[K] + UMOL) / (DZ[K] * DZZ[K - 1] * DH * DH)

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-DT2 * WUSURF / (-DZ[0] * DH) - UF[0]) / (A[0] - 1.)

    for K in range(1, KB - 2):
        VHP[K] = 1. / (A[K] + C[K] * (1. - VH[K - 1]) - 1.)
        VH[K] = A[K] * VHP[K]
        VHP[K] = (C[K] * VHP[K - 1] - UF[K]) * VHP(K)

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-DT2 * WUSURF / (-DZ[0] * DH) - UF[0]) / (A[0] - 1.)

    # IF(NO_BOT_STRESS)THEN
    CBC = 0.0
    # ELSE
    #    CBC = MAX(.0025,.16/LOG((ZZ(KB-1)-Z(KB))*DH/.01)**2)
    #    CBC = CBC*SQRT(UB(KB-1)**2+ (.25* (VB(KB-1)+VB(KB-1)+VB(KB-1)+ &
    #        VB(KB-1)))**2)
    # ENDIF
    #
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    UF[KB - 2] = (C[KB - 2] * VHP[KB - 3] - UF[KB - 2]) / (
                CBC * DT2 / (-DZ[KB - 2] * DH) - 1. - (VH[KB - 3] - 1.) * C[KB - 2])
    for K in range(1, KB - 1):
        KI = KB - K
        UF[KI - 1] = VH[KI - 1] * UF[KI] + VHP[KI - 1]

    WUBOT = -CBC * UF[KB - 2]  # 92
    for K in range(0, KB):
        VH[K] = 0.
        VHP[K] = 0.
        A[K] = 0.
        C[K] = 0.

    return

# EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
