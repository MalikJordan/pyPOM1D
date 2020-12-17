# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: PROFU
#
# DESCRIPTION
#
# This subroutine solves for the conservative (Temperature and Salinity)
# and non-conservative (BFM state var's) scalars of BFM-POM1D.
# It handles the surface and bottom boundary condition.
# When used to compute temperature it handles also the solar radiation
# penetration along the water column, Based on:
#
# Paulson C. A., Simpson J.J. (1977)
# Irradiance measurements in the upper ocean.
# Journal of Physical Oceanography, 7, 952-956.
#
# Note that when the system is run in diagnostic mode (prescribed
# Temperature and salinity values), the soutine is used only to compute
# the vertical profiles of the non-conservative BFM scalars.
#
# The routine dummy arguments are:
# FF:     Property to be computed
# WFSURF: Property surface flux (for temperature it lacks the incoming surface solar radiation).
# WFBOT:  Property bottom flux.
# SWRAD:  Incoming solar radiation
# FSURF:  Prescribed surface property value
# NBC:    Flag for definition of the surface boundary condition
# DT2:    Twice the Time step.
# NTP:    Flag to choose the Optical (Jerlov) Water type
# UMOL:   Background diffusivity.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
import numpy as np
from decimal import *


def PROFTS(FF, WFSURF, WFBOT, SWRAD, FSURF, NBC, DT2, NTP, UMOL):

    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleGlobalMem import RLEN, ZERO, ONE
    from BFM17_POM1D_VrsFnl.src.pom.phys.POMModule import H, KB, A, C, KH, DZ, DZZ, VH, VHP, Z

    getcontext().prec = 12  # 12-digit precision (ilong)

    # TWICE THE TIME STEP
    # DT2 = Decimal()

    # SURFACE TEMPERATURE / SALINITY / TRACER
    # FSURF = Decimal()

    # SURFACE INCIDENT SHORT WAVE RADIATION
    # SWRAD = Decimal()

    # SURFACE/BOTTOM HEAT FLUX LOSS TERM OR SALINITY / TRACER FLUX
    # WFSURF, WFBOT = Decimal()

    # FLAG FOR BOUNDARY CONDITION DEFINITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   NBC=1: SURF. B.C. IS WFSURF+SWRAD. NO RADIATIVE PENETRATION.
    #   NBC=2; SURF. B.C. IS WFSURF. SWRAD PENETRATES WATER COLUMN.
    #   NBC=3; SURF. B.C. IS TSURF. NO SWRAD RADIATIVE PENETRATION.
    #   NBC=4; SURF. B.C. IS TSURF. SWRAD PENETRATES WATER COLUMN.
    #
    #   NOTE THAT WTSURF (=WFSURF) AND SWRAD ARE NEGATIVE VALUES WHEN FLUX IS "IN" THE WATER COLUMN
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # NBC = Decimal()

    # FLAG FOR JERLOV WATER TYPE CHOICE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   JERLOV WATER TYPE CHOICE IS RELEVANT ONLY WHEN NBC = 2 OR NBC = 4.
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # NTP = Decimal()

    # BACKGROUND DIFFUSIVITY
    # UMOL = Decimal()

    # TEMPERATURE/SALINITY/TRACER
    # FF = np.empty(KB,dtype=float)

    # LOOP COUNTERS
    K, KI = Decimal()

    # SW PROFILE
    RAD = np.empty(KB,dtype=float)

    # IRRADIANCE PARAMETERS AFTER PAULSON & SIMPSON JPO 1977, 952-956
    # RCP, AD1, AD2 = np.empty(5,dtype=float)
    RP, AD1, AD2 = np.empty(5,dtype=float)

    # JERLOV WATER TYPES
    # NTP         = 1           2            3           4          5
    # JERLOV TYPE = I           IA           IB          II         III

    RP = [0.58, 0.62, 0.67, 0.77, 0.78]
    AD1 = [0.35, 0.60, 1.00, 1.50, 1.40]
    AD2 = [23.00, 20.00, 17.00, 14.00, 7.90]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   START COMPUTATION OF VERTICAL PROFILE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    for K in range(1, KB - 1):
        A[K - 1] = -DT2 * (KH[K] + UMOL) / (DZ[K - 1] * DZZ[K - 1] * H * H)
        C[K] = -DT2 * (KH[K] + UMOL) / (DZ[K] * DZZ[K - 1] * H * H)

    RAD[:] = ZERO

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SURFACE BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # *** PENETRATIVE RADIATION CALCULATION. AT THE BOTTOM ANY UNATTENUATED IS DEPOSITED IN THE BOTTOM LAYER.

    if NBC == 1:
        VH[0] = A[0] / (A[0] - ONE)
        VHP[0] = -DT2 * (WFSURF + SWRAD) / (-DZ[0] * H) - FF[0]
        VHP[0] = VHP[0] / (A[0] - ONE)

    elif NBC == 2:
        RAD[:] = SWRAD * (RP[NTP] * np.exp(Z[:] * H / AD1[NTP]) + (ONE - RP[NTP] * np.exp(Z[:] * H / AD2[NTP])))  # ***
        RAD[KB-1] = ZERO

        VH[0] = A[0] / (A[0] - ONE)
        VHP[0] = DT2 * (WFSURF + RAD[0] - RAD[1]) / (DZ[0] * H) - FF[0]
        VHP[0] = VHP[0] / (A[0] - ONE)

    elif NBC == 3:
        VH[0] = ZERO
        VHP[0] = FSURF

    elif NBC == 4:
        RAD[:] = SWRAD * (RP[NTP] * np.exp(Z[:] * H / AD1[NTP]) + (ONE - RP[NTP] * np.exp(Z[:] * H / AD2[NTP])))  # ***
        RAD[KB-1] = 0

        VH[0] = ZERO
        VHP[0] = FSURF

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SECTION SOLVES THE EQUATION
    #   DT2*(KH*FF')' -FF = -FB
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for K in range(1, KB - 2):
        VHP[K] = 1 / (A[K] + C[K] * (1 - VH[K - 1]) - 1)
        VH[K] = A[K] * VHP[K]
        VHP[K] = (C[K] * VHP[K - 1] - FF[K] + DT2 * (RAD[K] - RAD[K + 1]) / (H * DZ[K])) * VHP[K]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   APPLY A NON ADIABATIC BOTTOM BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    FF[KB - 2] = (C[KB - 2] * VHP[KB - 3] - FF[KB - 2] + (WFBOT * DT2 / (DZ[KB - 2] * H))) / (
                C[KB - 2] * (1 - VH[KB - 3]) - 1)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   APPLY A NON ADIABATIC BOTTOM BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for K in range(1, KB - 1):
        KI = KB - K
        FF[KI] = VH[KI] * FF[KI + 1] + VHP[KI]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   ZEROING
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    VH[:] = ZERO
    VHP[:] = ZERO
    A[:] = ZERO
    C[:] = ZERO

    return

# EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
