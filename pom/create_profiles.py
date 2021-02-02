from cppdefs import *
from pom.modules import *

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: PROFE
#
# DESCRIPTION:  This subroutine solves for vertical diffusivity.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def PROFE(FF, WFSURF, FSURF, NBC, DT2):

    UMOLPR = 1.E-05
    DH = H

    for K in range(1, KB - 1):
        A[K - 1] = -DT2 * (KH[K] + UMOLPR) / (DZ[K - 1] * DZZ[K - 1] * DH * DH)
        C[K] = -DT2 * (KH[K] + UMOLPR) / (DZ[K] * DZZ[K - 1] * DH * DH)

    VH[0] = A[0] / (A[0]-1.)
    VHP[0] = -DT2 * WFSURF / (-DZ[0]*DH) - FF[0]
    VHP[0] = VHP[0] / (A[0]-1.)

    VH[0] = 0.
    VHP[0] = FSURF

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SECTION SOLVES THE EQUATION
    #   DT2*(KH*FF')' -FF = FB
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for K in range(1, KB - 2):
        VHP[K] = 1. / (A[K] + C[K] * (1. - VH[K - 1]) - 1.)
        VH[K] = A[K] * VHP[K]
        VHP[K] = (C[K] * VHP[K - 1] - FF[K]) * VHP[K]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INSTEAD OF MATCHING SOLUTION TO A LOWER LAYER VALUE OF F(KB-1),
    #   ONE MAY IMPOSE AN ADIABATIC BOTTOM B.C. BY REMOVING C'S FROM
    #   C OL. 1 OF THE NEXT TWO LINES.
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    FF[KB - 2] = (C[KB - 2] * VHP[KB - 3] - FF[KB - 2]) / (C[KB - 2] * (1. - VH[KB - 1]) - 1.)

    for K in range(1, KB - 1):
        KI = KB - K
        FF[KI] = VH[KI] * FF[KI + 1] + VHP[KI]

    for K in range(0, KB):
        VH[K] = 0.0
        VHP[K] = 0.0
        A[K] = 0.0
        C[K] = 0.0

    return FF


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: PROFQ1D
#
# DESCRIPTION:  This subroutine solves for the turbulent closure.
#               Turbulent kinetic energy (Q2/2)
#               Turbulent length scale (Q2l)
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def PROFQ(DT2):


    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL ARRAYS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    BOYGR = np.empty(KB,dtype=float); CC = np.empty(KB,dtype=float)
    TEMP1 = np.empty(KB,dtype=float); TEMP2 = np.empty(KB,dtype=float); TEMP3 = np.empty(KB,dtype=float)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   DATA STATEMENTS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    A1 = 0.92
    B1 = 16.6
    A2 = 0.74
    B2 = 10.1
    C1 = 0.08
    E1 = 1.8
    E2 = 1.33
    E3 = 1.0
    KAPPA = 0.40
    SQ = 0.2
    CIWC = 1.0
    GEE = 9.806
    # SM = KB*0.39
    # SH = KB*0.49
    # GM = KB*0.154
    # GH = KB*0.154

    DH = H
    for K in range(1, KB - 1):
        A[K] = -DT2 * (KQ[K + 1] + KQ[K] + 2 * UMOL) * 0.5 / (DZZ[K - 1] * DZ[K] * DH * DH)
        C[K] = -DT2 * (KQ[K - 1] + KQ[K] + 2 * UMOL) * 0.5 / (DZZ[K - 1] * DZ[K - 1] * DH * DH)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SECTION SOLVES FOR THE EQUATION
    #   DT2*(KQ*Q2')' - Q2*(2.*DT2*DTEF+1.) = -Q2B
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    CONST1 = 16.6 ** 0.6666667 * CIWC

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   BOUNDARY CONDITIONS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    VH[0] = 0.0
    VHP[0] = np.sqrt(WUSURF ** 2 + WVSURF ** 2) * CONST1
    Q2F[KB - 1] = 0.5 * np.sqrt((WUBOT + WUBOT) ** 2 + (WVBOT + WVBOT) ** 2) * CONST1

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   CALCULATE PRESSURE IN UNITS OF DECIBARS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    for K in range(0, KB - 1):
        P = -GEE * 1.025 * ZZ[K] * DH * .1
        CC[K] = 1449.2 + 1.34 * (S[K]-35.) + 4.55*T[K] - 0.045*T[K]**2 + 0.00821*P + (15.0**1.e-9*P**2)
        TEMP1[K] = 2./CC[K]
        TEMP2[K] = (0.00821*P)
        TEMP3[K] = (1.-0.40 * (P/CC[K]**2))

    for K in range(0, KB - 1):
        P = -GEE * 1.025 * ZZ[K] * DH * .1
        CC[K] = CC[K] * (1. - TEMP1[K] * (TEMP2[K] + 15. * 1.e-9 * P ** 2) * TEMP3[K]) ** (-0.5)
        CC[K] = 1449.1 + .00821 * P + 4.55 * T[K] - .045 * T[K] ** 2 + 1.34 * (S[K] - 35.)
        CC[K] = CC[K] / np.sqrt((1. - .01642 * P / CC[K]) * (1. - 0.40 * P / CC[K] ** 2))

    for K in range(1, KB - 1):
        Q2B[K] = np.abs(Q2B[K])
        Q2LB[K] = np.abs(Q2LB[K])
        BOYGR[K] = GEE * (RHO[K - 1] - RHO[K]) / (
                    DZZ[K - 1] * DH)  # & (G) +GEE ** 2 * 2. * 1.025 / (CC(K - 1) ** 2 + CC(K) ** 2)(G)
        DTEF[K] = Q2B[K] * np.sqrt(Q2B[K]) / (B1 * Q2LB[K] + SMALL)
        SPROD[K] = .25 * KM[K] * ((U[K] + U[K] - U[K - 1] - U[K - 1]) ** 2 + (V[K] + V[K] - V[K - 1] - V[K - 1]) ** 2) / (DZZ[K - 1] * DH) ** 2 * CIWC ** 2
        BPROD[K] = KH[K] * BOYGR[K]
        PROD[K] = SPROD[K] + BPROD[K]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP DOWNWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for K in range(1, KB - 1):
        VHP[K] = 1. / (A[K] + C[K] * (1. - VH[K - 1]) - (2. * DT2 * DTEF[K] + 1.))
        VH[K] = A[K] * VHP[K]
        VHP[K] = (-2. * DT2 * PROD[K] + C[K] * VHP[K - 1] - Q2B[K]) * VHP[K]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP UPWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for K in range(0, KB - 1):  # 104
        KI = KB - K
        Q2F[KI] = VH[KI] * Q2F[KI + 1] + VHP[KI]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SEECTION SOLVES FOR TEH EQUATION
    #   DT2(KQ*Q2L')' - Q2L*(DT2*DTEF+1.) = -Q2LB
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   BOUNDARY CONDITIONS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    VH[0] = 0.
    VHP[0] = 0.
    Q2LF[KB - 1] = 0.

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP DOWNWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for K in range(1, KB - 1):
        DTEF[K] = DTEF[K] * (
                    1. + E2 * ((1. / np.abs(Z[K] - Z[0]) + 1. / np.abs(Z[K] - Z[KB])) * L[K] / (DH * KAPPA)) ** 2)
        VHP[K] = 1. / (A[K] + C[K] * (1. - VH[K - 1]) - (DT2 * DTEF[K] + 1.))
        VH[K] = A[K] * VHP[K]
        VHP[K] = (DT2 * (- (SPROD[K] + E3 * BPROD[K]) * L[K] * E1) + C[K] * VHP[K - 1] - Q2LB[K]) * VHP[K]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP UPWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for K in range(0, KB - 1):
        KI = KB - K
        Q2LF[KI] = VH[KI] * Q2LF[KI + 1] + VHP[KI]

    for K in range(1, KB - 1):
        if Q2F[K] > SMALL or Q2LF[K] > SMALL:
            break
        Q2F[K] = SMALL
        Q2LF[K] = SMALL

    for K in range(0, KB - 1):
        Q2F[K] = np.abs(Q2F[K])
        Q2LF[K] = np.abs(Q2LF[K])

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SECTION SOLVES FOR KM AND KH
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    COEF1 = A2 * (1. - 6. * A1 / B1)
    COEF2 = 3. * A2 * B2 + 18. * A1 * A2
    COEF3 = A1 * (1. - 3. * C1 - 6. * A1 / B1)
    COEF4 = 18. * A1 * A1 + 9. * A1 * A2
    COEF5 = 9. * A1 * A2

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   NOTE THAT SM AND SH LIMIT TO INFINITY WHEN GH APPROACHES 0.0288
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    L[0] = 0.
    L[KB - 1] = 0.
    GH[0] = 0.
    GH[KB - 1] = 0.

    for K in range(1, KB - 1):
        L[K] = Q2LF[K] / Q2F[K]
        GH[K] = L[K] ** 2 / Q2F[K] * BOYGR[K]

    for K in range(0, KB):
        GH[K] = np.mininimum(GH[K], .028)
        SH[K] = COEF1 / (1. - COEF2 * GH[K])
        SM[K] = COEF3 + SH(K) * COEF4 * GH[K]
        SM[K] = SM[K] / (1. - COEF5 * GH[K])

    for K in range(0, KB):
        KN[K] = L[K] * np.sqrt(np.abs(Q2[K]))
        KQ[K] = (KN[K] * .41 * SM[K] + KQ[K]) * .5
        #   KQ[K]= (KN[K] * .41 * SH[K] + KQ[K]) * .5
        KM[K] = (KN[K] * SM[K] + KM[K]) * .5
        KH[K] = (KN[K] * SH[K] + KH[K]) * .5

    return Q2F, Q2LF, KM, KH


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: PROFTS
#
# DESCRIPTION:  This subroutine solves for the conservative (Temperature and Salinity)
#               and non-conservative (BFM state var's) scalars of BFM-POM1D.
#               It handles the surface and bottom boundary condition.
#               When used to compute temperature it handles also the solar radiation
#               penetration along the water column, Based on:
#
#               Paulson C. A., Simpson J.J. (1977)
#               Irradiance measurements in the upper ocean.
#               Journal of Physical Oceanography, 7, 952-956.
#
#               Note that when the system is run in diagnostic mode (prescribed
#               Temperature and salinity values), the soutine is used only to compute
#               the vertical profiles of the non-conservative BFM scalars.
#
# THE ROUTINE DUMMY ARGUMENTS ARE:
#               FF:     Property to be computed
#               WFSURF: Property surface flux (for temperature it lacks the incoming surface solar radiation).
#               WFBOT:  Property bottom flux.
#               SWRAD:  Incoming solar radiation
#               FSURF:  Prescribed surface property value
#               NBC:    Flag for definition of the surface boundary condition
#               DT2:    Twice the Time step.
#               NTP:    Flag to choose the Optical (Jerlov) Water type
#               UMOL:   Background diffusivity.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def PROFTS(FF, WFSURF, WFBOT, SWRAD, FSURF, NBC, DT2, NTP, UMOL):

    # FLAG FOR BOUNDARY CONDITION DEFINITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   NBC=1: SURF. B.C. IS WFSURF+SWRAD. NO RADIATIVE PENETRATION.
    #   NBC=2; SURF. B.C. IS WFSURF. SWRAD PENETRATES WATER COLUMN.
    #   NBC=3; SURF. B.C. IS TSURF. NO SWRAD RADIATIVE PENETRATION.
    #   NBC=4; SURF. B.C. IS TSURF. SWRAD PENETRATES WATER COLUMN.
    #
    #   NOTE THAT WTSURF (=WFSURF) AND SWRAD ARE NEGATIVE VALUES WHEN FLUX IS "IN" THE WATER COLUMN
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # FLAG FOR JERLOV WATER TYPE CHOICE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   JERLOV WATER TYPE CHOICE IS RELEVANT ONLY WHEN NBC = 2 OR NBC = 4.
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # SW PROFILE
    RAD = np.empty(KB,dtype=float)

    # IRRADIANCE PARAMETERS AFTER PAULSON & SIMPSON JPO 1977, 952-956
    RP = [0.58, 0.62, 0.67, 0.77, 0.78]
    AD1 = [0.35, 0.60, 1.00, 1.50, 1.40]
    AD2 = [23.00, 20.00, 17.00, 14.00, 7.90]

    # JERLOV WATER TYPES
    # NTP         = 1           2            3           4          5
    # JERLOV TYPE = I           IA           IB          II         III

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

    return FF


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: PROFU
#
# DESCRIPTION:  This subroutine solves for the equation
#               DT2*(KM*U')' - U= -UB
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def PROFU(DT2):

    # UMOL1 = 0.0007

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

    # WUBOT = -CBC * UF[KB - 2]  # 92
    for K in range(0, KB):
        VH[K] = 0.
        VHP[K] = 0.
        A[K] = 0.
        C[K] = 0.

    return UF


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: PROFV
#
# DESCRIPTION:  This subroutine solves for the equation
#               DT2*(KM*V')' - V= -VB
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def PROFV(DT2):

    DH = H

    for K in range(1, KB - 1):
        A[K - 1] = -DT2 * (KM[K] + UMOL) / (DZ[K - 1] * DZZ[K - 1] * DH * DH)
        C[K] = -DT2 * (KM[K] + UMOL) / (DZ[K] * DZZ[K - 1] * DH * DH)

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-DT2 * WVSURF / (-DZ[0] * DH) - VF[0]) / (A[0] - 1.)

    # 98 CONTINUE

    for K in range(1, KB - 2):
        VHP[K] = 1. / (A[K] + C[K] * (1. - VH[K - 1]) - 1.)
        VH[K] = A[K] * VHP[K]
        VHP[K] = (C[K] * VHP[K - 1] - VF[K]) * VHP[K]

    # IF(NO_BOT_STRESS)THEN
    CBC = 0.0
    # ELSE
    # 104      CBC = MAX(.0025,.16/LOG((ZZ(KB-1)-Z(KB))*DH/.01)**2)
    #          CBC = CBC*SQRT((.25* (UB(KB-1)+UB(KB-1)+UB(KB-1)+UB(KB-1)))**2 + VB(KB-1)**2)
    # ENDIF
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   TO RESTORE BOTTOM B.L. DELETE NEXT LINE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    VF[KB - 2] = (C[KB - 1] * VHP[KB - 3] - VF[KB - 2]) / (
                CBC * DT2 / (-DZ[KB - 2] * DH) - 1. - (VH[KB - 3] - 1.) * C[KB - 2])

    for K in range(1, KB - 1):
        KI = KB - K
        VF[KI] = VH[KI] * VF[KI + 1] + VHP[KI]

    # WVBOT = -CBC * VF[KB - 2]  # 92
    for K in range(0, KB):
        VH[K] = 0.
        VHP[K] = 0.
        A[K] = 0.
        C[K] = 0.

    return VF


# EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
