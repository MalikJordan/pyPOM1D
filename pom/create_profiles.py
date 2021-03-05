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
    DH = bottom_depth



    for K in range(1, vertical_layers - 1):
        A[K - 1] = -DT2 * (diffusion_coefficient_tracers[K] + UMOLPR) / (vertical_spacing[K - 1] * vertical_spacing_staggered[K - 1] * DH * DH)
        C[K] = -DT2 * (diffusion_coefficient_tracers[K] + UMOLPR) / (vertical_spacing[K] * vertical_spacing_staggered[K - 1] * DH * DH)

    VH[0] = A[0] / (A[0]-1.)
    VHP[0] = -DT2 * WFSURF / (-vertical_spacing[0] * DH) - FF[0]
    VHP[0] = VHP[0] / (A[0]-1.)

    VH[0] = 0.
    VHP[0] = FSURF

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SECTION SOLVES THE EQUATION
    #   DT2*(KH*FF')' -FF = FB
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for K in range(1, vertical_layers - 2):
        VHP[K] = 1. / (A[K] + C[K] * (1. - VH[K - 1]) - 1.)
        VH[K] = A[K] * VHP[K]
        VHP[K] = (C[K] * VHP[K - 1] - FF[K]) * VHP[K]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INSTEAD OF MATCHING SOLUTION TO A LOWER LAYER VALUE OF F(KB-1),
    #   ONE MAY IMPOSE AN ADIABATIC BOTTOM B.C. BY REMOVING C'S FROM
    #   C OL. 1 OF THE NEXT TWO LINES.
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    FF[vertical_layers - 2] = (C[vertical_layers - 2] * VHP[vertical_layers - 3] - FF[vertical_layers - 2]) / (C[vertical_layers - 2] * (1. - VH[vertical_layers - 1]) - 1.)

    for K in range(1, vertical_layers - 1):
        KI = vertical_layers - K
        FF[KI] = VH[KI] * FF[KI + 1] + VHP[KI]

    for K in range(0, vertical_layers):
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
def PROFQ(twice_the_timestep):


    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL ARRAYS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    BOYGR = np.empty(vertical_layers, dtype=float); CC = np.empty(vertical_layers, dtype=float)
    TEMP1 = np.empty(vertical_layers, dtype=float); TEMP2 = np.empty(vertical_layers, dtype=float); TEMP3 = np.empty(vertical_layers, dtype=float)

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
    von_karman_constant = 0.40
    SQ = 0.2
    CIWC = 1.0
    gravity = 9.806
    # SM = KB*0.39
    # SH = KB*0.49
    # GM = KB*0.154
    # GH = KB*0.154

    DH = bottom_depth
    for K in range(1, vertical_layers - 1):
        A[K] = -twice_the_timestep * (diffusion_coefficient_kinetic_energy[K + 1] +
                                      diffusion_coefficient_kinetic_energy[K] + 2 * background_diffusion_momentum) * \
               0.5 / (vertical_spacing_staggered[K - 1] * vertical_spacing[K] * DH * DH)
        C[K] = -twice_the_timestep * (diffusion_coefficient_kinetic_energy[K - 1] +
                                      diffusion_coefficient_kinetic_energy[K] + 2 * background_diffusion_momentum) * \
               0.5 / (vertical_spacing_staggered[K - 1] * vertical_spacing[K - 1] * DH * DH)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SECTION SOLVES FOR THE EQUATION
    #   DT2*(KQ*Q2')' - Q2*(2.*DT2*DTEF+1.) = -Q2B
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    CONST1 = 16.6 ** 0.6666667 * CIWC

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   BOUNDARY CONDITIONS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    VH[0] = 0.0
    VHP[0] = np.sqrt(wind_stress_zonal ** 2 + wind_stress_meridional ** 2) * CONST1
    kinetic_energy_forward_time_level[vertical_layers - 1] = \
        0.5 * np.sqrt((bottom_stress_zonal + bottom_stress_zonal) ** 2 +
                      (bottom_stress_meridional + bottom_stress_meridional) ** 2) * CONST1

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   CALCULATE PRESSURE IN UNITS OF DECIBARS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    for K in range(0, vertical_layers - 1):
        pressure = -gravity * 1.025 * vertical_coordinates_staggered[K] * DH * .1
        CC[K] = 1449.2 + 1.34 * (salinity[K] - 35.) + 4.55 * temperature[K] - 0.045 * temperature[K] ** 2 + 0.00821 * pressure + (15.0 ** 1.e-9 * pressure ** 2)
        TEMP1[K] = 2./CC[K]
        TEMP2[K] = (0.00821*pressure)
        TEMP3[K] = (1.-0.40 * (pressure/CC[K]**2))

    for K in range(0, vertical_layers - 1):
        pressure = -gravity * 1.025 * vertical_coordinates_staggered[K] * DH * .1
        CC[K] = CC[K] * (1. - TEMP1[K] * (TEMP2[K] + 15. * 1.e-9 * pressure ** 2) * TEMP3[K]) ** (-0.5)
        CC[K] = 1449.1 + .00821 * pressure + 4.55 * temperature[K] - .045 * temperature[K] ** 2 + 1.34 * (salinity[K] - 35.)
        CC[K] = CC[K] / np.sqrt((1. - .01642 * pressure / CC[K]) * (1. - 0.40 * pressure / CC[K] ** 2))

    for K in range(1, vertical_layers - 1):
        kinetic_energy_backward_time_level[K] = np.abs(kinetic_energy_backward_time_level[K])
        kinetic_energy_times_length_backward_time_level[K] = np.abs(kinetic_energy_times_length_backward_time_level[K])
        BOYGR[K] = gravity * (density[K - 1] - density[K]) / (vertical_spacing_staggered[K - 1] * DH)  # & (G) +GEE ** 2 * 2. * 1.025 / (CC(K - 1) ** 2 + CC(K) ** 2)(G)
        DTEF[K] = kinetic_energy_backward_time_level[K] * np.sqrt(kinetic_energy_backward_time_level[K]) / (B1 * kinetic_energy_times_length_backward_time_level[K] + SMALL)
        SPROD[K] = .25 * diffusion_coefficient_momentum[K] * ((velocity_zonal[K] + velocity_zonal[K] - velocity_zonal[K - 1] - velocity_zonal[K - 1]) ** 2 + (velocity_meridional[K] + velocity_meridional[K] - velocity_meridional[K - 1] - velocity_meridional[K - 1]) ** 2) / (vertical_spacing_staggered[K - 1] * DH) ** 2 * CIWC ** 2
        BPROD[K] = diffusion_coefficient_tracers[K] * BOYGR[K]
        PROD[K] = SPROD[K] + BPROD[K]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP DOWNWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for K in range(1, vertical_layers - 1):
        VHP[K] = 1. / (A[K] + C[K] * (1. - VH[K - 1]) - (2. * twice_the_timestep * DTEF[K] + 1.))
        VH[K] = A[K] * VHP[K]
        VHP[K] = (-2. * twice_the_timestep * PROD[K] + C[K] * VHP[K - 1] - kinetic_energy_backward_time_level[K]) * VHP[K]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP UPWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for K in range(0, vertical_layers - 1):  # 104
        KI = vertical_layers - K
        kinetic_energy_forward_time_level[KI] = VH[KI] * kinetic_energy_forward_time_level[KI + 1] + VHP[KI]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SEECTION SOLVES FOR TEH EQUATION
    #   DT2(KQ*Q2L')' - Q2L*(DT2*DTEF+1.) = -Q2LB
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   BOUNDARY CONDITIONS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    VH[0] = 0.
    VHP[0] = 0.
    kinetic_energy_times_length_forward_time_level[vertical_layers - 1] = 0.

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP DOWNWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for K in range(1, vertical_layers - 1):
        DTEF[K] = DTEF[K] * (1. + E2 * ((1. / np.abs(Z[K] - Z[0]) + 1. / np.abs(Z[K] - Z[vertical_layers])) * L[K] / (DH * von_karman_constant)) ** 2)
        VHP[K] = 1. / (A[K] + C[K] * (1. - VH[K - 1]) - (twice_the_timestep * DTEF[K] + 1.))
        VH[K] = A[K] * VHP[K]
        VHP[K] = (twice_the_timestep * (- (SPROD[K] + E3 * BPROD[K]) * L[K] * E1) + C[K] * VHP[K - 1] - kinetic_energy_times_length_backward_time_level[K]) * VHP[K]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP UPWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for K in range(0, vertical_layers - 1):
        KI = vertical_layers - K
        kinetic_energy_times_length_forward_time_level[KI] = VH[KI] * kinetic_energy_times_length_forward_time_level[KI + 1] + VHP[KI]

    for K in range(1, vertical_layers - 1):
        if kinetic_energy_forward_time_level[K] > SMALL or kinetic_energy_times_length_forward_time_level[K] > SMALL:
            break
        kinetic_energy_forward_time_level[K] = SMALL
        kinetic_energy_times_length_forward_time_level[K] = SMALL

    for K in range(0, vertical_layers - 1):
        kinetic_energy_forward_time_level[K] = np.abs(kinetic_energy_forward_time_level[K])
        kinetic_energy_times_length_forward_time_level[K] = np.abs(kinetic_energy_times_length_forward_time_level[K])

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
    L[vertical_layers - 1] = 0.
    richardson_number[0] = 0.
    richardson_number[vertical_layers - 1] = 0.

    for K in range(1, vertical_layers - 1):
        L[K] = kinetic_energy_times_length_forward_time_level[K] / kinetic_energy_forward_time_level[K]
        richardson_number[K] = L[K] ** 2 / kinetic_energy_forward_time_level[K] * BOYGR[K]

    for K in range(0, vertical_layers):
        richardson_number[K] = np.mininimum(richardson_number[K], .028)
        SH[K] = COEF1 / (1. - COEF2 * richardson_number[K])
        SM[K] = COEF3 + SH(K) * COEF4 * richardson_number[K]
        SM[K] = SM[K] / (1. - COEF5 * richardson_number[K])

    for K in range(0, vertical_layers):
        KN[K] = L[K] * np.sqrt(np.abs(kinetic_energy_current_time_level[K]))
        diffusion_coefficient_kinetic_energy[K] = (KN[K] * .41 * SM[K] + diffusion_coefficient_kinetic_energy[K]) * .5
        #   KQ[K]= (KN[K] * .41 * SH[K] + KQ[K]) * .5
        diffusion_coefficient_momentum[K] = (KN[K] * SM[K] + diffusion_coefficient_momentum[K]) * .5
        diffusion_coefficient_tracers[K] = (KN[K] * SH[K] + diffusion_coefficient_tracers[K]) * .5

    return kinetic_energy_forward_time_level, kinetic_energy_times_length_forward_time_level, diffusion_coefficient_momentum, diffusion_coefficient_tracers, diffusion_coefficient_kinetic_energy


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
    RAD = np.empty(vertical_layers, dtype=float)

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
    for K in range(1, vertical_layers - 1):
        A[K - 1] = -DT2 * (diffusion_coefficient_tracers[K] + UMOL) / (vertical_spacing[K - 1] * vertical_spacing_staggered[K - 1] * bottom_depth * bottom_depth)
        C[K] = -DT2 * (diffusion_coefficient_tracers[K] + UMOL) / (vertical_spacing[K] * vertical_spacing_staggered[K - 1] * bottom_depth * bottom_depth)

    RAD[:] = ZERO

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SURFACE BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # *** PENETRATIVE RADIATION CALCULATION. AT THE BOTTOM ANY UNATTENUATED IS DEPOSITED IN THE BOTTOM LAYER.

    if NBC == 1:
        VH[0] = A[0] / (A[0] - ONE)
        VHP[0] = -DT2 * (WFSURF + SWRAD) / (-vertical_spacing[0] * bottom_depth) - FF[0]
        VHP[0] = VHP[0] / (A[0] - ONE)

    elif NBC == 2:
        RAD[:] = SWRAD * (RP[NTP] * np.exp(Z[:] * bottom_depth / AD1[NTP]) + (ONE - RP[NTP] * np.exp(Z[:] * bottom_depth / AD2[NTP])))  # ***
        RAD[vertical_layers - 1] = ZERO

        VH[0] = A[0] / (A[0] - ONE)
        VHP[0] = DT2 * (WFSURF + RAD[0] - RAD[1]) / (vertical_spacing[0] * bottom_depth) - FF[0]
        VHP[0] = VHP[0] / (A[0] - ONE)

    elif NBC == 3:
        VH[0] = ZERO
        VHP[0] = FSURF

    elif NBC == 4:
        RAD[:] = SWRAD * (RP[NTP] * np.exp(Z[:] * bottom_depth / AD1[NTP]) + (ONE - RP[NTP] * np.exp(Z[:] * bottom_depth / AD2[NTP])))  # ***
        RAD[vertical_layers - 1] = 0

        VH[0] = ZERO
        VHP[0] = FSURF

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SECTION SOLVES THE EQUATION
    #   DT2*(KH*FF')' -FF = -FB
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for K in range(1, vertical_layers - 2):
        VHP[K] = 1 / (A[K] + C[K] * (1 - VH[K - 1]) - 1)
        VH[K] = A[K] * VHP[K]
        VHP[K] = (C[K] * VHP[K - 1] - FF[K] + DT2 * (RAD[K] - RAD[K + 1]) / (bottom_depth * vertical_spacing[K])) * VHP[K]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   APPLY A NON ADIABATIC BOTTOM BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    FF[vertical_layers - 2] = (C[vertical_layers - 2] * VHP[vertical_layers - 3] - FF[vertical_layers - 2] + (WFBOT * DT2 / (vertical_spacing[vertical_layers - 2] * bottom_depth))) / (
            C[vertical_layers - 2] * (1 - VH[vertical_layers - 3]) - 1)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   APPLY A NON ADIABATIC BOTTOM BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for K in range(1, vertical_layers - 1):
        KI = vertical_layers - K
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

    DH = bottom_depth  # 85
    for K in range(1, vertical_layers - 1):
        A[K - 1] = -DT2 * (diffusion_coefficient_momentum[K] + background_diffusion_momentum) / (vertical_spacing[K - 1] * vertical_spacing_staggered[K - 1] * DH * DH)
        C[K] = -DT2 * (diffusion_coefficient_momentum[K] + background_diffusion_momentum) / (vertical_spacing[K] * vertical_spacing_staggered[K - 1] * DH * DH)

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-DT2 * wind_stress_zonal / (-vertical_spacing[0] * DH) - UF[0]) / (A[0] - 1.)

    for K in range(1, vertical_layers - 2):
        VHP[K] = 1. / (A[K] + C[K] * (1. - VH[K - 1]) - 1.)
        VH[K] = A[K] * VHP[K]
        VHP[K] = (C[K] * VHP[K - 1] - UF[K]) * VHP(K)

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-DT2 * wind_stress_zonal / (-vertical_spacing[0] * DH) - UF[0]) / (A[0] - 1.)

    # IF(NO_BOT_STRESS)THEN
    CBC = 0.0
    # ELSE
    #    CBC = MAX(.0025,.16/LOG((ZZ(KB-1)-Z(KB))*DH/.01)**2)
    #    CBC = CBC*SQRT(UB(KB-1)**2+ (.25* (VB(KB-1)+VB(KB-1)+VB(KB-1)+ &
    #        VB(KB-1)))**2)
    # ENDIF
    #
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    UF[vertical_layers - 2] = (C[vertical_layers - 2] * VHP[vertical_layers - 3] - UF[vertical_layers - 2]) / (
            CBC * DT2 / (-vertical_spacing[vertical_layers - 2] * DH) - 1. - (VH[vertical_layers - 3] - 1.) * C[vertical_layers - 2])
    for K in range(1, vertical_layers - 1):
        KI = vertical_layers - K
        UF[KI - 1] = VH[KI - 1] * UF[KI] + VHP[KI - 1]

    # WUBOT = -CBC * UF[KB - 2]  # 92
    for K in range(0, vertical_layers):
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

    DH = bottom_depth

    for K in range(1, vertical_layers - 1):
        A[K - 1] = -DT2 * (diffusion_coefficient_momentum[K] + background_diffusion_momentum) / (vertical_spacing[K - 1] * vertical_spacing_staggered[K - 1] * DH * DH)
        C[K] = -DT2 * (diffusion_coefficient_momentum[K] + background_diffusion_momentum) / (vertical_spacing[K] * vertical_spacing_staggered[K - 1] * DH * DH)

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-DT2 * wind_stress_meridional / (-vertical_spacing[0] * DH) - VF[0]) / (A[0] - 1.)

    # 98 CONTINUE

    for K in range(1, vertical_layers - 2):
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

    VF[vertical_layers - 2] = (C[vertical_layers - 1] * VHP[vertical_layers - 3] - VF[vertical_layers - 2]) / (
            CBC * DT2 / (-vertical_spacing[vertical_layers - 2] * DH) - 1. - (VH[vertical_layers - 3] - 1.) * C[vertical_layers - 2])

    for K in range(1, vertical_layers - 1):
        KI = vertical_layers - K
        VF[KI] = VH[KI] * VF[KI + 1] + VHP[KI]

    # WVBOT = -CBC * VF[KB - 2]  # 92
    for K in range(0, vertical_layers):
        VH[K] = 0.
        VHP[K] = 0.
        A[K] = 0.
        C[K] = 0.

    return VF


# EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
