from cppdefs import *
from pom.modules import vertical_layers, zero, one
from inputs.params_POMBFM import *
import numpy as np

twice_the_timestep = 2. * dti

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: PROFE
#
# DESCRIPTION:  This subroutine solves for vertical diffusivity.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def create_vertical_diffusivity_profile(vertical_spacing, vertical_spacing_staggered, WFSURF, FSURF, nbc):

    A = np.zeros(vertical_layers, dtype=float)
    C = np.zeros(vertical_layers, dtype=float)
    VH = np.zeros(vertical_layers, dtype=float)
    VHP = np.zeros(vertical_layers, dtype=float)
    diffusion_coefficient_tracers = np.zeros(vertical_layers, dtype=float)                      # KH

    vertical_diffusivity = np.zeros(vertical_layers, dtype=float)                               # FF

    umolpr = 1.E-05

    for i in range(1, vertical_layers - 1):
        A[i - 1] = -twice_the_timestep * (diffusion_coefficient_tracers[i] + umolpr) / (vertical_spacing[i - 1] * vertical_spacing_staggered[i - 1] * h * h)
        C[i] = -twice_the_timestep * (diffusion_coefficient_tracers[i] + umolpr) / (vertical_spacing[i] * vertical_spacing_staggered[i - 1] * h * h)

    # need condition for nb for lines 32-37
    VH[0] = A[0] / (A[0]-1.)
    VHP[0] = -twice_the_timestep * WFSURF / (-vertical_spacing[0] * h) - vertical_diffusivity[0]
    VHP[0] = VHP[0] / (A[0]-1.)

    VH[0] = 0.
    VHP[0] = FSURF

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SECTION SOLVES THE EQUATION
    #   DT2*(KH*FF')' -FF = FB
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(1, vertical_layers - 2):
        VHP[i] = 1. / (A[i] + C[i] * (1. - VH[i - 1]) - 1.)
        VH[i] = A[i] * VHP[i]
        VHP[i] = (C[i] * VHP[i - 1] - vertical_diffusivity[i]) * VHP[i]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INSTEAD OF MATCHING SOLUTION TO A LOWER LAYER VALUE OF F(KB-1),
    #   ONE MAY IMPOSE AN ADIABATIC BOTTOM B.C. BY REMOVING C'S FROM
    #   C OL. 1 OF THE NEXT TWO LINES.
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    vertical_diffusivity[vertical_layers - 2] = (C[vertical_layers - 2] * VHP[vertical_layers - 3] - vertical_diffusivity[vertical_layers - 2]) / (C[vertical_layers - 2] * (1. - VH[vertical_layers - 1]) - 1.)

    for i in range(1, vertical_layers - 1):
        j = vertical_layers - i
        vertical_diffusivity[j] = VH[j] * vertical_diffusivity[j + 1] + VHP[j]

    for i in range(0, vertical_layers):
        VH[i] = 0.0
        VHP[i] = 0.0
        A[i] = 0.0
        C[i] = 0.0

    return vertical_diffusivity


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
def PROFQ(vertical_spacing, vertical_spacing_staggered, vertical_coordinates, vertical_coordinates_staggered,
          length_scale, temperature, salinity, density_profile, velocity_zonal, velocity_meridional,
          kinetic_energy_backward, kinetic_energy, kinetic_energy_forward,
          kinetic_energy_times_length_backward, kinetic_energy_times_length, kinetic_energy_times_length_forward,
          wind_stress_zonal, wind_stress_meridional, bottom_stress_zonal, bottom_stress_meridional):

    # INITIALIZE VARIABLES
    diffusion_coefficient_momentum = np.zeros(vertical_layers,dtype=float)
    diffusion_coefficient_tracers = np.zeros(vertical_layers,dtype=float)
    diffusion_coefficient_kinetic_energy = np.zeros(vertical_layers,dtype=float)

    A = np.zeros(vertical_layers, dtype=float)
    C = np.zeros(vertical_layers, dtype=float)
    VH = np.zeros(vertical_layers, dtype=float)
    VHP = np.zeros(vertical_layers, dtype=float)

    KN = np.zeros(vertical_layers, dtype=float)
    GH = np.zeros(vertical_layers, dtype=float)
    SH = np.zeros(vertical_layers, dtype=float)
    SM = np.zeros(vertical_layers, dtype=float)

    DTEF = np.zeros(vertical_layers, dtype=float)
    BPROD = np.zeros(vertical_layers, dtype=float)
    PROD = np.zeros(vertical_layers, dtype=float)
    SPROD = np.zeros(vertical_layers, dtype=float)

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
    von_karman_constant = 0.40  # KAPPA
    SQ = 0.2
    CIWC = 1.0
    gravity = 9.806
    # SM = KB*0.39
    # SH = KB*0.49
    # GM = KB*0.154
    # GH = KB*0.154

    DH = h
    for i in range(1, vertical_layers - 1):
        A[i] = -twice_the_timestep * (diffusion_coefficient_kinetic_energy[i + 1] +
                                      diffusion_coefficient_kinetic_energy[i] + 2 * umol) * \
               0.5 / (vertical_spacing_staggered[i - 1] * vertical_spacing[i] * h * h)
        C[i] = -twice_the_timestep * (diffusion_coefficient_kinetic_energy[i - 1] +
                                      diffusion_coefficient_kinetic_energy[i] + 2 * umol) * \
               0.5 / (vertical_spacing_staggered[i - 1] * vertical_spacing[i - 1] * h * h)

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
    kinetic_energy_forward[vertical_layers - 1] = \
        0.5 * np.sqrt((bottom_stress_zonal + bottom_stress_zonal) ** 2 +
                      (bottom_stress_meridional + bottom_stress_meridional) ** 2) * CONST1

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   CALCULATE PRESSURE IN UNITS OF DECIBARS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    for i in range(0, vertical_layers - 1):
        pressure = -gravity * 1.025 * vertical_coordinates_staggered[i] * h * .1
        CC[i] = 1449.2 + 1.34 * (salinity[i] - 35.) + 4.55 * temperature[i] - 0.045 * temperature[i] ** 2 + \
                0.00821 * pressure + (15.0 ** 1.e-9 * pressure ** 2)
        TEMP1[i] = 2./CC[i]
        TEMP2[i] = (0.00821*pressure)
        TEMP3[i] = (1.-0.40 * (pressure/CC[i]**2))

    for i in range(0, vertical_layers - 1):
        pressure = -gravity * 1.025 * vertical_coordinates_staggered[i] * h * .1
        CC[i] = CC[i] * (1. - TEMP1[i] * (TEMP2[i] + 15. * 1.e-9 * pressure ** 2) * TEMP3[i]) ** (-0.5)
        CC[i] = 1449.1 + .00821 * pressure + 4.55 * temperature[i] - .045 * temperature[i] ** 2 + 1.34 * (salinity[i] - 35.)
        CC[i] = CC[i] / np.sqrt((1. - .01642 * pressure / CC[i]) * (1. - 0.40 * pressure / CC[i] ** 2))

    for i in range(1, vertical_layers - 1):
        kinetic_energy_backward[i] = np.abs(kinetic_energy_backward[i])
        kinetic_energy_times_length_backward[i] = np.abs(kinetic_energy_times_length_backward[i])
        BOYGR[i] = gravity * (density_profile[i - 1] - density_profile[i]) / (vertical_spacing_staggered[i - 1] * h)
        DTEF[i] = kinetic_energy_backward[i] * np.sqrt(kinetic_energy_backward[i]) / \
                  (B1 * kinetic_energy_times_length_backward[i] + SMALL)
        SPROD[i] = .25 * diffusion_coefficient_momentum[i] * \
                   ((velocity_zonal[i] + velocity_zonal[i] - velocity_zonal[i - 1] - velocity_zonal[i - 1]) ** 2
                    + (velocity_meridional[i] + velocity_meridional[i] - velocity_meridional[i - 1] - velocity_meridional[i - 1]) ** 2) / \
                   (vertical_spacing_staggered[i - 1] * h) ** 2 * CIWC ** 2
        BPROD[i] = diffusion_coefficient_tracers[i] * BOYGR[i]
        PROD[i] = SPROD[i] + BPROD[i]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP DOWNWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(1, vertical_layers - 1):
        VHP[i] = 1. / (A[i] + C[i] * (1. - VH[i - 1]) - (2. * twice_the_timestep * DTEF[i] + 1.))
        VH[i] = A[i] * VHP[i]
        VHP[i] = (-2. * twice_the_timestep * PROD[i] + C[i] * VHP[i - 1] - kinetic_energy_backward[i]) * VHP[i]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP UPWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(0, vertical_layers - 1):  # 104
        k = vertical_layers - i
        kinetic_energy_forward[k] = VH[k] * kinetic_energy_forward[k + 1] + VHP[k]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SEECTION SOLVES FOR TEH EQUATION
    #   DT2(KQ*Q2L')' - Q2L*(DT2*DTEF+1.) = -Q2LB
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   BOUNDARY CONDITIONS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    VH[0] = 0.
    VHP[0] = 0.
    kinetic_energy_times_length_forward[vertical_layers - 1] = 0.

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP DOWNWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(1, vertical_layers - 1):
        DTEF[i] = DTEF[i] * (1. + E2 * ((1. / np.abs(vertical_coordinates[i] - vertical_coordinates[0]) + 1. / np.abs(vertical_coordinates[i] - vertical_coordinates[vertical_layers])) * length_scale[i] / (h * von_karman_constant)) ** 2)
        VHP[i] = 1. / (A[i] + C[i] * (1. - VH[i - 1]) - (twice_the_timestep * DTEF[i] + 1.))
        VH[i] = A[i] * VHP[i]
        VHP[i] = (twice_the_timestep * (- (SPROD[i] + E3 * BPROD[i]) * length_scale[i] * E1)
                  + C[i] * VHP[i - 1] - kinetic_energy_times_length_backward[i]) * VHP[i]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP UPWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(0, vertical_layers - 1):
        k = vertical_layers - i
        kinetic_energy_times_length_forward[k] = VH[k] * kinetic_energy_times_length_forward[k + 1] + VHP[k]

    for i in range(1, vertical_layers - 1):
        if kinetic_energy_forward[i] > SMALL or kinetic_energy_times_length_forward[i] > SMALL:
            break
        kinetic_energy_forward[i] = SMALL
        kinetic_energy_times_length_forward[i] = SMALL

    for i in range(0, vertical_layers - 1):
        kinetic_energy_forward[i] = np.abs(kinetic_energy_forward[i])
        kinetic_energy_times_length_forward[i] = np.abs(kinetic_energy_times_length_forward[i])

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

    length_scale[0] = 0.
    length_scale[vertical_layers - 1] = 0.
    GH[0] = 0.
    GH[vertical_layers - 1] = 0.

    for i in range(1, vertical_layers - 1):
        length_scale[i] = kinetic_energy_times_length_forward[i] / kinetic_energy_forward[i]
        GH[i] = length_scale[i] ** 2 / kinetic_energy_forward[i] * BOYGR[i]

    for i in range(0, vertical_layers):
        GH[i] = np.mininimum(GH[i], .028)
        SH[i] = COEF1 / (1. - COEF2 * GH[i])
        SM[i] = COEF3 + SH[i] * COEF4 * GH[i]
        SM[i] = SM[i] / (1. - COEF5 * GH[i])

    for i in range(0, vertical_layers):
        KN[i] = length_scale[i] * np.sqrt(np.abs(kinetic_energy[i]))
        diffusion_coefficient_kinetic_energy[i] = (KN[i] * .41 * SM[i] + diffusion_coefficient_kinetic_energy[i]) * .5
        #   KQ[K]= (KN[K] * .41 * SH[K] + KQ[K]) * .5
        diffusion_coefficient_momentum[i] = (KN[i] * SM[i] + diffusion_coefficient_momentum[i]) * .5
        diffusion_coefficient_tracers[i] = (KN[i] * SH[i] + diffusion_coefficient_tracers[i]) * .5

    return kinetic_energy_forward, kinetic_energy_times_length_forward, \
           diffusion_coefficient_momentum, diffusion_coefficient_tracers, diffusion_coefficient_kinetic_energy


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
def calculate_vertical_temperature_and_salinity_profiles(vertical_spacing, vertical_spacing_staggered, vertical_coordinates, vertical_coordinates_staggered,
                                                         diffusion_coefficient_tracers, FF, WFSURF, WFBOT, shortwave_radiation, FSURF, nbc):

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

    A = np.zeros(vertical_layers, dtype=float)
    C = np.zeros(vertical_layers, dtype=float)
    VH = np.zeros(vertical_layers, dtype=float)
    VHP = np.zeros(vertical_layers, dtype=float)

    # SW PROFILE
    vertical_radiation_profile = np.empty(vertical_layers, dtype=float)

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
    for i in range(1, vertical_layers - 1):
        A[i - 1] = -twice_the_timestep * (diffusion_coefficient_tracers[i] + umol) / (vertical_spacing[i - 1] * vertical_spacing_staggered[i - 1] * h * h)
        C[i] = -twice_the_timestep * (diffusion_coefficient_tracers[i] + umol) / (vertical_spacing[i] * vertical_spacing_staggered[i - 1] * h * h)

    vertical_radiation_profile[:] = zero

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SURFACE BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # *** PENETRATIVE RADIATION CALCULATION. AT THE BOTTOM ANY UNATTENUATED IS DEPOSITED IN THE BOTTOM LAYER.

    if nbc == 1:
        VH[0] = A[0] / (A[0] - one)
        VHP[0] = -twice_the_timestep * (WFSURF + shortwave_radiation) / (-vertical_spacing[0] * h) - FF[0]
        VHP[0] = VHP[0] / (A[0] - one)

    elif nbc == 2:
        vertical_radiation_profile[:] = shortwave_radiation * (RP[ntp] * np.exp(vertical_coordinates[:] * h / AD1[ntp]) + (one - RP[ntp] * np.exp(vertical_coordinates[:] * h / AD2[ntp])))  # ***
        vertical_radiation_profile[vertical_layers - 1] = zero

        VH[0] = A[0] / (A[0] - one)
        VHP[0] = twice_the_timestep * (WFSURF + vertical_radiation_profile[0] - vertical_radiation_profile[1]) / (vertical_spacing[0] * h) - FF[0]
        VHP[0] = VHP[0] / (A[0] - one)

    elif nbc == 3:
        VH[0] = zero
        VHP[0] = FSURF

    elif nbc == 4:
        vertical_radiation_profile[:] = shortwave_radiation * (RP[ntp] * np.exp(vertical_coordinates[:] * h / AD1[ntp]) + (one - RP[ntp] * np.exp(vertical_coordinates[:] * h / AD2[ntp])))  # ***
        vertical_radiation_profile[vertical_layers - 1] = 0

        VH[0] = zero
        VHP[0] = FSURF

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SECTION SOLVES THE EQUATION
    #   DT2*(KH*FF')' -FF = -FB
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(1, vertical_layers - 2):
        VHP[i] = 1 / (A[i] + C[i] * (1 - VH[i - 1]) - 1)
        VH[i] = A[i] * VHP[i]
        VHP[i] = (C[i] * VHP[i - 1] - FF[i] + twice_the_timestep * (vertical_radiation_profile[i] - vertical_radiation_profile[i + 1]) / (h * vertical_spacing[i])) * VHP[i]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   APPLY A NON ADIABATIC BOTTOM BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    FF[vertical_layers - 2] = (C[vertical_layers - 2] * VHP[vertical_layers - 3] - FF[vertical_layers - 2] + (WFBOT * twice_the_timestep / (vertical_spacing[vertical_layers - 2] * h))) / (
            C[vertical_layers - 2] * (1 - VH[vertical_layers - 3]) - 1)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   APPLY A NON ADIABATIC BOTTOM BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(1, vertical_layers - 1):
        k = vertical_layers - i
        FF[k] = VH[k] * FF[k + 1] + VHP[k]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   ZEROING
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    VH[:] = zero
    VHP[:] = zero
    A[:] = zero
    C[:] = zero

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
def calculate_vertical_zonal_velocity_profile(vertical_spacing, vertical_spacing_staggered, wind_stress_zonal,
                                              diffusion_coefficient_momentum, velocity_zonal_forward):

    A = np.zeros(vertical_layers, dtype=float)
    C = np.zeros(vertical_layers, dtype=float)
    VH = np.zeros(vertical_layers, dtype=float)
    VHP = np.zeros(vertical_layers, dtype=float)
    # UMOL1 = 0.0007

    for i in range(1, vertical_layers - 1):
        A[i - 1] = -twice_the_timestep * (diffusion_coefficient_momentum[i] + umol) / \
                   (vertical_spacing[i - 1] * vertical_spacing_staggered[i - 1] * h * h)
        C[i] = -twice_the_timestep * (diffusion_coefficient_momentum[i] + umol) / \
               (vertical_spacing[i] * vertical_spacing_staggered[i - 1] * h * h)

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-twice_the_timestep * wind_stress_zonal / (-vertical_spacing[0] * h) - velocity_zonal_forward[0]) / (A[0] - 1.)

    for i in range(1, vertical_layers - 2):
        VHP[i] = 1. / (A[i] + C[i] * (1. - VH[i - 1]) - 1.)
        VH[i] = A[i] * VHP[i]
        VHP[i] = (C[i] * VHP[i - 1] - velocity_zonal_forward[i]) * VHP[i]

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-twice_the_timestep * wind_stress_zonal / (-vertical_spacing[0] * h) - velocity_zonal_forward[0]) / (A[0] - 1.)

    # IF(NO_BOT_STRESS)THEN
    CBC = 0.0
    # ELSE
    #    CBC = MAX(.0025,.16/LOG((ZZ(KB-1)-Z(KB))*h/.01)**2)
    #    CBC = CBC*SQRT(UB(KB-1)**2+ (.25* (VB(KB-1)+VB(KB-1)+VB(KB-1)+ &
    #        VB(KB-1)))**2)
    # ENDIF
    #
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    velocity_zonal_forward[vertical_layers - 2] = (C[vertical_layers - 2] * VHP[vertical_layers - 3] - velocity_zonal_forward[vertical_layers - 2]) / (
            CBC * twice_the_timestep / (-vertical_spacing[vertical_layers - 2] * h) - 1. - (VH[vertical_layers - 3] - 1.) * C[vertical_layers - 2])
    for i in range(1, vertical_layers - 1):
        KI = vertical_layers - i
        velocity_zonal_forward[KI - 1] = VH[KI - 1] * velocity_zonal_forward[KI] + VHP[KI - 1]

    # WUBOT = -CBC * UF[KB - 2]  # 92
    for i in range(0, vertical_layers):
        VH[i] = 0.
        VHP[i] = 0.
        A[i] = 0.
        C[i] = 0.

    return velocity_zonal_forward


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
def calculate_vertical_meridional_velocity_profile(vertical_spacing, vertical_spacing_staggered, wind_stress_meridional,
                                                   diffusion_coefficient_momentum, velocity_meridional_forward):

    A = np.zeros(vertical_layers, dtype=float)
    C = np.zeros(vertical_layers, dtype=float)
    VH = np.zeros(vertical_layers, dtype=float)
    VHP = np.zeros(vertical_layers, dtype=float)

    for i in range(1, vertical_layers - 1):
        A[i - 1] = -twice_the_timestep * (diffusion_coefficient_momentum[i] + umol) / \
                   (vertical_spacing[i - 1] * vertical_spacing_staggered[i - 1] * h * h)
        C[i] = -twice_the_timestep * (diffusion_coefficient_momentum[i] + umol) / \
               (vertical_spacing[i] * vertical_spacing_staggered[i - 1] * h * h)

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-twice_the_timestep * wind_stress_meridional / (-vertical_spacing[0] * h) - velocity_meridional_forward[0]) / (A[0] - 1.)

    # 98 CONTINUE

    for i in range(1, vertical_layers - 2):
        VHP[i] = 1. / (A[i] + C[i] * (1. - VH[i - 1]) - 1.)
        VH[i] = A[i] * VHP[i]
        VHP[i] = (C[i] * VHP[i - 1] - velocity_meridional_forward[i]) * VHP[i]

    # IF(NO_BOT_STRESS)THEN
    CBC = 0.0
    # ELSE
    # 104      CBC = MAX(.0025,.16/LOG((ZZ(KB-1)-Z(KB))*h/.01)**2)
    #          CBC = CBC*SQRT((.25* (UB(KB-1)+UB(KB-1)+UB(KB-1)+UB(KB-1)))**2 + VB(KB-1)**2)
    # ENDIF
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   TO RESTORE BOTTOM B.L. DELETE NEXT LINE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    velocity_meridional_forward[vertical_layers - 2] = (C[vertical_layers - 1] * VHP[vertical_layers - 3] - velocity_meridional_forward[vertical_layers - 2]) / (
            CBC * twice_the_timestep / (-vertical_spacing[vertical_layers - 2] * h) - 1. - (VH[vertical_layers - 3] - 1.) * C[vertical_layers - 2])

    for i in range(1, vertical_layers - 1):
        k = vertical_layers - i
        velocity_meridional_forward[k] = VH[k] * velocity_meridional_forward[k + 1] + VHP[k]

    # WVBOT = -CBC * VF[KB - 2]  # 92
    for i in range(0, vertical_layers):
        VH[i] = 0.
        VHP[i] = 0.
        A[i] = 0.
        C[i] = 0.

    return velocity_meridional_forward


# EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
