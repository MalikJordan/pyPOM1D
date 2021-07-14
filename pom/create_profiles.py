from cppdefs import *
from inputs import params_POMBFM
import numpy as np
from pom.pom_constants import vertical_layers, twice_the_timestep

np.set_printoptions(precision=16)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: PROFE
#
# DESCRIPTION:  This subroutine solves for vertical diffusivity.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def create_vertical_diffusivity_profile(vertical_grid, diffusion, property):

    A = np.zeros(vertical_layers)
    C = np.zeros(vertical_layers)
    VH = np.zeros(vertical_layers)
    VHP = np.zeros(vertical_layers)

    vertical_diffusivity = np.zeros(vertical_layers)                               # FF

    umolpr = 1.E-05

    for i in range(1, vertical_layers - 1):
        A[i - 1] = -twice_the_timestep * (diffusion.tracers[i] + umolpr) / (vertical_grid.vertical_spacing[i - 1] * vertical_grid.vertical_spacing_staggered[i - 1] * params_POMBFM.h * params_POMBFM.h)
        print(A[i-1])
        C[i] = -twice_the_timestep * (diffusion.tracers[i] + umolpr) / (vertical_grid.vertical_spacing[i] * vertical_grid.vertical_spacing_staggered[i - 1] * params_POMBFM.h * params_POMBFM.h)

    # need condition for nbc for VH[0] and VHP[0]
    VH[0] = A[0] / (A[0]-1.)
    VHP[0] = -twice_the_timestep * property.surface_flux / (-vertical_grid.vertical_spacing[0] * params_POMBFM.h) - vertical_diffusivity[0]
    VHP[0] = VHP[0] / (A[0]-1.)

    VH[0] = 0.
    VHP[0] = property.surface_value

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
def create_kinetic_energy_profile(vertical_grid, diffusion, temperature, salinity, density_profile, velocity_zonal, velocity_meridional,
                                  kinetic_energy, kinetic_energy_times_length, wind_stress, bottom_stress):

    # INITIALIZE VARIABLES
    A = np.zeros(vertical_layers)
    C = np.zeros(vertical_layers)
    VH = np.zeros(vertical_layers)
    VHP = np.zeros(vertical_layers)

    KN = np.zeros(vertical_layers)
    GH = np.zeros(vertical_layers)
    SH = np.zeros(vertical_layers)
    SM = np.zeros(vertical_layers)

    DTEF = np.zeros(vertical_layers)
    BPROD = np.zeros(vertical_layers)
    PROD = np.zeros(vertical_layers)
    SPROD = np.zeros(vertical_layers)
    pressure = 0.

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL ARRAYS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    BOYGR = np.empty(vertical_layers); CC = np.empty(vertical_layers)
    TEMP1 = np.empty(vertical_layers); TEMP2 = np.empty(vertical_layers); TEMP3 = np.empty(vertical_layers)

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
    SMALL = 1.E-08
    # SM = KB*0.39
    # SH = KB*0.49
    # GM = KB*0.154
    # GH = KB*0.154

    for i in range(1, vertical_layers - 1):
        A[i] = -twice_the_timestep * (diffusion.kinetic_energy[i + 1] +
                                      diffusion.kinetic_energy[i] + 2 * params_POMBFM.umol) * \
               0.5 / (vertical_grid.vertical_spacing_staggered[i - 1] * vertical_grid.vertical_spacing[i] * params_POMBFM.h * params_POMBFM.h)
        C[i] = -twice_the_timestep * (diffusion.kinetic_energy[i - 1] +
                                      diffusion.kinetic_energy[i] + 2 * params_POMBFM.umol) * \
               0.5 / (vertical_grid.vertical_spacing_staggered[i - 1] * vertical_grid.vertical_spacing[i - 1] * params_POMBFM.h * params_POMBFM.h)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SECTION SOLVES FOR THE EQUATION
    #   DT2*(KQ*Q2')' - Q2*(2.*DT2*DTEF+1.) = -Q2B
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    CONST1 = 16.6 ** 0.6666667 * CIWC

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   BOUNDARY CONDITIONS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    VH[0] = 0.0
    VHP[0] = np.sqrt(wind_stress.zonal ** 2 + wind_stress.meridional ** 2) * CONST1
    kinetic_energy.forward[vertical_layers - 1] = 0.5 * np.sqrt((bottom_stress.zonal + bottom_stress.zonal) ** 2 +
                                                                (bottom_stress.meridional + bottom_stress.meridional) ** 2) * CONST1

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   CALCULATE PRESSURE IN UNITS OF DECIBARS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # for i in range(0, vertical_layers - 1):
    #     pressure = -gravity * 1.025 * vertical_grid.vertical_coordinates_staggered[i] * params_POMBFM.h * .1
    #     CC[i] = 1449.2 + 1.34 * (salinity.current[i] - 35.) + 4.55 * temperature.current[i] - 0.045 * temperature.current[i] ** 2 + \
    #             0.00821 * pressure + (15.0 ** 1.e-9 * pressure ** 2)
    #     TEMP1[i] = 2./CC[i]
    #     TEMP2[i] = (0.00821*pressure)
    #     TEMP3[i] = (1.-0.40 * (pressure/CC[i]**2))
    #
    for i in range(0, vertical_layers - 1):
        CC[i] = CC[i] * (1. - TEMP1[i] * (TEMP2[i] + 15. * 1.E-9 * pressure ** 2) * TEMP3[i]) ** (-0.5)
        pressure = -gravity * 1.025 * vertical_grid.vertical_coordinates_staggered[i] * params_POMBFM.h * .1
        CC[i] = 1449.1 + .00821 * pressure + 4.55 * temperature.current[i] - .045 * temperature.current[i] ** 2 + 1.34 * (salinity.current[i] - 35.)
        CC[i] = CC[i] / np.sqrt((1. - .01642 * pressure / CC[i]) * (1. - 0.40 * pressure / CC[i] ** 2))

    for i in range(1, vertical_layers - 1):
        kinetic_energy.backward[i] = np.abs(kinetic_energy.backward[i])
        kinetic_energy_times_length.backward[i] = np.abs(kinetic_energy_times_length.backward[i])
        BOYGR[i] = gravity * (density_profile[i - 1] - density_profile[i]) / (vertical_grid.vertical_spacing_staggered[i - 1] * params_POMBFM.h)
        DTEF[i] = kinetic_energy.backward[i] * np.sqrt(kinetic_energy.backward[i]) / \
                  (B1 * kinetic_energy_times_length.backward[i] + SMALL)
        SPROD[i] = .25 * diffusion.momentum[i] * \
                   ((velocity_zonal.current[i] + velocity_zonal.current[i] - velocity_zonal.current[i - 1] - velocity_zonal.current[i - 1]) ** 2
                    + (velocity_meridional.current[i] + velocity_meridional.current[i] - velocity_meridional.current[i - 1] - velocity_meridional.current[i - 1]) ** 2) / \
                   (vertical_grid.vertical_spacing_staggered[i - 1] * params_POMBFM.h) ** 2 * CIWC ** 2
        BPROD[i] = diffusion.tracers[i] * BOYGR[i]
        PROD[i] = SPROD[i] + BPROD[i]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP DOWNWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(1, vertical_layers - 1):
        VHP[i] = 1. / (A[i] + C[i] * (1. - VH[i - 1]) - (2. * twice_the_timestep * DTEF[i] + 1.))
        VH[i] = A[i] * VHP[i]
        VHP[i] = (-2. * twice_the_timestep * PROD[i] + C[i] * VHP[i - 1] - kinetic_energy.backward[i]) * VHP[i]
    # print(velocity_meridional.current[1])
    # print()
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP UPWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(2, vertical_layers+1):  # 104
        k = vertical_layers - i
        kinetic_energy.forward[k] = VH[k] * kinetic_energy.forward[k + 1] + VHP[k]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SEECTION SOLVES FOR TEH EQUATION
    #   DT2(KQ*Q2L')' - Q2L*(DT2*DTEF+1.) = -Q2LB
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   BOUNDARY CONDITIONS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    VH[0] = 0.
    VHP[0] = 0.
    kinetic_energy_times_length.forward[vertical_layers - 1] = 0.

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP DOWNWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(1, vertical_layers - 1):
        DTEF[i] = DTEF[i] * (1. + E2 * ((1. / np.abs(vertical_grid.vertical_coordinates[i] - vertical_grid.vertical_coordinates[0])
                                         + 1. / np.abs(vertical_grid.vertical_coordinates[i] - vertical_grid.vertical_coordinates[vertical_layers-1]))
                                        * vertical_grid.length_scale[i] / (params_POMBFM.h * von_karman_constant)) ** 2)
        VHP[i] = 1. / (A[i] + C[i] * (1. - VH[i - 1]) - (twice_the_timestep * DTEF[i] + 1.))
        VH[i] = A[i] * VHP[i]
        VHP[i] = (twice_the_timestep * (- (SPROD[i] + E3 * BPROD[i]) * vertical_grid.length_scale[i] * E1)
                  + C[i] * VHP[i - 1] - kinetic_energy_times_length.backward[i]) * VHP[i]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP UPWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(2, vertical_layers+1):
        k = vertical_layers - i
        kinetic_energy_times_length.forward[k] = VH[k] * kinetic_energy_times_length.forward[k + 1] + VHP[k]
    # print(kinetic_energy_times_length.forward[1])
    for i in range(1, vertical_layers - 1):
        if kinetic_energy.forward[i] > SMALL or kinetic_energy_times_length.forward[i] > SMALL:
            continue
        else:
            kinetic_energy.forward[i] = SMALL
            kinetic_energy_times_length.forward[i] = SMALL

    for i in range(0, vertical_layers - 1):
        kinetic_energy.forward[i] = np.abs(kinetic_energy.forward[i])
        kinetic_energy_times_length.forward[i] = np.abs(kinetic_energy_times_length.forward[i])

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

    vertical_grid.length_scale[0] = 0.
    vertical_grid.length_scale[vertical_layers - 1] = 0.
    GH[0] = 0.
    GH[vertical_layers - 1] = 0.

    for i in range(1, vertical_layers - 1):
        vertical_grid.length_scale[i] = kinetic_energy_times_length.forward[i] / kinetic_energy.forward[i]
        GH[i] = vertical_grid.length_scale[i] ** 2 / kinetic_energy.forward[i] * BOYGR[i]

    for i in range(0, vertical_layers):
        GH[i] = min(GH[i], .028)
        SH[i] = COEF1 / (1. - COEF2 * GH[i])
        SM[i] = COEF3 + SH[i] * COEF4 * GH[i]
        SM[i] = SM[i] / (1. - COEF5 * GH[i])

    for i in range(0, vertical_layers):
        KN[i] = vertical_grid.length_scale[i] * np.sqrt(np.abs(kinetic_energy.current[i]))
        diffusion.kinetic_energy[i] = (KN[i] * .41 * SM[i] + diffusion.kinetic_energy[i]) * .5
        #   KQ[K]= (KN[K] * .41 * SH[K] + KQ[K]) * .5
        diffusion.momentum[i] = (KN[i] * SM[i] + diffusion.momentum[i]) * .5
        diffusion.tracers[i] = (KN[i] * SH[i] + diffusion.tracers[i]) * .5

    return kinetic_energy, kinetic_energy_times_length, diffusion, vertical_grid


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
def calculate_vertical_temperature_and_salinity_profiles(vertical_grid,diffusion, property, shortwave_radiation, nbc):

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

    A = np.zeros(vertical_layers)
    C = np.zeros(vertical_layers)
    VH = np.zeros(vertical_layers)
    VHP = np.zeros(vertical_layers)

    # SW PROFILE
    vertical_radiation_profile = np.empty(vertical_layers)

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
        A[i - 1] = -twice_the_timestep * (diffusion.tracers[i] + params_POMBFM.umol) / (vertical_grid.vertical_spacing[i - 1] * vertical_grid.vertical_spacing_staggered[i - 1] * params_POMBFM.h * params_POMBFM.h)
        C[i] = -twice_the_timestep * (diffusion.tracers[i] + params_POMBFM.umol) / (vertical_grid.vertical_spacing[i] * vertical_grid.vertical_spacing_staggered[i - 1] * params_POMBFM.h * params_POMBFM.h)

    vertical_radiation_profile[:] = 0.

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SURFACE BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # *** PENETRATIVE RADIATION CALCULATION. AT THE BOTTOM ANY UNATTENUATED IS DEPOSITED IN THE BOTTOM LAYER.

    if nbc == 1:
        VH[0] = A[0] / (A[0] - 1.)
        VHP[0] = -twice_the_timestep * (property.surface_flux + shortwave_radiation) / (-vertical_grid.vertical_spacing[0] * params_POMBFM.h) - property.forward[0]
        VHP[0] = VHP[0] / (A[0] - 1.)

    elif nbc == 2:
        vertical_radiation_profile[:] = shortwave_radiation * (RP[params_POMBFM.ntp] * np.exp(vertical_grid.vertical_coordinates[:] * params_POMBFM.h / AD1[params_POMBFM.ntp]) + (1. - RP[params_POMBFM.ntp] * np.exp(vertical_grid.vertical_coordinates[:] * params_POMBFM.h / AD2[params_POMBFM.ntp])))  # ***
        vertical_radiation_profile[vertical_layers - 1] = 0.

        VH[0] = A[0] / (A[0] - 1.)
        VHP[0] = twice_the_timestep * (property.surface_flux + vertical_radiation_profile[0] - vertical_radiation_profile[1]) / (vertical_grid.vertical_spacing[0] * params_POMBFM.h) - property.forward[0]
        VHP[0] = VHP[0] / (A[0] - 1.)

    elif nbc == 3:
        VH[0] = 0.
        VHP[0] = property.surface_value

    elif nbc == 4:
        vertical_radiation_profile[:] = shortwave_radiation * (RP[params_POMBFM.ntp] * np.exp(vertical_grid.vertical_coordinates[:] * params_POMBFM.h / AD1[params_POMBFM.ntp]) + (1. - RP[params_POMBFM.ntp] * np.exp(vertical_grid.vertical_coordinates[:] * params_POMBFM.h / AD2[params_POMBFM.ntp])))  # ***
        vertical_radiation_profile[vertical_layers - 1] = 0.

        VH[0] = 0.
        VHP[0] = property.surface_value

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SECTION SOLVES THE EQUATION
    #   DT2*(KH*FF')' -FF = -FB
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(1, vertical_layers - 2):
        VHP[i] = 1 / (A[i] + C[i] * (1 - VH[i - 1]) - 1)
        VH[i] = A[i] * VHP[i]
        VHP[i] = (C[i] * VHP[i - 1] - property.forward[i] + twice_the_timestep * (vertical_radiation_profile[i] - vertical_radiation_profile[i + 1]) / (params_POMBFM.h * vertical_grid.vertical_spacing[i])) * VHP[i]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   APPLY A NON ADIABATIC BOTTOM BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    property.forward[vertical_layers - 2] = (C[vertical_layers - 2] * VHP[vertical_layers - 3] - property.forward[vertical_layers - 2] + (property.bottom_flux * twice_the_timestep / (vertical_grid.vertical_spacing[vertical_layers - 2] * params_POMBFM.h))) / (
            C[vertical_layers - 2] * (1 - VH[vertical_layers - 3]) - 1)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   APPLY A NON ADIABATIC BOTTOM BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(vertical_layers-2, 0, -1):
        # k = (vertical_layers) - i
        property.forward[i] = VH[i] * property.forward[i+1] + VHP[i]


    return property


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
def calculate_vertical_zonal_velocity_profile(vertical_grid, wind_stress, bottom_stress, diffusion, velocity_zonal):

    A = np.zeros(vertical_layers)
    C = np.zeros(vertical_layers)
    VH = np.zeros(vertical_layers)
    VHP = np.zeros(vertical_layers)
    # UMOL1 = 0.0007

    for i in range(1, vertical_layers - 1):
        A[i - 1] = -twice_the_timestep * (diffusion.momentum[i] + params_POMBFM.umol) / \
                   (vertical_grid.vertical_spacing[i - 1] * vertical_grid.vertical_spacing_staggered[i - 1] * params_POMBFM.h * params_POMBFM.h)
        C[i] = -twice_the_timestep * (diffusion.momentum[i] + params_POMBFM.umol) / \
               (vertical_grid.vertical_spacing[i] * vertical_grid.vertical_spacing_staggered[i - 1] * params_POMBFM.h * params_POMBFM.h)

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-twice_the_timestep * wind_stress.zonal / (-vertical_grid.vertical_spacing[0] * params_POMBFM.h) - velocity_zonal.forward[0]) / (A[0] - 1.)

    for i in range(1, vertical_layers - 2):
        VHP[i] = 1. / (A[i] + C[i] * (1. - VH[i - 1]) - 1.)
        VH[i] = A[i] * VHP[i]
        VHP[i] = (C[i] * VHP[i - 1] - velocity_zonal.forward[i]) * VHP[i]

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-twice_the_timestep * wind_stress.zonal / (-vertical_grid.vertical_spacing[0] * params_POMBFM.h) - velocity_zonal.forward[0]) / (A[0] - 1.)

    # IF(NO_BOT_STRESS)THEN
    CBC = 0.0
    # ELSE
    #    CBC = MAX(.0025,.16/LOG((ZZ(KB-1)-Z(KB))*h/.01)**2)
    #    CBC = CBC*SQRT(UB(KB-1)**2+ (.25* (VB(KB-1)+VB(KB-1)+VB(KB-1)+ &
    #        VB(KB-1)))**2)
    # ENDIF
    #
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    velocity_zonal.forward[vertical_layers - 2] = (C[vertical_layers - 2] * VHP[vertical_layers - 3] - velocity_zonal.forward[vertical_layers - 2]) / (
            CBC * twice_the_timestep / (-vertical_grid.vertical_spacing[vertical_layers - 2] * params_POMBFM.h) - 1. - (VH[vertical_layers - 3] - 1.) * C[vertical_layers - 2])
    for i in range(1, vertical_layers - 1):
        k = vertical_layers - 1 - i
        velocity_zonal.forward[k - 1] = VH[k - 1] * velocity_zonal.forward[k] + VHP[k - 1]
    # print(diffusion.momentum[1])
    bottom_stress.zonal = -CBC * velocity_zonal.forward[vertical_layers - 2]  # 92

    return velocity_zonal, bottom_stress


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
def calculate_vertical_meridional_velocity_profile(vertical_grid, wind_stress, bottom_stress, diffusion, velocity_meridional):

    A = np.zeros(vertical_layers)
    C = np.zeros(vertical_layers)
    VH = np.zeros(vertical_layers)
    VHP = np.zeros(vertical_layers)

    for i in range(1, vertical_layers - 1):
        A[i - 1] = -twice_the_timestep * (diffusion.momentum[i] + params_POMBFM.umol) / \
                   (vertical_grid.vertical_spacing[i - 1] * vertical_grid.vertical_spacing_staggered[i - 1] * params_POMBFM.h * params_POMBFM.h)
        C[i] = -twice_the_timestep * (diffusion.momentum[i] + params_POMBFM.umol) / \
               (vertical_grid.vertical_spacing[i] * vertical_grid.vertical_spacing_staggered[i - 1] * params_POMBFM.h * params_POMBFM.h)

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-twice_the_timestep * wind_stress.meridional / (-vertical_grid.vertical_spacing[0] * params_POMBFM.h) - velocity_meridional.forward[0]) / (A[0] - 1.)

    # 98 CONTINUE

    for i in range(1, vertical_layers - 2):
        VHP[i] = 1. / (A[i] + C[i] * (1. - VH[i - 1]) - 1.)
        VH[i] = A[i] * VHP[i]
        VHP[i] = (C[i] * VHP[i - 1] - velocity_meridional.forward[i]) * VHP[i]

    # IF(NO_BOT_STRESS)THEN
    CBC = 0.0
    # ELSE
    # 104      CBC = MAX(.0025,.16/LOG((ZZ(KB-1)-Z(KB))*h/.01)**2)
    #          CBC = CBC*SQRT((.25* (UB(KB-1)+UB(KB-1)+UB(KB-1)+UB(KB-1)))**2 + VB(KB-1)**2)
    # ENDIF
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   TO RESTORE BOTTOM B.L. DELETE NEXT LINE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    velocity_meridional.forward[vertical_layers - 2] = (C[vertical_layers - 1] * VHP[vertical_layers - 3] - velocity_meridional.forward[vertical_layers - 2]) / (
            CBC * twice_the_timestep / (-vertical_grid.vertical_spacing[vertical_layers - 2] * params_POMBFM.h) - 1. - (VH[vertical_layers - 3] - 1.) * C[vertical_layers - 2])

    for i in range(1, vertical_layers - 1):
        k = vertical_layers - 1 - i
        velocity_meridional.forward[k - 1] = VH[k - 1] * velocity_meridional.forward[k] + VHP[k - 1]

    bottom_stress.meridional = -CBC * velocity_meridional.forward[vertical_layers - 2]  # 92

    return velocity_meridional, bottom_stress


# EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
