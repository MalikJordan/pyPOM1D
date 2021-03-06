# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from cppdefs import *
# from pom.modules import vertical_layers, zero, one, seconds_per_day
import numpy as np
from inputs.params_POMBFM import *
twice_the_timestep = 2. * dti

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# ROUTINE: adverte
#
# DESCRIPTION:  SUBROUTINE TO HANDLE THE SINKING OF BFM STATE VAR'S SINKING IS TREATED AS DOWNWARD VERTICAL ADVECTION
#               COMPUTED WITH UPSTREAM FINITE DIFFERENCES.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def adverte(FB,F,FF,W,vertical_spacing_reciprocal):

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL VARIABLES
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # FB, F, FF = np.empty(KB,dtype=float)
    # W = np.empty(KB,dtype=float)
    # DTI2 = float()
    # K = int()

    F[vertical_layers-1] = F[vertical_layers-2]
    FB[vertical_layers-1] = FB[vertical_layers-2]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Calculate vertical advection. Mind downward velocities are negative
    #   Upwind scheme:
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    FF[0] = vertical_spacing_reciprocal[0] * F[0] * W[1]

    for i in range(1,vertical_layers-1):
        FF[i] = vertical_spacing_reciprocal[i] * (F[i]*W[i+1] - F[i-1]*W[i])

    return FF


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# ROUTINE: CALCDEPTH
#
# DESCRIPTION:  This subroutine establishes the vertical resolution log distributions
#               at the top and bottom, and a linear distribution between KL1 and KL2.
#               Default values: KL1 = .3*KB AND KL2 = KB-2.
#               Yields a log distribution at the top and none at the bottom.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def create_vertical_coordinate_system(vertical_layers, surface_layers_with_log_distribution, bottom_layers_with_log_distribution):

    # KL1 = surface_log_distribution
    # KL2 = bottom_log_distribution
    # KB = vertical_layers

    vertical_coordinates = np.zeros(vertical_layers,dtype=float)
    vertical_coordinates_staggered = np.zeros(vertical_layers,dtype=float)
    vertical_spacing = np.zeros(vertical_layers, dtype=float)
    vertical_spacing_staggered = np.zeros(vertical_layers, dtype=float)

    surface_logspace_layers = surface_layers_with_log_distribution - 2.
    bottom_logspace_layers = vertical_layers - bottom_layers_with_log_distribution - 1.

    initial_spacing = (bottom_layers_with_log_distribution - surface_layers_with_log_distribution + 4.) * 2 ** (-6.3)

    vertical_coordinates_staggered[0] = 0.5 * initial_spacing

    for i in range(1, int(surface_layers_with_log_distribution) - 1):
        vertical_coordinates[i] = -initial_spacing * 2**(i-2)
        vertical_coordinates_staggered[i] = -initial_spacing * 2**(i-1.5)

    for i in range(int(surface_layers_with_log_distribution) - 1, vertical_layers):
        vertical_coordinates[i]  = -(i - surface_logspace_layers) / (bottom_layers_with_log_distribution
                                                                     - surface_layers_with_log_distribution + 4.)
        vertical_coordinates_staggered[i] = -(i - surface_logspace_layers + 0.5) / (bottom_layers_with_log_distribution
                                                                                    - surface_layers_with_log_distribution + 4.)

    for i in range(0,vertical_layers-1):
        vertical_spacing[i] = vertical_coordinates[i] - vertical_coordinates[i+1]
        vertical_spacing_staggered[i] = vertical_coordinates_staggered[i] - vertical_coordinates_staggered[i+1]

    return vertical_coordinates, vertical_coordinates_staggered, vertical_spacing, vertical_spacing_staggered



# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# ROUTINE: CalcLightDistribution
#
# DESCRIPTION
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def CalcLightDistribution():

    EIR[0] = EIR[0] * p_PAR / E2W
    for BoxNumberZ in range(1, NO_BOXES_Z):
        box_no = BoxNumberZ
        EIR[box_no] = EIR[box_no - 1] * np.exp(- 1.0 * xEPS[box_no - 1] * Depth[box_no - 1])

    return EIR


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# ROUTINE: DENS
#
# DESCRIPTION:  This subroutine computes density.
#               T = Potential temperature
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def calculate_vertical_density_profile(temperature,salinity,vertical_spacing_staggered,timestep,vertical_layers):

    vertical_density_profile = np.zeros(vertical_layers,dtype=float)
    gravity = 9.806

    for i in range(0,vertical_layers-1):

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   APPROXIMATE PRESSURE IN UNITS OF BARS
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        pressure = -gravity * 1.025 * vertical_spacing_staggered[i] * timestep * 0.01
        density = 999.842594 + 6.793952E-2 * temperature[i] - 9.095290E-3 * temperature[i] ** 2 + \
               1.001685E-4 * temperature[i] ** 3 - 1.120083E-6 * temperature[i] ** 4 + 6.536332E-9 * temperature[i] ** 5
        density = density + (0.824493 - 4.0899E-3 * temperature[i] + 7.6438E-5 * temperature[i] ** 2 -
                             8.2467E-7 * temperature[i] ** 3 + 5.3875E-9 * temperature[i] ** 4) * salinity[i] + \
                  (-5.72466E-3 + 1.0227E-4 * temperature[i] - 1.6546E-6 * temperature[i] ** 2) * \
                  (np.abs(salinity[i])) ** 1.5 + 4.8314E-4 * salinity[i] ** 2
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   FOR SHALLOW WATER THE PRESSURE DEPENDENCY CAN BE NEGLECTED
        #   IN WHICH IT SHOULD ALSO BE OMITTED IN PROFQ
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        vertical_density_profile[i] = (density - 1000.) * 1.E-3

    vertical_density_profile[vertical_layers-1] = vertical_density_profile[vertical_layers-2]

    return vertical_density_profile




# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# ROUTINE: LF1D
#
# DESCRIPTION:
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# ROUTINE: MLDPTH
#
# DESCRIPTION:  This subroutine computes calculates the mixed layer depth.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def calculate_mixed_layer_depth(vertical_coordinates_staggered, temperature, vertical_layers):

    mixed_layer_depth = np.zeros(vertical_layers,dtype=float)
    zero = 1.E-06

    for i in range(0,vertical_layers-1):

        if temperature[0] > temperature[i]+0.2:
            break

        mixed_layer_depth[i] = vertical_coordinates_staggered[i] - (temperature[i] + 0.2 - temperature[0]) * \
                            (vertical_coordinates_staggered[i] - vertical_coordinates_staggered[i+1]) / \
                               (temperature[i] - temperature[i+1] + zero)

    return mixed_layer_depth


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# ROUTINE: POM_DIA_BFM
#
# DESCRIPTION:  This routine calculates means and writes the output in diagnostic mode writes also the restart
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def pom_dia_bfm(kt,TT):

    # need netcdf_bfm stuff

    # TIME IN SECONDS
    localtime = float()

    # SAVING FREQUENCY
    time_to_save = int()  # time in seconds

    # TIME ELLAPSED (IN SECONDS) SINCE THE BEGINNING OF THE SIMULATION
    localtime = TIME * seconds_per_day

    # SAVING FREQUENCY IN TIME MARCHING LOOP ITERATIONS
    time_to_save = int(out_delta*seconds_per_day)

    # SUMMING UP THE FIELDS TO BE SAVED
    calcmean_bfm(ACCUMULATE)

    # WRITE OUTPUT
    if TT + (dti/seconds_per_day) > out_delta:

        calcmean_bfm(MEAN)
        save_bfm(localtime)

        # RESET TIME COUNTER
        TT = TT - int(TT)

    # WHEN RUNNING IN COUPLING WITH POM RESTART IS WRITTEN IN MAIN AFTER THE END
    # OF THE TIME MARCHING LOOP

    # WRITE RESTART
    if kt >= IEND:

        if -TT < dti/seconds_per_day:
            save_rst_bfm(localtime)
            close_ncdf(ncid_rst)
            close_ncdf(ncid_bfm)
            print('POM_DIA: NETCDF RESTART WRITTEN, TIME--> ', TIME, kt, IEND, DTI, TT)

    return


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# ROUTINE: vdiff_SOS
#
# DESCRIPTION:  This routine calculates the vertical diffusivity of BFM biochemical components and
#               integrates BFM state var's with Source Splitting (SoS) method.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def vdiff_SOS():

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL VARIABLES
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # COUNTER & FLAGS
    # K, M, N, NBC = float()

    # BFM STATE VAR. @ time t-DTI, t, T+DTI RESPECTIVELY
    # fbio, ffbio, fbbio = np.empty(KB,dtype=float)

    # SURFACE FLUX STORAGE FOR NUT'S O2 & CO2
    # surflux = float()
    # botflux = float()

    # RELAXATION VELOCITY FOR NUT'S
    # trelax_o2o = float()
    # trelax_n1p = float()
    # trelax_n3n = float()
    # trelax_n4n = float()

    # SEDIMENTATION VELOCITY
    # sink = np.empty(KB,dtype=float)
    # POCsink = float()
    # W1R6 = float()

    # The input general cir. vertical vel. is suppose to be in m/s
    W_ON = 1.0
    # The input eddy vertical vel. is provided in m/d
    Weddy_ON = 0.1/86400.0  # to m/s

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL VARIABLES
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    trelax_o2o = nrt_o2o / seconds_per_day
    trelax_n1p = nrt_n1p / seconds_per_day
    trelax_n3n = nrt_n3n / seconds_per_day
    trelax_n4n = nrt_n4n

    # LOOP OVER BFM STATE VAR'S
    for M in range(0,NO_D3_BOX_STATES):

        # ZEROING

        surflux = zero
        botflux = zero
        fbio[:] = zero
        fbbio[:] = zero
        ffbio[:] = zero
        sink[:] = zero
        POCsink = zero

        # LOAD BFM STATE VAR.
        for i in range(0,vertical_layers-1):
            fbio[i] = D3STATE[M][i]
            fbbio[i] = D3STATEB[M][i]

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   NUTRIENTS SURFACE AND BOTTOM FLUXES
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        if M == ppO2o:
            surflux = -(jsurO2o[0] / seconds_per_day)
            botflux = (o2o[vertical_layers-2] - O2BOTT) * trelax_o2o
        elif M == ppO3c:
            surflux = zero
        elif M == ppN1p:
            surflux = zero
            botflux = (n1p[vertical_layers-2] - PO4BOTT) * trelax_n1p
        elif M == ppN3n:
            surflux = zero
            botflux = (n3n[vertical_layers-2] - NO3BOTT) * trelax_n3n
        elif M == ppN5s:
            surflux = zero
        else:
            surflux = zero

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   BOTTOM FLUX
        #   R1: Dissolved Organic Matter
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        # The botflux for Dissolved Organic Matter is left equal to ZERO

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   BOTTOM FLUX
        #   R6: Particulate Organic Matter
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        # The botflux for Particulate Organic Matter is left equal to ZERO

        if M >= ppR6c and M <= ppR6s:

            for i in range(0,vertical_layers-1):
                sink[i] = sink[i] - sediR6[i]/seconds_per_day

            # FINAL SINK VALUE
            sink[vertical_layers-1] = sink[vertical_layers-2]

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   SEDIMENTATION PHYTOPLANKTON
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        # The botflux for Phytoplankton is left equal to ZERO

        if M >= ppP1c and M <= ppP4l:

            if M in range(ppP1c,ppP1s):
                N = iiP1
            elif M in range(ppP2c,ppP2l):
                N = iiP2
            elif M in range(ppP3c,ppP3l):
                N = iiP3
            elif M in range(ppP4c,ppP4l):
                N = iiP4

            for i in range(0,vertical_layers-1):
                sink[i] = sink[i] - sediPPY[i]/seconds_per_day

            # FINAL SINK VALUE
            sink[vertical_layers - 1] = sink[vertical_layers - 2]

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   BOTTOM FLUX
        #   Z: Zooplankton
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        # The bot flux for Zooplankton is left equal to ZERO

        # SINKING: UPSTREAM VERTICAL ADVECTION
        adverte(fbbio,fbio,ffbio,sink)

        # SOURCE SPLITTING LEAPFROG INTEGRATION
        for i in range(0,vertical_layers-1):
            ffbio[i] = ffbio[i] + twice_the_timestep*((ffbio[i]/H) + D3SOURCE[M][i])

        # COMPUTE VERTICAL DIFFUSION AND TERMINATE INTEGRATION
        # IMPLICIT LEAPFROGGING
        PROFTS(ffbio,surflux,botflux,zero,zero,nbcbfm,twice_the_timestep,ntp,umolbfm)

        # CLIPPING......IF NEEDED
        for i in range(0,vertical_layers-1):
            ffbio[i] = max(p_small,ffbio[i])

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   Mix the time step and restore time sequence
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        for N in range(0,vertical_layers-1):
            D3STATEB[M][N] = fbio[N] + 0.5*smoth*(ffbio[N] + fbbio[N] - 2.*fbio[N])
            D3STATE[M][N] = ffbio[N]

    if AssignAirPelFluxesInBFMFlag is False:
        jsurO2o[:] = zero
        jsurO3c[:] = zero

    return


# EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
