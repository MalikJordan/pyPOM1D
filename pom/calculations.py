# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from cppdefs import *
# from pom.modules import *
import numpy as np

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
def adverte(FB,F,FF,W):

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL VARIABLES
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # FB, F, FF = np.empty(KB,dtype=float)
    # W = np.empty(KB,dtype=float)
    # DTI2 = float()
    # K = int()

    F[KB-1] = F[KB-2]
    FB[KB-1] = FB[KB-2]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Calculate vertical advection. Mind downward velocities are negative
    #   Upwind scheme:
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    FF[0] = DZR[0] * F[0] * W[1]

    for K in range(1,KB-1):
        FF[K] = DZR[K] * (F[K]*W[K+1] - F[K-1]*W[K])

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
def CALCDEPTH(Z,ZZ,DZ,DZZ,KB,KL1,KL2):

    BB = float(KL2-KL1) + 4.
    CC = float(KL1) - 2.
    DEL1 = 2./BB/np.exp(.693147*float(KL1-2))
    DEL2 = 2./BB/np.exp(.693147*float(KB-KL2-1))
    Z[0] = 0.
    ZZ[0] = -DEL1/2.

    for K in range(1,int(KL1) - 1):
        Z[K] = -DEL1 * np.exp(.693147 * float(K - 2))
        ZZ[K] = -DEL1 * np.exp(.693147 * (float(K) - 1.5))

    for K in range(int(KL1) - 1,int(KL2) + 2):
        Z[K] = - (float(K) - CC) / BB
        ZZ[K] = - (float(K) - CC + 0.5) / BB

    for K in range(0,int(KB) - 1):
        DZ[K] = Z[K] - Z[K+1]
        DZZ[K] = ZZ[K] - ZZ[K+1]


    return Z, ZZ, DZ, DZZ


def calculate_vertical_grid_spacing(vertical_layers,surface_log_distribution,bottom_log_distribution):

    # KL1 = surface_log_distribution
    # KL2 = bottom_log_distribution
    # KB = vertical_layers

    vertical_coordinates = np.zeros(vertical_layers,dtype=float)
    vertical_coordinates_staggered = np.zeros(vertical_layers,dtype=float)
    vertical_spacing = np.zeros(vertical_layers, dtype=float)
    vertical_spacing_staggered = np.zeros(vertical_layers, dtype=float)

    surface_logspace_layers = surface_log_distribution - 2.
    bottom_logspace_layers = vertical_layers - bottom_log_distribution - 1.

    initial_spacing = (bottom_log_distribution - surface_log_distribution + 4.) * 2**(-6.3)

    vertical_coordinates_staggered[0] = 0.5 * initial_spacing

    for i in range(1,int(surface_log_distribution)-1):
        vertical_coordinates[i] = -initial_spacing * 2**(i-2)
        vertical_coordinates_staggered[i] = -initial_spacing * 2**(i-1.5)

    for i in range(int(surface_log_distribution)-1,vertical_layers):
        vertical_coordinates[i]  = -(i - surface_logspace_layers) / (bottom_log_distribution
                                                                     - surface_log_distribution + 4.)
        vertical_coordinates_staggered[i] = -(i - surface_logspace_layers + 0.5) / (bottom_log_distribution
                                                                                    - surface_log_distribution + 4.)

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
def DENS(T, S, ZZ, DT, RHO, KB):

    GRAV = 9.806

    for K in range(0, KB - 1):
        TR = T[K]
        SR = S[K]

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   HERE, THE (APPROXIMATE) PRESSURE IS IN UNITS OF BARS
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        P = -GRAV * 1.025 * ZZ[K] * DT * 0.01
        RHOR = 999.842594 + 6.793952E-2 * TR - 9.095290E-3 * TR ** 2 + \
               1.001685E-4 * TR ** 3 - 1.120083E-6 * TR ** 4 + 6.536332E-9 * TR ** 5
        RHOR = RHOR + (0.824493 - 4.0899E-3 * TR + 7.6438E-5 * TR ** 2 -
                             8.2467E-7 * TR ** 3 + 5.3875E-9 * TR ** 4) * SR + \
                  (-5.72466E-3 + 1.0227E-4 * TR - 1.6546E-6 * TR ** 2) * (np.abs(SR)) ** 1.5 + 4.8314E-4 * SR ** 2

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   FOR SHALLOW WATER THE PRESSURE DEPENDENCY CAN BE NEGLECTED
        #   IN WHICH IT SHOULD ALSO BE OMITTED IN PROFQ
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        CR = 1449.1 + .0821 * P + 4.55 * TR - .045 * TR ** 2 + 1.34 * (SR - 35.)

        RHO[K] = (RHOR - 1000.) * 1.E-3

    RHO[KB - 1] = RHO[KB - 2]

    return RHO

def calculate_density_profile(temperature,salinity,vertical_coordinates,timestep,vertical_layers):

    density_profile = np.zeros(vertical_layers,dtype=float)
    gravity = 9.806

    for i in range(0,vertical_layers-1):

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   APPROXIMATE PRESSURE IN UNITS OF BARS
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        pressure = -gravity * 1.025 * vertical_coordinates[i] * timestep * 0.01
        rho = 999.842594 + 6.793952E-2 * temperature[i] - 9.095290E-3 * temperature[i] ** 2 + \
               1.001685E-4 * temperature[i] ** 3 - 1.120083E-6 * temperature[i] ** 4 + 6.536332E-9 * temperature[i] ** 5
        rho = rho + (0.824493 - 4.0899E-3 * temperature[i] + 7.6438E-5 * temperature[i] ** 2 -
                             8.2467E-7 * temperature[i] ** 3 + 5.3875E-9 * temperature[i] ** 4) * salinity[i] + \
                  (-5.72466E-3 + 1.0227E-4 * temperature[i] - 1.6546E-6 * temperature[i] ** 2) * \
                  (np.abs(salinity[i])) ** 1.5 + 4.8314E-4 * salinity[i] ** 2
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   FOR SHALLOW WATER THE PRESSURE DEPENDENCY CAN BE NEGLECTED
        #   IN WHICH IT SHOULD ALSO BE OMITTED IN PROFQ
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        density_profile[i] = (rho - 1000.) * 1.E-3

    density_profile[vertical_layers-1] = density_profile[vertical_layers-2]

    return density_profile




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
def MLDPTH(ZZ,T,KB,ZZMLD):

    ZERO = 1.E-06
    T = T

    for K in range(0,KB - 1):

        if T[0] > T[K] + 0.2:
            break

        ZZMLD[K] = ZZ[K] - (T[K] + 0.2 - T[0]) * (ZZ[K] - ZZ[K + 1]) / (T[K] - T[K + 1] + ZERO)

    return ZZMLD

def calculate_mixed_layer_depth(vertical_coordinates,temperature,vertical_layers):

    mixed_layer_depth = np.zeros(vertical_layers,dtype=float)
    zero = 1.E-06

    for i in range(0,vertical_layers-1):

        if temperature[0] > temperature[i]:
            break

        mixed_layer_depth[i] = vertical_coordinates[i] - (temperature[i] + 0.2 - temperature[0]) * \
                            (vertical_coordinates[i] - vertical_coordinates[i+1]) / (temperature[i]
                                                                                        - temperature[i+1] + zero)

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
    localtime = TIME * SEC_PER_DAY

    # SAVING FREQUENCY IN TIME MARCHING LOOP ITERATIONS
    time_to_save = int(out_delta*SEC_PER_DAY)

    # SUMMING UP THE FIELDS TO BE SAVED
    calcmean_bfm(ACCUMULATE)

    # WRITE OUTPUT
    if TT + (DTI/SEC_PER_DAY) > out_delta:

        calcmean_bfm(MEAN)
        save_bfm(localtime)

        # RESET TIME COUNTER
        TT = TT - int(TT)

    # WHEN RUNNING IN COUPLING WITH POM RESTART IS WRITTEN IN MAIN AFTER THE END
    # OF THE TIME MARCHING LOOP

    # WRITE RESTART
    if kt >= IEND:

        if -TT < DTI/SEC_PER_DAY:
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

    # TWICE THE TIME STEP
    # DTI2 = float()
    # The input general cir. vertical vel. is suppose to be in m/s
    W_ON = 1.0
    # The input eddy vertical vel. is provided in m/d
    Weddy_ON = 0.1/86400.0  # to m/s

    DTI2 = DTI * 2.

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL VARIABLES
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    trelax_o2o = NRT_o2o / SEC_PER_DAY
    trelax_n1p = NRT_n1p / SEC_PER_DAY
    trelax_n3n = NRT_n3n / SEC_PER_DAY
    trelax_n4n = NRT_n4n

    # LOOP OVER BFM STATE VAR'S
    for M in range(0,NO_D3_BOX_STATES):

        # ZEROING

        surflux = ZERO
        botflux = ZERO
        fbio[:] = ZERO
        fbbio[:] = ZERO
        ffbio[:] = ZERO
        sink[:] = ZERO
        POCsink = ZERO

        # LOAD BFM STATE VAR.
        for K in range(0,KB-1):
            fbio[K] = D3STATE[M][K]
            fbbio[K] = D3STATEB[M][K]

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   NUTRIENTS SURFACE AND BOTTOM FLUXES
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        if M == ppO2o:
            surflux = -(jsurO2o[0] / SEC_PER_DAY)
            botflux = (o2o[KB-2] - O2BOTT) * trelax_o2o
        elif M == ppO3c:
            surflux = ZERO
        elif M == ppN1p:
            surflux = ZERO
            botflux = (n1p[KB-2] - PO4BOTT) * trelax_n1p
        elif M == ppN3n:
            surflux = ZERO
            botflux = (n3n[KB-2] - NO3BOTT) * trelax_n3n
        elif M == ppN5s:
            surflux = ZERO
        else:
            surflux = ZERO

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

            for K in range(0,KB-1):
                sink[K] = sink[K] - sediR6[K]/SEC_PER_DAY

            # FINAL SINK VALUE
            sink[KB-1] = sink[KB-2]

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

            for K in range(0,KB-1):
                sink[K] = sink[K] - sediPPY[K]/SEC_PER_DAY

            # FINAL SINK VALUE
            sink[KB - 1] = sink[KB - 2]

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   BOTTOM FLUX
        #   Z: Zooplankton
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        # The bot flux for Zooplankton is left equal to ZERO

        # SINKING: UPSTREAM VERTICAL ADVECTION
        adverte(fbbio,fbio,ffbio,sink)

        # SOURCE SPLITTING LEAPFROG INTEGRATION
        for K in range(0,KB-1):
            ffbio[K] = ffbio[K] + DTI2*((ffbio[K]/H) + D3SOURCE[M][K])

        # COMPUTE VERTICAL DIFFUSION AND TERMINATE INTEGRATION
        # IMPLICIT LEAPFROGGING
        PROFTS(ffbio,surflux,botflux,ZERO,ZERO,NBCBFM,DTI2,NTP,UMOLBFM)

        # CLIPPING......IF NEEDED
        for K in range(0,KB-1):
            ffbio[K] = max(p_small,ffbio[K])

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   Mix the time step and restore time sequence
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        for N in range(0,KB-1):
            D3STATEB[M][N] = fbio[N] + Decimal(0.5)*smoth*(ffbio[N] + fbbio[N] - Decimal(2.0)*fbio[N])
            D3STATE[M][N] = ffbio[N]

    if AssignAirPelFluxesInBFMFlag is False:
        jsurO2o[:] = ZERO
        jsurO3c[:] = ZERO

    return


# EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
