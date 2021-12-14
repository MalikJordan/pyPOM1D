# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from cppdefs import *
import numpy as np
from inputs import params_POMBFM
from pom.constants import twice_the_timestep, vertical_layers, seconds_per_day
from pom.data_classes import VerticalGridData


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
def create_vertical_coordinate_system(surface_layers_with_log_distribution, bottom_layers_with_log_distribution):

    # KL1 = surface_log_distribution
    # KL2 = bottom_log_distribution
    # KB = vertical_layers

    vertical_coordinates = np.zeros(vertical_layers)
    vertical_coordinates_staggered = np.zeros(vertical_layers)
    vertical_spacing = np.zeros(vertical_layers)
    vertical_spacing_staggered = np.zeros(vertical_layers)
    vertical_spacing_reciprocal = np.zeros(vertical_layers)

    surface_logspace_layers = surface_layers_with_log_distribution - 2.
    bottom_logspace_layers = vertical_layers - bottom_layers_with_log_distribution - 1.

    BB = (bottom_layers_with_log_distribution - surface_layers_with_log_distribution) + 4.
    CC = surface_layers_with_log_distribution - 2.
    initial_spacing = 2. / BB / np.exp(.693147 * (surface_layers_with_log_distribution - 2))

    vertical_coordinates_staggered[0] = -0.5 * initial_spacing

    for i in range(1, int(surface_layers_with_log_distribution) - 1):
        vertical_coordinates[i-1] = -initial_spacing * 2**(i-2)
        vertical_coordinates_staggered[i-1] = -initial_spacing * 2**(i-1.5)

    for i in range(int(surface_layers_with_log_distribution) - 1, vertical_layers+1):
        vertical_coordinates[i-1]  = -(i - surface_logspace_layers) / (bottom_layers_with_log_distribution
                                                                       - surface_layers_with_log_distribution + 4.)
        vertical_coordinates_staggered[i-1] = -(i - surface_logspace_layers + 0.5) / (bottom_layers_with_log_distribution
                                                                                      - surface_layers_with_log_distribution + 4.)

    for i in range(0,vertical_layers-1):
        vertical_spacing[i] = vertical_coordinates[i] - vertical_coordinates[i+1]
        vertical_spacing_staggered[i] = vertical_coordinates_staggered[i] - vertical_coordinates_staggered[i+1]

    # VERTICAL SPACING RECIROCAL (DZR)
    vertical_spacing[vertical_layers-1] = 1.E-6  # Small value to avoid division by zero for vertical_spacing_reciprocal
    vertical_spacing_reciprocal[:] = 1. / vertical_spacing[:]
    vertical_spacing[vertical_layers-1] = 0.

    vertical_grid = VerticalGridData()
    vertical_grid.vertical_coordinates = vertical_coordinates
    vertical_grid.vertical_coordinates_staggered = vertical_coordinates_staggered
    vertical_grid.vertical_spacing = vertical_spacing
    vertical_grid.vertical_spacing_staggered = vertical_spacing_staggered
    vertical_grid.vertical_spacing_reciprocal = vertical_spacing_reciprocal

    return vertical_grid



# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# ROUTINE: CalcLightDistribution
#
# DESCRIPTION
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# def CalcLightDistribution():
#
#     EIR[0] = EIR[0] * p_PAR / E2W
#     for BoxNumberZ in range(1, NO_BOXES_Z):
#         box_no = BoxNumberZ
#         EIR[box_no] = EIR[box_no - 1] * np.exp(- 1.0 * xEPS[box_no - 1] * Depth[box_no - 1])
#
#     return EIR


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
def calculate_vertical_density_profile(temperature,salinity,vertical_grid):

    vertical_density_profile = np.zeros(vertical_layers)
    gravity = 9.806
    for i in range(0,vertical_layers-1):

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   APPROXIMATE PRESSURE IN UNITS OF BARS
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        pressure = -gravity * 1.025 * vertical_grid.vertical_spacing_staggered[i] * params_POMBFM.dti * 0.01
        density = 999.842594 + 6.793952E-2 * temperature.current[i] - 9.095290E-3 * temperature.current[i] ** 2 + \
               1.001685E-4 * temperature.current[i] ** 3 - 1.120083E-6 * temperature.current[i] ** 4 + 6.536332E-9 * temperature.current[i] ** 5
        density = density + (0.824493 - 4.0899E-3 * temperature.current[i] + 7.6438E-5 * temperature.current[i] ** 2 -
                             8.2467E-7 * temperature.current[i] ** 3 + 5.3875E-9 * temperature.current[i] ** 4) * salinity.current[i] + \
                  (-5.72466E-3 + 1.0227E-4 * temperature.current[i] - 1.6546E-6 * temperature.current[i] ** 2) * \
                  (np.abs(salinity.current[i])) ** 1.5 + 4.8314E-4 * salinity.current[i] ** 2
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
def calculate_mixed_layer_depth(vertical_coordinates_staggered, temperature):

    mixed_layer_depth = np.zeros(vertical_layers)
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
# def pom_dia_bfm(kt,TT):
#
#     # need netcdf_bfm stuff
#
#     # TIME IN SECONDS
#     localtime = float()
#
#     # SAVING FREQUENCY
#     time_to_save = int()  # time in seconds
#
#     # TIME ELLAPSED (IN SECONDS) SINCE THE BEGINNING OF THE SIMULATION
#     localtime = TIME * seconds_per_day
#
#     # SAVING FREQUENCY IN TIME MARCHING LOOP ITERATIONS
#     time_to_save = int(out_delta*seconds_per_day)
#
#     # SUMMING UP THE FIELDS TO BE SAVED
#     calcmean_bfm(ACCUMULATE)
#
#     # WRITE OUTPUT
#     if TT + (dti/seconds_per_day) > out_delta:
#
#         calcmean_bfm(MEAN)
#         save_bfm(localtime)
#
#         # RESET TIME COUNTER
#         TT = TT - int(TT)
#
#     # WHEN RUNNING IN COUPLING WITH POM RESTART IS WRITTEN IN MAIN AFTER THE END
#     # OF THE TIME MARCHING LOOP
#
#     # WRITE RESTART
#     if kt >= IEND:
#
#         if -TT < dti/seconds_per_day:
#             save_rst_bfm(localtime)
#             close_ncdf(ncid_rst)
#             close_ncdf(ncid_bfm)
#             print('POM_DIA: NETCDF RESTART WRITTEN, TIME--> ', TIME, kt, IEND, DTI, TT)
#
#     return


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
# def vdiff_SOS():
#
#     # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#     #   GLOBAL DEFINITION OF PELAGIC (D3/D2) STATE VARIABLES (From ModuleMem)
#     # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#     ppdisOxygen_IO_O = 1; ppphospate_IO_P = 2; ppnitrate_IO_N = 3
#     ppammonium_IO_N = 4; ppO4n = 5; ppsilicate_IO_Si = 6; ppreductEquiv_IO_R = 7; pppelBacteria_LO_C = 8; pppelBacteria_LO_N = 9; pppelBacteria_LO_P = 10; ppdiatoms_LO_C = 11
#     ppdiatoms_LO_N = 12; ppdiatoms_LO_P = 13; ppdiatoms_LO_Chl = 14; ppdiatoms_LO_Si = 15; ppnanoflagellates_LO_C = 16; ppnanoflagellates_LO_N = 17; ppnanoflagellates_LO_P = 18
#     ppnanoflagellates_LO_Chl = 19; ppP2s = 0; pppicophyto_LO_C = 20; pppicophyto_LO_N = 21; pppicophyto_LO_P = 22; pppicophyto_LO_Chl = 23; ppP3s = 0
#     pplargephyto_LO_C = 24; pplargephyto_LO_N = 25; pplargephyto_LO_P = 26; pplargephyto_LO_Chl = 27; ppP4s = 0; ppcarnivMesozoo_LO_C = 28; ppcarnivMesozoo_LO_N = 29
#     ppcarnivMesozoo_LO_P = 30; ppomnivMesozoo_LO_C = 31; ppomnivMesozoo_LO_N = 32; ppomnivMesozoo_LO_P = 33; ppmicrozoo_LO_C = 34; ppmicrozoo_LO_N = 35; ppmicrozoo_LO_P = 36
#     ppheteroFlagellates_LO_C = 37; ppheteroFlagellates_LO_N = 38; ppheteroFlagellates_LO_P = 39; pplabileDOM_NO_C = 40; pplabileDOM_NO_N = 41; pplabileDOM_NO_P = 42; ppR1s = 0
#     ppsemilabileDOC_NO_C = 43; ppR2n = 0; ppR2p = 0; ppR2s = 0; ppsemirefractDOC_NO_C = 44; ppR3n = 0; ppR3p = 0; ppR3s = 0
#     ppparticOrganDetritus_NO_C = 45; ppparticOrganDetritus_NO_N = 46; ppparticOrganDetritus_NO_P = 47; ppparticOrganDetritus_NO_Si = 48; ppdisInorgCarbon_IO_C = 49; pptotalAlkalinity_IO = 50
#
#     # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#     #   LOCAL VARIABLES
#     # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
#     # COUNTER & FLAGS
#     # K, M, N, NBC = float()
#
#     # BFM STATE VAR. @ time t-DTI, t, T+DTI RESPECTIVELY
#     fbio, ffbio, fbbio = np.zeros(vertical_layers)
#
#     # SURFACE FLUX STORAGE FOR NUT'S O2 & CO2
#     # surflux = float()
#     # botflux = float()
#
#     # RELAXATION VELOCITY FOR NUT'S
#     # trelax_disOxygen_IO_O = float()
#     # trelax_phospate_IO_P = float()
#     # trelax_nitrate_IO_N = float()
#     # trelax_ammonium_IO_N = float()
#
#     # SEDIMENTATION VELOCITY
#     sink = np.zeros(vertical_layers)
#     # POCsink = float()
#     # W1R6 = float()
#
#     # The input general cir. vertical vel. is suppose to be in m/s
#     W_ON = 1.0
#     # The input eddy vertical vel. is provided in m/d
#     Weddy_ON = 0.1/86400.0  # to m/s
#
#     # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#     #   LOCAL VARIABLES
#     # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
#     trelax_disOxygen_IO_O = nrt_disOxygen_IO_O / seconds_per_day
#     trelax_phospate_IO_P = nrt_phospate_IO_P / seconds_per_day
#     trelax_nitrate_IO_N = nrt_nitrate_IO_N / seconds_per_day
#     trelax_ammonium_IO_N = nrt_ammonium_IO_N
#
#     # LOOP OVER BFM STATE VAR'S
#     for M in range(0,no_d3_box_states):
#
#         # ZEROING
#
#         surflux = zero
#         botflux = zero
#         fbio[:] = zero
#         fbbio[:] = zero
#         ffbio[:] = zero
#         sink[:] = zero
#         POCsink = zero
#
#         # LOAD BFM STATE VAR.
#         for i in range(0,vertical_layers-1):
#             fbio[i] = D3STATE[M][i]
#             fbbio[i] = D3STATEB[M][i]
#
#         # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#         #   NUTRIENTS SURFACE AND BOTTOM FLUXES
#         # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
#         if M == ppdisOxygen_IO_O:
#             surflux = -(jsurdisOxygen_IO_O[0] / seconds_per_day)
#             botflux = (disOxygen_IO_O[vertical_layers-2] - O2BOTT) * trelax_disOxygen_IO_O
#         elif M == ppdisInorgCarbon_IO_C:
#             surflux = zero
#         elif M == ppphospate_IO_P:
#             surflux = zero
#             botflux = (phospate_IO_P[vertical_layers-2] - PO4BOTT) * trelax_phospate_IO_P
#         elif M == ppnitrate_IO_N:
#             surflux = zero
#             botflux = (nitrate_IO_N[vertical_layers-2] - NO3BOTT) * trelax_nitrate_IO_N
#         elif M == ppsilicate_IO_Si:
#             surflux = zero
#         else:
#             surflux = zero
#
#         # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#         #   BOTTOM FLUX
#         #   R1: Dissolved Organic Matter
#         # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
#         # The botflux for Dissolved Organic Matter is left equal to ZERO
#
#         # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#         #   BOTTOM FLUX
#         #   R6: Particulate Organic Matter
#         # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
#         # The botflux for Particulate Organic Matter is left equal to ZERO
#
#         if M >= ppparticOrganDetritus_NO_C and M <= ppparticOrganDetritus_NO_Si:
#
#             for i in range(0,vertical_layers-1):
#                 sink[i] = sink[i] - sediR6[i]/seconds_per_day
#
#             # FINAL SINK VALUE
#             sink[vertical_layers-1] = sink[vertical_layers-2]
#
#         # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#         #   SEDIMENTATION PHYTOPLANKTON
#         # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
#         # The botflux for Phytoplankton is left equal to ZERO
#
#         if M >= ppdiatoms_LO_C and M <= pplargephyto_LO_Chl:
#
#             if M in range(ppdiatoms_LO_C,ppdiatoms_LO_Si):
#                 N = iiP1
#             elif M in range(ppnanoflagellates_LO_C,ppnanoflagellates_LO_Chl):
#                 N = iiP2
#             elif M in range(pppicophyto_LO_C,pppicophyto_LO_Chl):
#                 N = iiP3
#             elif M in range(pplargephyto_LO_C,pplargephyto_LO_Chl):
#                 N = iiP4
#
#             for i in range(0,vertical_layers-1):
#                 sink[i] = sink[i] - sediPPY[i]/seconds_per_day
#
#             # FINAL SINK VALUE
#             sink[vertical_layers - 1] = sink[vertical_layers - 2]
#
#         # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#         #   BOTTOM FLUX
#         #   Z: Zooplankton
#         # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
#         # The bot flux for Zooplankton is left equal to ZERO
#
#         # SINKING: UPSTREAM VERTICAL ADVECTION
#         adverte(fbbio,fbio,ffbio,sink)
#
#         # SOURCE SPLITTING LEAPFROG INTEGRATION
#         for i in range(0,vertical_layers-1):
#             ffbio[i] = ffbio[i] + twice_the_timestep*((ffbio[i]/h) + D3SOURCE[M][i])
#
#         # COMPUTE VERTICAL DIFFUSION AND TERMINATE INTEGRATION
#         # IMPLICIT LEAPFROGGING
#         PROFTS(ffbio,surflux,botflux,zero,zero,nbcbfm,twice_the_timestep,ntp,umolbfm)
#
#         # CLIPPING......IF NEEDED
#         for i in range(0,vertical_layers-1):
#             ffbio[i] = max(p_small,ffbio[i])
#
#         # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#         #   Mix the time step and restore time sequence
#         # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
#         for N in range(0,vertical_layers-1):
#             D3STATEB[M][N] = fbio[N] + 0.5*smoth*(ffbio[N] + fbbio[N] - 2.*fbio[N])
#             D3STATE[M][N] = ffbio[N]
#
#     if AssignAirPelFluxesInBFMFlag is False:
#         jsurdisOxygen_IO_O[:] = zero
#         jsurdisInorgCarbon_IO_C[:] = zero
#
#     return


# EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
