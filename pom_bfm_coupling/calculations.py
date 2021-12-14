import numpy as np
from pom.constants import vertical_layers, twice_the_timestep, seconds_per_day
from inputs import params_POMBFM
from pom_bfm_coupling.data_classes import BfmStateVariableData
from pom.create_profiles import calculate_vertical_temperature_and_salinity_profiles
from bfm.global_parameters import AssignAirPelFluxesInBFMFlag, p_small
from bfm.constants import num_d3_box_states
from main_pombfm1d import d3state, d3stateb

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
def calculate_vertical_advection(property, sinking_velocity, vertical_grid):
    # sinking velocity input from vdiff_SOS
    property.current[vertical_layers-1] = property.current[vertical_layers-2]
    property.backward[vertical_layers-1] = property.backward[vertical_layers-2]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Calculate vertical advection. Mind downward velocities are negative
    #   Upwind scheme:
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    property.forward[0] = vertical_grid.vertical_spacing_reciprocal[0] * property.current[0] * sinking_velocity[1]

    for i in range(1,vertical_layers-1):
        property.forward[i] = vertical_grid.vertical_spacing_reciprocal[i] * (property.current[i] * sinking_velocity[i + 1] - property.current[i - 1] * sinking_velocity[i])

    return property


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
def calculate_vertical_diffusivity(vertical_grid, diffusion, nutrients, disOxygen_IO_O, phospate_IO_P, nitrate_IO_N):

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   GLOBAL DEFINITION OF PELAGIC (D3/D2) STATE VARIABLES (From ModuleMem)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Values correspond to index from bfm.variable_info.py
    ppdisOxygen_IO_O = 0; ppphospate_IO_P = 1; ppnitrate_IO_N = 2
    ppammonium_IO_N = 3; ppO4n = 4; ppsilicate_IO_Si = 5; ppreductEquiv_IO_R = 6; pppelBacteria_LO_C = 7; pppelBacteria_LO_N = 8; pppelBacteria_LO_P = 9; ppdiatoms_LO_C = 10
    ppdiatoms_LO_N = 11; ppdiatoms_LO_P = 12; ppdiatoms_LO_Chl = 13; ppdiatoms_LO_Si = 14; ppnanoflagellates_LO_C = 15; ppnanoflagellates_LO_N = 16; ppnanoflagellates_LO_P = 17
    ppnanoflagellates_LO_Chl = 18; ppP2s = 0; pppicophyto_LO_C = 19; pppicophyto_LO_N = 20; pppicophyto_LO_P = 21; pppicophyto_LO_Chl = 22; ppP3s = 0
    pplargephyto_LO_C = 23; pplargephyto_LO_N = 24; pplargephyto_LO_P = 25; pplargephyto_LO_Chl = 26; ppP4s = 0; ppcarnivMesozoo_LO_C = 27; ppcarnivMesozoo_LO_N = 28
    ppcarnivMesozoo_LO_P = 29; ppomnivMesozoo_LO_C = 30; ppomnivMesozoo_LO_N = 31; ppomnivMesozoo_LO_P = 32; ppmicrozoo_LO_C = 33; ppmicrozoo_LO_N = 34; ppmicrozoo_LO_P = 35
    ppZ6c = 36; ppZ6n = 37; ppZ6p = 38; pplabileDOM_NO_C = 39; pplabileDOM_NO_N = 40; pplabileDOM_NO_P = 41; ppR1s = 0
    ppsemilabileDOC_NO_C = 42; ppR2n = 0; ppR2p = 0; ppR2s = 0; ppsemirefractDOC_NO_C = 43; ppR3n = 0; ppR3p = 0; ppR3s = 0
    ppparticOrganDetritus_NO_C = 44; ppparticOrganDetritus_NO_N = 45; ppparticOrganDetritus_NO_P = 46; ppparticOrganDetritus_NO_Si = 47; ppdisInorgCarbon_IO_C = 48; ppO3h = 49

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL VARIABLES
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # BFM STATE VAR. @ time t-DTI, t, T+DTI RESPECTIVELY
    bfm_state_var = BfmStateVariableData()

    # SEDIMENTATION VELOCITY
    sinking_velocity = np.zeros(vertical_layers)

    # The input general cir. vertical vel. is suppose to be in m/s
    W_ON = 1.0
    # The input eddy vertical vel. is provided in m/d
    Weddy_ON = 0.1/86400.0  # to m/s

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL VARIABLES
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    trelax_disOxygen_IO_O = params_POMBFM.nrt_disOxygen_IO_O / seconds_per_day
    trelax_phospate_IO_P = params_POMBFM.nrt_phospate_IO_P / seconds_per_day
    trelax_nitrate_IO_N = params_POMBFM.nrt_nitrate_IO_N / seconds_per_day
    trelax_ammonium_IO_N = params_POMBFM.nrt_ammonium_IO_N

    # LOOP OVER BFM STATE VAR'S
    for M in range(0,num_d3_box_states):

        # ZEROING
        bfm_state_var.surface_flux = 0.
        bfm_state_var.bottom_flux = 0.
        bfm_state_var.current[:] = 0.
        bfm_state_var.backward[:] = 0.
        bfm_state_var.forward[:] = 0.
        sinking_velocity[:] = 0.
        POCsink = 0.

        # LOAD BFM STATE VAR.
        for i in range(0,vertical_layers-1):
            bfm_state_var.current[i] = d3state[M][i]
            bfm_state_var.backward[i] = d3stateb[M][i]

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   NUTRIENTS SURFACE AND BOTTOM FLUXES
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        if M == ppdisOxygen_IO_O:
            bfm_state_var.surface_flux = -(jsurdisOxygen_IO_O[0] / seconds_per_day)
            bfm_state_var.bottom_flux = (disOxygen_IO_O[vertical_layers-2] - nutrients.O2BOTT) * trelax_disOxygen_IO_O
        elif M == ppdisInorgCarbon_IO_C:
            bfm_state_var.surface_flux = 0.
        elif M == ppphospate_IO_P:
            bfm_state_var.surface_flux = 0.
            bfm_state_var.bottom_flux = (phospate_IO_P[vertical_layers-2] - nutrients.PO4BOTT) * trelax_phospate_IO_P
        elif M == ppnitrate_IO_N:
            bfm_state_var.surface_flux = 0.
            bfm_state_var.bottom_flux = (nitrate_IO_N[vertical_layers-2] - nutrients.NO3BOTT) * trelax_nitrate_IO_N
        elif M == ppsilicate_IO_Si:
            bfm_state_var.surface_flux = 0.
        else:
            bfm_state_var.surface_flux = 0.

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
        if ppparticOrganDetritus_NO_C <= M <= ppparticOrganDetritus_NO_Si:

            for i in range(0,vertical_layers-1):
                sinking_velocity[i] = sinking_velocity[i] - sediR6[i]/seconds_per_day

            # FINAL SINK VALUE
            sinking_velocity[vertical_layers-1] = sinking_velocity[vertical_layers-2]

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   SEDIMENTATION PHYTOPLANKTON
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        # The botflux for Phytoplankton is left equal to ZERO
        if ppdiatoms_LO_C <= M <= pplargephyto_LO_Chl:

            # FROM MODULEMEM --> iiP1 = 1, iiP2 = 2, 11P3 = 3, iiP4 = 4
            if M in range(ppdiatoms_LO_C,ppdiatoms_LO_Si):
                N = 1   # iiP1
            elif M in range(ppnanoflagellates_LO_C,ppnanoflagellates_LO_Chl):
                N = 2   # iiP2
            elif M in range(pppicophyto_LO_C,pppicophyto_LO_Chl):
                N = 3   # iiP3
            elif M in range(pplargephyto_LO_C,pplargephyto_LO_Chl):
                N = 4   # iiP4

            for i in range(0,vertical_layers-1):
                sinking_velocity[i] = sinking_velocity[i] - sediPPY[i]/seconds_per_day

            # FINAL SINK VALUE
            sinking_velocity[vertical_layers - 1] = sinking_velocity[vertical_layers - 2]

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   BOTTOM FLUX
        #   Z: Zooplankton
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        # The bot flux for Zooplankton is left equal to ZERO

        # SINKING: UPSTREAM VERTICAL ADVECTION
        bfm_state_var = calculate_vertical_advection(bfm_state_var,sinking_velocity,vertical_grid)

        # SOURCE SPLITTING LEAPFROG INTEGRATION
        for i in range(0,vertical_layers-1):
            bfm_state_var.forward[i] = bfm_state_var.forward[i] + twice_the_timestep*((bfm_state_var.forward[i]/params_POMBFM.h) + d3source[M][i])

        # COMPUTE VERTICAL DIFFUSION AND TERMINATE INTEGRATION
        # IMPLICIT LEAPFROGGING
        calculate_vertical_temperature_and_salinity_profiles(vertical_grid, diffusion, bfm_state_var, 0, params_POMBFM.nbcbfm, params_POMBFM.umolbfm)

        # CLIPPING......IF NEEDED
        for i in range(0,vertical_layers-1):
            bfm_state_var.forward[i] = max(p_small,bfm_state_var.forward[i])

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   Mix the time step and restore time sequence
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        for N in range(0,vertical_layers-1):
            d3stateb[M][N] = bfm_state_var.current[N] + 0.5*params_POMBFM.smoth*(bfm_state_var.forward[N] + bfm_state_var.backward[N] - 2.*bfm_state_var.current[N])
            d3state[M][N] = bfm_state_var.forward[N]

    if not AssignAirPelFluxesInBFMFlag:
        jsurdisOxygen_IO_O[:] = 0.
        jsurdisInorgCarbon_IO_C[:] = 0.


    return




