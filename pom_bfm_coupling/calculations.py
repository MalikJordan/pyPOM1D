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
def calculate_vertical_diffusivity(vertical_grid, diffusion, nutrients, o2o, n1p, n3n):

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   GLOBAL DEFINITION OF PELAGIC (D3/D2) STATE VARIABLES (From ModuleMem)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Values correspond to index from bfm.variable_info.py
    ppO2o = 0; ppN1p = 1; ppN3n = 2
    ppN4n = 3; ppO4n = 4; ppN5s = 5; ppN6r = 6; ppB1c = 7; ppB1n = 8; ppB1p = 9; ppP1c = 10
    ppP1n = 11; ppP1p = 12; ppP1l = 13; ppP1s = 14; ppP2c = 15; ppP2n = 16; ppP2p = 17
    ppP2l = 18; ppP2s = 0; ppP3c = 19; ppP3n = 20; ppP3p = 21; ppP3l = 22; ppP3s = 0
    ppP4c = 23; ppP4n = 24; ppP4p = 25; ppP4l = 26; ppP4s = 0; ppZ3c = 27; ppZ3n = 28
    ppZ3p = 29; ppZ4c = 30; ppZ4n = 31; ppZ4p = 32; ppZ5c = 33; ppZ5n = 34; ppZ5p = 35
    ppZ6c = 36; ppZ6n = 37; ppZ6p = 38; ppR1c = 39; ppR1n = 40; ppR1p = 41; ppR1s = 0
    ppR2c = 42; ppR2n = 0; ppR2p = 0; ppR2s = 0; ppR3c = 43; ppR3n = 0; ppR3p = 0; ppR3s = 0
    ppR6c = 44; ppR6n = 45; ppR6p = 46; ppR6s = 47; ppO3c = 48; ppO3h = 49

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

    trelax_o2o = params_POMBFM.nrt_o2o / seconds_per_day
    trelax_n1p = params_POMBFM.nrt_n1p / seconds_per_day
    trelax_n3n = params_POMBFM.nrt_n3n / seconds_per_day
    trelax_n4n = params_POMBFM.nrt_n4n

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

        if M == ppO2o:
            bfm_state_var.surface_flux = -(jsurO2o[0] / seconds_per_day)
            bfm_state_var.bottom_flux = (o2o[vertical_layers-2] - nutrients.O2BOTT) * trelax_o2o
        elif M == ppO3c:
            bfm_state_var.surface_flux = 0.
        elif M == ppN1p:
            bfm_state_var.surface_flux = 0.
            bfm_state_var.bottom_flux = (n1p[vertical_layers-2] - nutrients.PO4BOTT) * trelax_n1p
        elif M == ppN3n:
            bfm_state_var.surface_flux = 0.
            bfm_state_var.bottom_flux = (n3n[vertical_layers-2] - nutrients.NO3BOTT) * trelax_n3n
        elif M == ppN5s:
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
        if ppR6c <= M <= ppR6s:

            for i in range(0,vertical_layers-1):
                sinking_velocity[i] = sinking_velocity[i] - sediR6[i]/seconds_per_day

            # FINAL SINK VALUE
            sinking_velocity[vertical_layers-1] = sinking_velocity[vertical_layers-2]

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   SEDIMENTATION PHYTOPLANKTON
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        # The botflux for Phytoplankton is left equal to ZERO
        if ppP1c <= M <= ppP4l:

            # FROM MODULEMEM --> iiP1 = 1, iiP2 = 2, 11P3 = 3, iiP4 = 4
            if M in range(ppP1c,ppP1s):
                N = 1   # iiP1
            elif M in range(ppP2c,ppP2l):
                N = 2   # iiP2
            elif M in range(ppP3c,ppP3l):
                N = 3   # iiP3
            elif M in range(ppP4c,ppP4l):
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
        calculate_vertical_advection(bfm_state_var,sinking_velocity,vertical_grid)

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
        jsurO2o[:] = 0.
        jsurO3c[:] = 0.


    return




