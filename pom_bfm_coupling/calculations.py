import numpy as np
from include import BFM_POM
from pom.constants import vertical_layers, twice_the_timestep, seconds_per_day
from inputs import params_POMBFM
from pom_bfm_coupling.data_classes import BfmStateVariableData
from pom.create_profiles import calculate_vertical_temperature_and_salinity_profiles
from bfm.global_parameters import AssignAirPelFluxesInBFMFlag, p_small
from bfm.constants import num_d3_box_states

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
def calculate_vertical_diffusivity(vertical_grid, diffusion, nutrients, d3state, d3stateb, bfm_rates, bfm_phys_vars, dOdt_wind, do3cdt_air_sea_flux):

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   GLOBAL DEFINITION OF PELAGIC (D3/D2) STATE VARIABLES (From ModuleMem)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Values correspond to index from bfm.variable_info.py
    ppo2o = 0; ppn1p = 1; ppn3n = 2
    ppn4n = 3; ppo4n = 4; ppn5s = 5; ppn6r = 6; ppb1c = 7; ppb1n = 8; ppb1p = 9; ppp1c = 10
    ppp1n = 11; ppp1p = 12; ppp1l = 13; ppp1s = 14; ppp2c = 15; ppp2n = 16; ppp2p = 17
    ppp2l = 18; ppP2s = 0; ppp3c = 19; ppp3n = 20; ppp3p = 21; ppp3l = 22; ppP3s = 0
    ppp4c = 23; ppp4n = 24; ppp4p = 25; ppp4l = 26; ppP4s = 0; ppz3c = 27; ppz3n = 28
    ppz3p = 29; ppz4c = 30; ppz4n = 31; ppz4p = 32; ppz5c = 33; ppz5n = 34; ppz5p = 35
    ppz6c = 36; ppz6n = 37; ppz6p = 38; ppr1c = 39; ppr1n = 40; ppr1p = 41; ppR1s = 0
    ppr2c = 42; ppR2n = 0; ppR2p = 0; ppR2s = 0; ppr3c = 43; ppR3n = 0; ppR3p = 0; ppR3s = 0
    ppr6c = 44; ppr6n = 45; ppr6p = 46; ppr6s = 47; ppo3c = 48; ppo3h = 49

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
            bfm_state_var.current[i] = d3state[i,M]
            bfm_state_var.backward[i] = d3stateb[i,M]

        bfm_state_var.current[vertical_layers-1] = bfm_state_var.current[vertical_layers-2]
        bfm_state_var.backward[vertical_layers-1] = bfm_state_var.backward[vertical_layers-2]

        for i in range(0,vertical_layers):
            sinking_velocity[i] = W_ON*bfm_phys_vars.wgen[i] + Weddy_ON*bfm_phys_vars.weddy[i]
        
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   NUTRIENTS SURFACE AND BOTTOM FLUXES
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        if M == ppo2o:   # Dissolved Oxygen (o2o)
            bfm_state_var.surface_flux = -(dOdt_wind[0] / seconds_per_day)
            bfm_state_var.bottom_flux = (d3state[vertical_layers-2,0] - nutrients.O2bott) * trelax_o2o
        elif M == ppo3c:    # Dissolved Inorganic Carbon (o3c)
            bfm_state_var.surface_flux = 0.
        elif M == ppn1p:  # Phosphate (n1p)
            bfm_state_var.surface_flux = 0.
            bfm_state_var.bottom_flux = (d3state[vertical_layers-2,1] - nutrients.PO4bott) * trelax_n1p
        elif M == ppn3n:   # Nitrate (n3n)
            bfm_state_var.surface_flux = 0.
            bfm_state_var.bottom_flux = (d3state[vertical_layers-2,2] - nutrients.NO3bott) * trelax_n3n
        elif M == ppn5s: # Silicate (n5s)
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
        detritus_sedimentation_rate = detritus_sedimentation()
        # The botflux for Particulate Organic Matter is left equal to ZERO
        if ppr6c <= M <= ppr6s:

            for i in range(0,vertical_layers-1):
                sinking_velocity[i] = sinking_velocity[i] - detritus_sedimentation_rate[i]/seconds_per_day

            # FINAL SINK VALUE
            sinking_velocity[vertical_layers-1] = sinking_velocity[vertical_layers-2]

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   SEDIMENTATION PHYTOPLANKTON
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        phyto_sedimentation_rates = phyto_sedimentation()
        # The botflux for Phytoplankton is left equal to ZERO
        if ppp1c <= M <= ppp4l:

            # FROM MODULEMEM --> iiP1 = 1, iiP2 = 2, 11P3 = 3, iiP4 = 4
            if M in range(ppp1c,ppp1s):
                K = 1   # iiP1
            elif M in range(ppp2c,ppp2l):
                K = 2   # iiP2
            elif M in range(ppp3c,ppp3l):
                K = 3   # iiP3
            elif M in range(ppp4c,ppp4l):
                K = 4   # iiP4

            for i in range(0,vertical_layers-1):
                sinking_velocity[i] = sinking_velocity[i] - phyto_sedimentation_rates[i,K-1]/seconds_per_day

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
            bfm_state_var.forward[i] = bfm_state_var.backward[i] + twice_the_timestep*((bfm_state_var.forward[i]/params_POMBFM.h) + bfm_rates[i,M])
        # bfm_state_var.forward[vertical_layers-1] = bfm_state_var.forward[vertical_layers-2]

        # COMPUTE VERTICAL DIFFUSION AND TERMINATE INTEGRATION
        # IMPLICIT LEAPFROGGING
        bfm_state_var = calculate_vertical_temperature_and_salinity_profiles(vertical_grid, diffusion, bfm_state_var, 0, params_POMBFM.nbcbfm, params_POMBFM.umolbfm)

        # CLIPPING......IF NEEDED
        for i in range(0,vertical_layers-1):
            bfm_state_var.forward[i] = max(p_small,bfm_state_var.forward[i])

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   Mix the time step and restore time sequence
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        for N in range(0,vertical_layers-1):
            d3stateb[N,M] = bfm_state_var.current[N] + 0.5*params_POMBFM.smoth*(bfm_state_var.forward[N] + bfm_state_var.backward[N] - 2.*bfm_state_var.current[N])
            d3state[N,M] = bfm_state_var.forward[N]

    if not AssignAirPelFluxesInBFMFlag:
        dOdt_wind[:] = 0.
        do3cdt_air_sea_flux[:] = 0.

    return d3state, d3stateb


def detritus_sedimentation():

    # FROM namelists_bfm
    # p_rR6m        [m/d]   detritus sinking rate
    # p_burvel_R6   [m/d]   Bottom burial velocity for detritus

    p_rR6m = 1.
    p_burvel_R6 = 1.
    # p_burvel_R6 = 1.5

    # FROM PelGlobal.F90 (145-148)
    detritus_sedimentation_rate = p_rR6m * np.ones(vertical_layers-1)
    # try:
    #     BFM_POM
    # except NameError:
    #     BFM_POM = False
    if not BFM_POM:
        detritus_sedimentation_rate[vertical_layers-2] = p_burvel_R6

    return detritus_sedimentation_rate


def phyto_sedimentation():

    # FROM namelists_bfm
    # p_rPIm        [m/d]   phytoplanktom background sinking rate
    # p_burvel_PI   [m/d]   Botttom burial velocity for detritus
    p_rPIm = [0.0, 0.0, 0.0, 0.0]
    p_burvel_PI = 0.0

    # FROM MODULEMEM.F90 (338)
    iiPhytoPlankton = 4
    iiP1 = 1; iiP2 = 2; iiP3 = 3; iiP4 = 4

    # FROM PelGLobal.F90 (149-154)
    phyto_sedimentation_rates = np.zeros((vertical_layers-1,iiPhytoPlankton))
    for i in range(0,iiPhytoPlankton):
        phyto_sedimentation_rates[:,i] = p_rPIm[i]
        # try:
        #     BFM_POM
        # except NameError:
        #     BFM_POM = False
        if not BFM_POM:
            phyto_sedimentation_rates[vertical_layers-2,i] = p_burvel_PI

    return phyto_sedimentation_rates

