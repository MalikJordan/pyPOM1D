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
    ppdisOxygen_IO_O = 0; ppphospate_IO_P = 1; ppnitrate_IO_N = 2
    ppammonium_IO_N = 3; ppnitrogenSink = 4; ppsilicate_IO_Si = 5; ppreductEquiv_IO_R = 6; pppelBacteria_LO_C = 7; pppelBacteria_LO_N = 8; pppelBacteria_LO_P = 9; ppdiatoms_LO_C = 10
    ppdiatoms_LO_N = 11; ppdiatoms_LO_P = 12; ppdiatoms_LO_Chl = 13; ppdiatoms_LO_Si = 14; ppnanoflagellates_LO_C = 15; ppnanoflagellates_LO_N = 16; ppnanoflagellates_LO_P = 17
    ppnanoflagellates_LO_Chl = 18; ppP2s = 0; pppicophyto_LO_C = 19; pppicophyto_LO_N = 20; pppicophyto_LO_P = 21; pppicophyto_LO_Chl = 22; ppP3s = 0
    pplargephyto_LO_C = 23; pplargephyto_LO_N = 24; pplargephyto_LO_P = 25; pplargephyto_LO_Chl = 26; ppP4s = 0; ppcarnivMesozoo_LO_C = 27; ppcarnivMesozoo_LO_N = 28
    ppcarnivMesozoo_LO_P = 29; ppomnivMesozoo_LO_C = 30; ppomnivMesozoo_LO_N = 31; ppomnivMesozoo_LO_P = 32; ppmicrozoo_LO_C = 33; ppmicrozoo_LO_N = 34; ppmicrozoo_LO_P = 35
    ppheteroFlagellates_LO_C = 36; ppheteroFlagellates_LO_N = 37; ppheteroFlagellates_LO_P = 38; pplabileDOM_NO_C = 39; pplabileDOM_NO_N = 40; pplabileDOM_NO_P = 41; ppR1s = 0
    ppsemilabileDOC_NO_C = 42; ppR2n = 0; ppR2p = 0; ppR2s = 0; ppsemirefractDOC_NO_C = 43; ppR3n = 0; ppR3p = 0; ppR3s = 0
    ppparticOrganDetritus_NO_C = 44; ppparticOrganDetritus_NO_N = 45; ppparticOrganDetritus_NO_P = 46; ppparticOrganDetritus_NO_Si = 47; ppdisInorgCarbon_IO_C = 48; pptotalAlkalinity_IO = 49

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
            bfm_state_var.current[i] = d3state[i,M]
            bfm_state_var.backward[i] = d3stateb[i,M]

        bfm_state_var.current[vertical_layers-1] = bfm_state_var.current[vertical_layers-2]
        bfm_state_var.backward[vertical_layers-1] = bfm_state_var.backward[vertical_layers-2]

        for i in range(0,vertical_layers):
            sinking_velocity[i] = W_ON*bfm_phys_vars.wgen[i] + Weddy_ON*bfm_phys_vars.weddy[i]
        
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   NUTRIENTS SURFACE AND BOTTOM FLUXES
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        if M == ppdisOxygen_IO_O:   # Dissolved Oxygen (o2o)
            bfm_state_var.surface_flux = -(dOdt_wind[0] / seconds_per_day)
            bfm_state_var.bottom_flux = (d3state[vertical_layers-2,0] - nutrients.O2bott) * trelax_disOxygen_IO_O
        elif M == ppdisInorgCarbon_IO_C:    # Dissolved Inorganic Carbon (o3c)
            bfm_state_var.surface_flux = 0.
        elif M == ppphospate_IO_P:  # Phosphate (n1p)
            bfm_state_var.surface_flux = 0.
            bfm_state_var.bottom_flux = (d3state[vertical_layers-2,1] - nutrients.PO4bott) * trelax_phospate_IO_P
        elif M == ppnitrate_IO_N:   # Nitrate (n3n)
            bfm_state_var.surface_flux = 0.
            bfm_state_var.bottom_flux = (d3state[vertical_layers-2,2] - nutrients.NO3bott) * trelax_nitrate_IO_N
        elif M == ppsilicate_IO_Si: # Silicate (n5s)
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
        if ppparticOrganDetritus_NO_C <= M <= ppparticOrganDetritus_NO_Si:

            for i in range(0,vertical_layers-1):
                sinking_velocity[i] = sinking_velocity[i] - detritus_sedimentation_rate[i]/seconds_per_day

            # FINAL SINK VALUE
            sinking_velocity[vertical_layers-1] = sinking_velocity[vertical_layers-2]

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   SEDIMENTATION PHYTOPLANKTON
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        phyto_sedimentation_rates = phyto_sedimentation()
        # The botflux for Phytoplankton is left equal to ZERO
        if ppdiatoms_LO_C <= M <= pplargephyto_LO_Chl:

            # FROM MODULEMEM --> iiP1 = 1, iiP2 = 2, 11P3 = 3, iiP4 = 4
            if M in range(ppdiatoms_LO_C,ppdiatoms_LO_Si):
                K = 1   # iiP1
            elif M in range(ppnanoflagellates_LO_C,ppnanoflagellates_LO_Chl):
                K = 2   # iiP2
            elif M in range(pppicophyto_LO_C,pppicophyto_LO_Chl):
                K = 3   # iiP3
            elif M in range(pplargephyto_LO_C,pplargephyto_LO_Chl):
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

