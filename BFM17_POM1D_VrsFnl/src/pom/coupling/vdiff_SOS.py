# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: vdiff_SOS
#
# DESCRIPTION
#
# This routine calculates the vertical diffusivity of BFM biochemical components and
# integrates BFM state var's with Source Splitting (SoS) method.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from decimal import *
import numpy as np


def vdiff_SOS():

    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleGlobalMem import RLEN, ZERO
    import BFM17_POM1D_VrsFnl.src.BFM.General.ModuleService
    # api_bfm (still need to find and write this file)
    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleConstants import SEC_PER_DAY
    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleParam import AssignAirPelFluxesInBFMFlag, p_small
    # mem_Settling (still need to find and write this file)
    from BFM17_POM1D_VrsFnl.src.pom.phys.profTS import PROFTS

    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleMem import D3STATE, D3SOURCE, NO_D3_BOX_STATES, \
        ppO2o,                                                 \
        ppO3c,                                                 \
        ppN1p, ppN3n, ppN4n, ppN5s,                            \
        ppR1c, ppR1n, ppR1p,                                   \
        ppR6c, ppR6n, ppR6p, ppR6s,                            \
        ppP1c, ppP1n, ppP1p, ppP1s,ppP1l,                      \
        ppP2c, ppP2n, ppP2p, ppP2l,                            \
        ppP3c, ppP3n, ppP3p, ppP3l,                            \
        ppP4c, ppP4n, ppP4p, ppP4l,                            \
        ppZ3c, ppZ3n, ppZ3p,                                   \
        ppZ4c, ppZ4n, ppZ4p,                                   \
        ppZ5c, ppZ5n, ppZ5p,                                   \
        ppZ6c, ppZ6n, ppZ6p,                                   \
        sediR6, sediPPY,                                       \
        iiP1, iiP2, iiP3, iiP4,                                \
        n1p, n3n, n4n,n5s,                                     \
        jsurO2o, jsurO3c, o2o,                                 \
        Depth,                                                 \
        iiPhytoPlankton
    from BFM17_POM1D_VrsFnl.src.pom.phys.POMModule import SMOTH, KB, H, DTI, DZR, NRT_o2o, NRT_n1p, NRT_n3n, NRT_n4n, \
      NBCBFM, UMOLBFM, NTP, TIME

    getcontext().prec = 12  # 12-digit precision (ilong)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL VARIABLES
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # COUNTER & FLAGS
    K, M, N, NBC = Decimal()

    # BFM STATE VAR. @ time t-DTI, t, T+DTI RESPECTIVELY
    fbio, ffbio, fbbio = np.empty(KB,dtype=float)

    # SURFACE FLUX STORAGE FOR NUT'S O2 & CO2
    surflux = Decimal()
    botflux = Decimal()

    # RELAXATION VELOCITY FOR NUT'S
    trelax_o2o = Decimal()
    trelax_n1p = Decimal()
    trelax_n3n = Decimal()
    trelax_n4n = Decimal()

    # SEDIMENTATION VELOCITY
    sink = np.empty(KB,dtype=float)
    POCsink = Decimal()
    W1R6 = Decimal()

    # TWICE THE TIME STEP
    DTI2 = Decimal()
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


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: adverte
#
# DESCRIPTION
#
#   SUBROUTINE TO HANDLE THE SINKING OF BFM STATE VAR'S SINKING IS TREATED AS DOWNWARD VERTICAL ADVECTION
#   COMPUTED WITH UPSTREAM FINITE DIFFERENCES.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


def adverte(FB,F,FF,W):

    from BFM17_POM1D_VrsFnl.src.pom.phys.POMModule import KB, DZR, H
    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleGlobalMem import RLEN

    getcontext().prec = 12  # 12-digit precision (ilong)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL VARIABLES
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    FB, F, FF = np.empty(KB,dtype=float)
    W = np.empty(KB,dtype=float)
    DTI2 = Decimal()
    K = int()

    F[KB-1] = F[KB-2]
    FB[KB-1] = FB[KB-2]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Calculate vertical advection. Mind downward velocities are negative
    #   Upwind scheme:
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    FF[0] = DZR[0] * F[0] * W[1]

    for K in range(1,KB-1):
        FF[K] = DZR[K] * (F[K]*W[K+1] - F[K-1]*W[K])

    return




