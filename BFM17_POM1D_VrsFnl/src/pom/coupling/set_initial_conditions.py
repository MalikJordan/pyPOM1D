# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: set_initial_conditions
#
# DESCRIPTION
#
# This routine assigns initial conditioons of biochemical variables in POM
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from decimal import *
import numpy as np
from f90nml import namelist

from BFM17_POM1D_VrsFnl.src.pom.phys.POMModule import KB


def set_initial_conditions():

    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleGlobalMem import ZERO,NML_OPEN,NML_READ,NMLUNIT,error_msg_prn
    import BFM17_POM1D_VrsFnl.src.BFM.General.ModuleMem
    import BFM17_POM1D_VrsFnl.src.pom.phys.POMModule
    import BFM17_POM1D_VrsFnl.src.BFM.General.ModuleService
    import f90nml

    getcontext().prec = 12  # 12-digit precision (ilong)

    ib = int()

    dd, d1, d2, d1cc, d1cn, d1cp, d1cs, d1ci = Decimal()
    d2cc = np.empty(KB,dtype=float)

    p_nRc = Decimal(0.0126)
    p_pRc = Decimal(0.7862E-03)
    p_sRc = Decimal(0.0118)
    p_iRc = Decimal(1. / 25.)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   DEFINITION OF INITIAL PELAGIC (D3) STATE VARIABLES
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    O2o0, N1p0, N3n0, N4n0, O4n0, N5s0, N6r0, B1c0, B1n0, B1p0, \
        P1c0, P1n0, P1p0, P1l0, P1s0, P2c0, P2n0, P2p0, P2l0, P3c0, P3n0, P3p0, \
        P3l0, P4c0, P4n0, P4p0, P4l0, Z3c0, Z3n0, Z3p0, Z4c0, Z4n0, Z4p0, Z5c0, \
        Z5n0, Z5p0, Z6c0, Z6n0, Z6p0, R1c0, R1n0, R1p0, R2c0, R3c0, R6c0, R6n0, \
        R6p0, R6s0, O3c0, O3h0 = Decimal()

    bfm_init_nml = f90nml.read('bfm_init_nml.nml')

    if INCLUDE_BEN is not None:

        G3c0, G3h0, G13c0, G13h0, G23c0, G23h0, Y1c0, Y1n0, Y1p0, \
            Y2c0, Y2n0, Y2p0, Y3c0, Y3n0, Y3p0, Y4c0, Y4n0, Y4p0, Y5c0, Y5n0, Y5p0, \
            Q6c0, Q6n0, Q6p0, Q6s0, Q1c0, Q1n0, Q1p0, Q11c0, Q11n0, Q11p0, H1c0, H1n0, \
            H1p0, H2c0, H2n0, H2p0, K1p0, K11p0, K21p0, K4n0, K14n0, K24n0, K6r0, \
            K16r0, K26r0, K3n0, K5s0, G2o0, G4n0, D1m0, D2m0, D6m0, D7m0, D8m0, D9m0 = Decimal()

        bfm_init_nml_ben = f90nml.read('bfm_init_nml_ben.nml')

    NMLUNIT = f90nml.read('BFM_General.nml')

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   DEFINITION OF BIOGEOCHEMICAL GLOBAL VARIABLES
    #   IrrOPT  in the equation of Steele and light
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    eir[:] = ZERO

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   DEFINITION OF GENERAL PELAGIC STATE VARIABLES:
    #   PELAGIC GASES
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    get_IC()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   PELAGIC NUTRIENTS (mMol / m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    n5s[:] = N5s0
    O4n[:] = O4n0
    N6r[:] = N6r0

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   PELAGIC DETRITUS (respectively mg C/m3 mMol N/m3 mMOL P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    R6n[:] = r6c[:]*p_nRc
    R6p[:] = r6c[:]*p_pRc
    R6s[:] = r6c[:]*p_sRc
    R2c[:] = R2c0
    R3c[:] = R3c0

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   DISSOLVED ORGANIC MATTER
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    R1n[:] = R1c[:]*p_nRc * 0.5
    R1p[:] = R1c[:]*p_pRc * 0.5

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   STATE VARIABLES FOR PHYTOPLANKTON MODEL
    #   PELAGIC DIATOMS (respectively mg C/m3 mMol N/m3 mMOL P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    d1cc = P1c0
    d1cn = d1cc * p_nRc
    d1cp = d1cc * p_pRc
    d1cs = d1cc * p_sRc
    d1ci = d1cc * p_iRc

    if CalcPhytoPlankton[0] is not None:
        P1c[:] = d1cc
        P1n[:] = d1cn
        P1p[:] = d1cp
        P1s[:] = d1cs
        P1l[:] = d1ci
    else:
        P1c[:] = ZERO
        P1n[:] = ZERO
        P1p[:] = ZERO
        P1s[:] = ZERO
        P1l[:] = ZERO

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   PELAGIC FLAGELLATES (respectively mg C/m3 mMol N/m3 mMOL P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    d1cc = P2c0
    d1cn = d1cc * p_nRc
    d1cp = d1cc * p_pRc
    d1ci = d1cc * p_iRc * 0.5

    if CalcPhytoPlankton[1] is not None:
        d2cc[:] = P2c[:]
        P2n[:] = d2cc[:]*p_nRc
        P2p[:] = d2cc[:]*p_pRc
        P2l[:] = d2cc[:]*p_iRc * 0.5
    else:
        P2c[:] = ZERO
        P2n[:] = ZERO
        P2p[:] = ZERO
        P2l[:] = ZERO

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   PICOPHYTOPLANKTON (respectively mg C/m3 mMol N/m3 mMOL P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    d1cc = P3c0
    d1cn = d1cc * p_nRc
    d1cp = d1cc * p_pRc
    d1ci = d1cc * p_iRc * 0.5

    if CalcPhytoPlankton[2] is not None:
        P3c[:] = d1cc
        P3n[:] = d1cn
        P3p[:] = d1cp
        P3l[:] = d1ci
    else:
        P3c[:] = ZERO
        P3n[:] = ZERO
        P3p[:] = ZERO
        P3l[:] = ZERO

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LARGE PHYTOPLANKTON (respectively mg C/m3 mMol N/m3 mMOL P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    d1cc = P4c0
    d1cn = d1cc * p_nRc
    d1cp = d1cc * p_pRc
    d1ci = d1cc * p_iRc * 0.5

    if CalcPhytoPlankton[3] is not None:
        P4c[:] = d1cc
        P4n[:] = d1cn
        P4p[:] = d1cp
        P4l[:] = d1ci
    else:
        P4c[:] = ZERO
        P4n[:] = ZERO
        P4p[:] = ZERO
        P4l[:] = ZERO

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   STATE VARIABLES FOR MESOZOOPLANKTON MODEL
    #   CARNIVOROUS MESOZOOPLANKTON ( mg C/m3 )
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if CalcMesoZooPlankton[0] is not None:
        Z3c[:] = Z3c0
        Z3n[:] = Z3c[:]*p_nRc
        Z3p[:] = Z3c[:]*p_pRc
    else:
        Z3c[:] = ZERO
        Z3n[:] = ZERO
        Z3p[:] = ZERO

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   OMNIVOROUS MESOZOOPLANKTON ( mg C/m3 )
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if CalcMesoZooPlankton[1] is not None:
        Z4c[:] = Z4c0
        Z4n[:] = Z4c[:]*p_nRc
        Z4p[:] = Z4c[:]*p_pRc
    else:
        Z4c[:] = ZERO
        Z4n[:] = ZERO
        Z4p[:] = ZERO

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   STATE VARIABLES FOR MICROZOOPLANKTON MODEL
    #   PELAGIC MICROZOOPLANKTON (respectively mg C/m3 mMol N/m3 mMOL P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if CalcMicroZooPlankton[0] is not None:
        Z5n[:] = Z5c[:]*p_nRc
        Z5p[:] = Z5c[:]*p_pRc
    else:
        Z5c[:] = ZERO
        Z5n[:] = ZERO
        Z5p[:] = ZERO

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   HETEROTROPHIC FLAGELLATES (respectively mg C/m3 mMol N/m3 mMOL P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if CalcMicroZooPlankton[1] is not None:
        Z6c[:] = Z6c0
        Z6n[:] = Z6c[:]*p_nRc
        Z6p[:] = Z6c[:]*p_pRc
    else:
        Z6c[:] = ZERO
        Z6n[:] = ZERO
        Z6p[:] = ZERO

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   STATE VARIABLES FOR PELAGIC BACTERIA MODEL B1
    #   PELAGIC BACTERIA (respectively mg C/m3 mMol N/m3 mMOL P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if CalcPelBacteria[0] is not None:
        B1c[:] = B1c0
        B1n[:] = B1c[:]*p_nRc
        B1p[:] = B1c[:]*p_pRc
    else:
        B1c[:] = ZERO
        B1n[:] = ZERO
        B1p[:] = ZERO

    if INCLUDE_BEN is not None:
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   STATE VARIABLES FOR BENTHIC MODULES
        #   ZOOBENTHOS (respectively mg C/m3 mMol N/m3 mMOL P/m3)
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        if CalcBenOrganisms[0] is not None:
            Y1c[:] = Y1c0
        else:
            Y1c[:] = ZERO

        if CalcBenOrganisms[1] is not None:
            Y2c[:] = Y2c0
        else:
            Y2c[:] = ZERO

        if CalcBenOrganisms[2] is not None:
            Y3c[:] = Y3c0
        else:
            Y3c[:] = ZERO

        if CalcBenOrganisms[3] is not None:
            Y4c[:] = Y4c0
        else:
            Y4c[:] = ZERO

        if CalcBenOrganisms[4] is not None:
            Y5c[:] = Y5c0
        else:
            Y5c[:] = ZERO

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   BENTHIC NUTRIENTS
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        K5s[:] = 20.75
        K6r[:] = K6r0
        K4n[:] = K4n0
        K14n[:] = K14n0
        K24n[:] = K24n0
        K1p[:] = K1p0
        K11p[:] = K11p0
        K21p[:] = K21p0
        K3n[:] = K3n0

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   BENTHIC DETRITUS (respectively mg C/m3 mMol N/m3 mMOL P/m3)
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        Q1c[:] = Q1c0
        Q11c[:] = Q11c0

        if IMPFLUX is not None:
            Q6c[:] = 1.E9
            Q6n[:] = 1.E9
            Q6p[:] = 1.E9
            Q6s[:] = 1.E9
        else:
            Q6c[:] = Q6c0
            Q6n[:] = Q6n0
            Q6p[:] = Q6p0
            Q6s[:] = Q6s0

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   GASES
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        G2o[:] = G2o0
        G3c[:] = G3c0
        G4n[:] = 37.0

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   LAYERS
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        D1m[:] = D1m0
        D2m[:] = D2m0
        D6m[:] = D6m0
        D7m[:] = D7m0
        D8m[:] = D8m0
        D9m[:] = D9m0

    error_msg_prn(NML_OPEN, "InitParam.f90", "BFM_General.nml")
    error_msg_prn(NML_READ, "InitParam.f90", "Param_parameters")
    error_msg_prn(NML_READ, "InitParam.f90", "Param_parameters_ben")

    return

#   EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
