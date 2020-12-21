# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: get_rst
#
# DESCRIPTION
#
# This subroutine opens and reads the fort file containing the restart.
# The file path is specified in pom_input.nml
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from decimal import *
import numpy as np


def get_rst():

    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleGlobalMem import error_msg_prn, NML_OPEN, NML_READ
    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleService import wind_input,surfaceS_input,radiance_input, \
        ism_input,Sal_input,Temp_input,W_input,Weddy_input1, Weddy_input2,Sprofile_input,Tprofile_input,heat_input, \
        surfNut_input,bottNut_input,ISM,read_restart
    from BFM17_POM1D_VrsFnl.src.pom.phys.POMModule import TIME0, U, UB, V, VB, T, TB, S, SB, Q2, Q2B, Q2L, Q2LB, \
        KH, KM, KQ, L, WUBOT, WVBOT, RHO
    import f90nml

    getcontext().prec = 12  # 12-digit precision (ilong)

    # NAMELIST READING UNIT
    namlst = 10

    # RECORD LENGTH
    RLENGTH = Decimal()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   OPEN NML WITH FORCING DATA PATH SPECIFIED AND READ DATA
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    pom_input = f90nml.read('pom_input.nml')  # OPEN NAMLST

    TIME0 = pom_input[read_restart][0]
    U = pom_input[read_restart][1]
    UB = pom_input[read_restart][2]
    V = pom_input[read_restart][3]
    VB = pom_input[read_restart][4]
    T = pom_input[read_restart][5]
    TB = pom_input[read_restart][6]
    S = pom_input[read_restart][7]
    SB = pom_input[read_restart][8]
    Q2 = pom_input[read_restart][9]
    Q2B = pom_input[read_restart][10]
    Q2L = pom_input[read_restart][11]
    Q2LB = pom_input[read_restart][12]
    KH = pom_input[read_restart][13]
    KM = pom_input[read_restart][14]
    KQ = pom_input[read_restart][15]
    L = pom_input[read_restart][16]
    WUBOT = pom_input[read_restart][17]
    WVBOT = pom_input[read_restart][18]
    RHO = pom_input[read_restart][19]

    error_msg_prn(NML_OPEN, 'opendat.py', 'pom_input.nml')
    error_msg_prn(NML_READ, 'opendat.py', 'pom_input')

    return


#   EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
