# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: get_TS_IC
#
# DESCRIPTION
#
# This subroutine opens and reads files containing the T&S initial conditions
# Files are read in direct access mode reading path specified in pom_input nml
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from decimal import *
import numpy as np


def get_TS_IC():

    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleGlobalMem import error_msg_prn, NML_OPEN, NML_READ
    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleService import wind_input, surfaceS_input, radiance_input, ism_input, \
        Sal_input, Temp_input, W_input, Weddy_input1, Weddy_input2, Sprofile_input, Tprofile_input, heat_input, \
        surfNut_input, bottNut_input, ISM, read_restart
    from BFM17_POM1D_VrsFnl.src.pom.phys.POMModule import KB, T, TB, S, SB
    import f90nml

    # LOOP COUNTER
    K = int()

    # RECORD LENGTH
    RLENGTH = int()

    # NAMELIST READING UNIT
    namlst = 10

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   OPEN NML WITH FORCING DATA PATH SPECIFIED AND READ DATA
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    pom_input = f90nml.read('pom_input.nml')  # OPEN NAMLST

    for K in range(0,KB):
        SB[0] = pom_input[Sprofile_input]
        TB[0] = pom_input[Tprofile_input]

    T[:] = TB[:]
    S[:] = SB[:]

    error_msg_prn(NML_OPEN,"get_TS_IC.F90","pom_input.nml")
    error_msg_prn(NML_READ, "get_TS_IC.F90", "pom_input")

    return

#   EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

