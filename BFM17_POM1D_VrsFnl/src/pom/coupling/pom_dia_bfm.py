# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: pom_dia_bfm
#
# DESCRIPTION
#
# This routine calculates means and writes the output in diagnostic mode writes also the restart
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from decimal import *
import numpy as np


def pom_dia_bfm(kt,TT):

    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleGlobalMem import RLEN
    # need netcdf_bfm stuff
    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleConstants import SEC_PER_DAY
    from BFM17_POM1D_VrsFnl.src.pom.phys.POMModule import TIME, DTI, IEND

    getcontext().prec = 12  # 12-digit precision

    # TIME IN SECONDS
    localtime = Decimal()

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


#   EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
