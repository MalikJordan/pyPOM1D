# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: restart_BFM_inPOM
#
# DESCRIPTION
#
# This routine writes the file needed for the BFM restart
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from decimal import *
import numpy as np


def restart_BFM_inPOM():

    from BFM17_POM1D_VrsFnl.src.share.netcdf_bfm import save_bfm, close_ncdf, ncid_bfm
    from BFM17_POM1D_VrsFnl.src.share.netcdf_bfm import save_rst_bfm, ncid_rst
    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleConstants import SEC_PER_DAY
    from BFM17_POM1D_VrsFnl.src.pom.phys.POMModule import TIME
    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleGlobalMem import RLEN
    import BFM17_POM1D_VrsFnl.src.BFM.General.ModuleMem
    import BFM17_POM1D_VrsFnl.src.share.api_bfm

    getcontext().prec = 12  # 12-digit precision (ilong)

    localtime = Decimal()

    # TIME ELLAPSED (IN SECONDS) SINCE THE BEGINNING OF THE SIMULATION
    localtime = TIME*SEC_PER_DAY

    # WRITE RESTART
    save_rst_bfm(localtime)

    # CLOSE OUTPUT AND RESTART FILES
    close_ncdf(ncid_rst)
    close_ncdf(ncid_bfm)

    print('NETCDF RESTART WRITTE, TIME-->', TIME)

    return

#   EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
