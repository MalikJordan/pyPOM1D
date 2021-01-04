# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: pom_ini_bfm
#
# DESCRIPTION
#
# Initialize the BFM in POM
# Main communication of array dimensions between BFM and POM
# Initialisation of variables and netcdf output
# THIS IS THE 1D VERSION
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from decimal import *
import numpy as np


def pom_ini_bfm_1D():

    import BFM17_POM1D_VrsFnl.src.BFM.General.ModuleMem

    from BFM17_POM1D_VrsFnl.src.pom.phys.POMModule import DTI, IEND, KB, H, DZ, IHOTST, ALAT, ALON, ZZ, DZZ

    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleGlobalMem import ZERO

    import BFM17_POM1D_VrsFnl.src.share.api_bfm

    from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleService import savef, deltat, nitend

    from BFM17_POM1D_VrsFnl.src.share.netcdf_bfm import init_netcdf_bfm, init_save_bfm
    from BFM17_POM1D_VrsFnl.src.share.netcdf_bfm import init_netcdf_rst_bfm, read_rst_bfm

    getcontext().prec = 12  # 12-digit precision (ilong)

    # COUNTER
    K, i = int()

    # READING UNITS
    namlst = 10
    unit = 11

    # OCEAN WET POINTS (IN 1D = 1)
    ocepoint = int()
    # OCEAN SURFACE POINTS (IN 1D = 1)
    surfpoint = int()
    # OCEAN BOTTOM POINTS (IN 1D = 1)
    botpoint = int()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   MODEL TIMESTEP
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    deltat = DTI

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   ITERATIONS NEEDED FOR AN IDAYS SIMULATIONS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    nitend = IEND

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   out_delta = saving frequency
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    out_delta = savef

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SET THE BFM DIMENSIONS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #
    #   SINCE THIS IS A 1D IMPLEMENTATION THE SIZE OF THE ARRAYS IS SHRUNK
    #   TO A 1D VECTOR  ALONG the "VERTICAL" DIMENSION (1, 1, KB-1)
    #
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    NO_BOXES = KB - 1
    NO_BOXES_X = 1
    NO_BOXES_Y = 1
    NO_BOXES_Z = NO_BOXES
    NO_BOXES_XY = 1

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   ALLOCATE MASKS FOR  ARRAY PACKING
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    SEAmask = np.empty(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z)
    BOTmask = np.empty(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z)
    SRFmask = np.empty(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z)

    SEAmask = True
    BOTmask = True
    SRFmask = True

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   ALLOCATE ANCILLARY PACK MASK
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ZEROS = np.empty(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z)
    ZEROS = ZERO

    if INCLUDE_BEN is not None:

        NO_STATES = NO_D3_BOX_STATES * NO_BOXES + NO_D2_BOX_STATES_BEN
        NO_BOXES_Z_BEN = 0
        NO_BOXES_BEN = NO_BOXES_XY * NO_BOXES_Z_BEN
        NO_STATES_BEN = NO_BOXES_BEN * NO_D2_BOX_STATES_BEN

    else:

        NO_STATES = NO_D3_BOX_STATES * NO_BOXES

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   COMPRESSED COORDINATES FOR netcdf OUTPUT
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    lon_len = NO_BOXES_X
    lat_len = NO_BOXES_Y

    ocepoint = np.empty(NO_BOXES)
    ocepoint = 1

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   PREPARES THE ARRAY CONTAINING THE INDICES OF THE ELEMENTS IN PELAGIC
    #   BFM 1D ARRAYS THAT HAVE A BENTHIC LAYER
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    BOTindices = np.empty(NO_BOXES_XY)
    BOTindices = KB - 1

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   PREPARES THE ARRAY CONTAINING THE INDICES OF THE ELEMENTS IN PELAGIC
    #   BFM 1D ARRAYS THAT HAVE A SURFACE LAYER
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    SRFindices = np.empty(NO_BOXES_XY)
    SRFindices = 1

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INITIALIZE ANCILLARY ARRAYS FOR OUTPUT
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    init_bfm(namlst)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INITIALIZE STATE VARIABLE NAMES AND DIAGNOSTICS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    set_var_info_bfm

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   ALLOCATE MEMORY AD GIVE INITIAL VALUES
    #   THE ARGUMENT LIST IS MANDATORY WITH BFM
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    init_var_bfm(namlst,'bfm.nml',unit,bio_setup)

    # THE LEADING RESTART FLAG IS IHOTST, SO THE VALUE IN bfm_init IS OVERWRITTEN
    bfm_init = IHOTST

    # SET THE THICKNESS OF THE KB-1 LAYERS
    for K in range(0,NO_BOXES_Z):
        Depth[K] = DZ[K] * H

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INITIALIZE netcdf OUTPUT
    #   MAV: CHECK THE MAPPING OF OCEPOINT
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    from BFM17_POM1D_VrsFnl.src.share.calcmean_bfm import calcmean_bfm
    calcmean_bfm(INIT)

    init_netcdf_bfm((title=out_title,start_time='01-01-0000',
                    time_unit=0,lat=alat,lon=alon,Z=ZZ,DZ=DZZ,
                    oceanpoint=ocepoint,
                    surfacepoint=(/(i,i=1,NO_BOXES_XY)/),
                    bottompoint=(/(i,i=1,NO_BOXES_XY)/)))

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   ALLOCATE AND INITIALIZE ADDITIONAL INTEGRATION ARRAYS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    D3STATEB = np.empty(NO_D3_BOX_STATES,NO_BOXES)
    D3STATEB = ZERO

    if INCLUDE_BEN is not None:
        D2STATEB_BEN = np.empty(NO_D2_BOX_STATES_BEN,NO_BOXES_XY)
        D2STATEB_BEN = ZERO

    # DEFINE INITIAL CONDITIONS
    from BFM17_POM1D_VrsFnl.src.pom.coupling.set_initial_conditions import set_initial_conditions
    set_initial_conditions()

    init_save_bfm()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INITIALIZE PRIOR TIME STEP FOR LEAP-FROG
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    D3STATEB = D3STATE

    if INCLUDE_BEN is not None:
        D2STATEB_BEN = D2STATE_BEN

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   READ RESTART (Bfm_init = 1 in bfm.nml)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    init_netcdf_rst_bfm(rst_fname,start_time='01-01-0000',
                        time_unit=0,lat=alat,lon=alon,Z=ZZ,DZ=DZZ,
                        oceanpoint=ocepoint,
                        surfacepoint=(/(i,i=1,NO_BOXES_XY)/),
                        bottompoint=(/(i,i=1,NO_BOXES_XY)/))

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   READ RESTART FILE (IF FLAG)
    #   OVERWRITE PREVIOUS INITIALIZATION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    print('before read rst')
    if bfm_init == 1:
        read_rst_bfm(rst_fname)

    return

#   EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
