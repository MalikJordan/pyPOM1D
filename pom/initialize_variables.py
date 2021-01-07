# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
import numpy as np
from os import path
import f90nml

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# SUBROUTINE: read_pom_input
#
# DESCRIPTION: Opens forcing files reading path specified in pom_input nml.
# (formerly opendat)
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def read_pom_input():

    from main_pombfm1d import current_path
    from pom.modules import KB
    from pom.modules import path_error

    # PATHS TO INPUT DATA FILES
    wind_input = current_path + '/inputs/POM_BFM17/monthly_surf_wind_stress_bermuda_killworth2.da'
    surface_s_input = current_path + '/inputs/POM_BFM17/monthly_surf_salt_bermuda_150m_killworth2.da'
    radiance_input = current_path + '/inputs/POM_BFM17/monthly_surf_qs_bermuda_killworth2.da'
    ism_input = current_path + '/inputs/POM_BFM17/monthly_clima_ISM_150m_bermuda_killworth.da'
    sal_input = current_path + '/inputs/POM_BFM17/monthly_clima_salt_150m_bermuda_killworth2.da'
    temp_input = current_path + '/inputs/POM_BFM17/monthly_clima_temp_150m_bermuda_killworth2.da'
    w_input = current_path + '/inputs/POM_BFM17/monthly_clima_w_150m_bermuda_ekman.da'
    weddy_input1 = current_path + '/inputs/POM_BFM17/bimonthly_random_eddy_w_150m_bermuda_norm1.da'
    weddy_input2 = current_path + '/inputs/POM_BFM17/bimonthly_random_eddy_w_150m_bermuda_norm2.da'
    Sprofile_input = current_path + '/inputs/POM_BFM17/init_prof_S_150m_bermuda_killworth2.da'
    Tprofile_input = current_path + '/inputs/POM_BFM17/init_prof_T_150m_bermuda_killworth2.da'
    heat_input = current_path + '/inputs/POM_BFM17/monthly_surf_rad_bermuda_killworth2.da'
    surfNut_input = current_path + '/inputs/POM_BFM17/NutrientsARPAOGS.da'
    bottNut_input = current_path + '/inputs/POM_BFM17/monthly_bott_nut_bermuda_150m_killworth.da'

    # LENGTH OF INPUT ARRAYS
    array_length = 13

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   WIND SPEED (u,v)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(wind_input):
        path_error(wind_input)
    open11 = np.fromfile(wind_input,dtype=float)
    WSU1   = np.zeros(array_length)
    WSV1   = np.zeros(array_length)
    for k in range(0,array_length):
        WSU1[k] = open11[2*k + 0]
        WSV1[k] = open11[2*k + 1]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SURFACE SALINITY
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(surface_s_input):
        path_error(surface_s_input)
    open13 = np.fromfile(surface_s_input,dtype=float)
    SSS1   = np.zeros(array_length)
    for k in range(0,array_length):
        SSS1[k] = open13[k]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   RADIANCE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(radiance_input):
        path_error(radiance_input)
    open17 = np.fromfile(radiance_input,dtype=float)
    SLUX1  = np.zeros(array_length)
    for k in range(0,array_length):
        SLUX1[k] = open17[k]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INORGANIC SUSPENDED MATTER
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(ism_input):
        path_error(ism_input)
    open19 = np.fromfile(ism_input,dtype=float)
    ISM1   = np.zeros((array_length,KB))
    for k in range(0,array_length):
        for x in range(0,KB):
            ISM1[k][x] = open19[KB*k + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SALINITY CLIMATOLOGY (DIAGNOSTIC MODE)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(sal_input):
        path_error(sal_input)
    open20 = np.fromfile(sal_input,dtype=float)
    SCLIM1 = np.zeros((array_length,KB))
    for k in range(0,array_length):
        for x in range(0,KB):
            SCLIM1[k][x] = open20[KB*k + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   TEMPERATURE CLIMATOLOGY (DIAGNOSTIC MODE)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(temp_input):
        path_error(temp_input)
    open15 = np.fromfile(temp_input,dtype=float)
    TCLIM1 = np.zeros((array_length,KB))
    for k in range(0,array_length):
        for x in range(0,KB):
            TCLIM1[k][x] = open15[KB*k + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   GENERAL CIRCULATION W VELOITY CLIMATOLOGY
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(w_input):
        path_error(w_input)
    open215 = np.fromfile(w_input,dtype=float)
    WCLIM1  = np.zeros((array_length,KB))
    for k in range(0,array_length):
        for x in range(0,KB):
            WCLIM1[k][x] = open215[KB*k + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INTERMEDIATE EDDY W VELOCITY 1
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(weddy_input1):
        path_error(weddy_input1)
    open216 = np.fromfile(weddy_input1,dtype=float)
    WEDDY1  = np.zeros((array_length,KB))
    for k in range(0,array_length):
        for x in range(0,KB):
            WEDDY1[k][x] = open216[KB*k + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   INTERMEDIATE EDDY W VELOCITY 2
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(weddy_input2):
        path_error(weddy_input2)
    open217 = np.fromfile(weddy_input2,dtype=float)
    WEDDY2  = np.zeros((array_length,KB))
    for k in range(0,array_length):
        for x in range(0,KB):
            WEDDY2[k][x] = open217[KB*k + x]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SALINITY INITIAL PROFILE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(Sprofile_input):
        path_error(Sprofile_input)
    open29 = np.fromfile(Sprofile_input,dtype=float)
    SB     = np.zeros(KB)
    for k in range(0,KB):
        SB[k] = open29[k]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   TEMPERATURE INITIAL PROFILE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(Tprofile_input):
        path_error(Tprofile_input)
    open10 = np.fromfile(Tprofile_input,dtype=float)
    TB     = np.zeros(KB)
    for k in range(0,KB):
        TB[k] = open10[k]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   HEAT FLUX
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(heat_input):
        path_error(heat_input)
    open21 = np.fromfile(heat_input,dtype=float)
    SWRAD1  = np.zeros(array_length)
    WTSURF1 = np.zeros(array_length)
    QCORR1  = np.zeros(array_length)
    for k in range(0,array_length):
        SWRAD1[k]  = open21[3*k + 0]
        WTSURF1[k] = open21[3*k + 1]
        QCORR1[k]  = open21[3*k + 2]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SURFACE NUTRIENTS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(surfNut_input):
        path_error(surfNut_input)
    open18  = np.fromfile(surfNut_input,dtype=float)
    NO3_s1  = np.zeros(array_length)
    NH4_s1  = np.zeros(array_length)
    PO4_s1  = np.zeros(array_length)
    SIO4_s1 = np.zeros(array_length)
    for k in range(0,array_length):
        NO3_s1[k]  = open18[4*k + 0]
        NH4_s1[k]  = open18[4*k + 1]
        PO4_s1[k]  = open18[4*k + 2]
        SIO4_s1[k] = open18[4*k + 3]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   BOTTOM NUTRIENTS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if not path.exists(bottNut_input):
        path_error(bottNut_input)
    open300 = np.fromfile(bottNut_input,dtype=float)
    O2_b1   = np.zeros(array_length)
    NO3_b1  = np.zeros(array_length)
    PO4_b1  = np.zeros(array_length)
    PON_b1  = np.zeros(array_length)
    for k in range(0,array_length):
        O2_b1[k]  = open300[4*k + 0]
        NO3_b1[k] = open300[4*k + 1]
        PO4_b1[k] = open300[4*k + 2]
        PON_b1[k] = open300[4*k + 3]

    return WSU1, WSV1, SSS1, SLUX1, ISM1, SCLIM1, TCLIM1, WCLIM1, WEDDY1, WEDDY2, SB, TB, \
               SWRAD1, WTSURF1, QCORR1, NO3_s1, NH4_s1, PO4_s1, SIO4_s1, O2_b1, NO3_b1, PO4_b1, PON_b1


# WSU1, WSV1, SSS1, SLUX1, ISM1, SCLIM1, TCLIM1, WCLIM1, WEDDY1, WEDDY2, SB, TB, \
#     SWRAD1, WTSURF1, QCORR1, NO3_s1, NH4_s1, PO4_s1, SIO4_s1, O2_b1, NO3_b1, PO4_b1, PON_b1 = read_pom_input()

#   EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
