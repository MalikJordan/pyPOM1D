import numpy as np
from pom.constants import vertical_layers

def detritus_sedimentation():

    # FROM namelists_bfm
    # p_rR6m        [m/d]   detritus sinking rate
    # p_burvel_R6   [m/d]   Bottom burial velocity for detritus

    p_rR6m = 1.
    p_burvel_R6 = 1.
    # p_burvel_R6 = 1.5

    # FROM PelGlobal.F90 (145-148)
    detritus_sedimentation_rate = p_rR6m * np.ones(vertical_layers-1)
    try:
        BFM_POM
    except NameError:
        BFM_POM = False
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
        phyto_sedimentation_rates[:,iiPhytoPlankton] = p_rPIm[i]
        try:
            BFM_POM
        except NameError:
            BFM_POM = False
        if not BFM_POM:
            phyto_sedimentation_rates[vertical_layers-2,i] = p_burvel_PI

    return phyto_sedimentation_rates
