# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: CalcLightDistribution
#
# DESCRIPTION
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
import numpy as np


def CalcLightDistribution():

    BoxNumberZ, NO_BOXES_Z, BoxNumberY, NO_BOXES_Y, BoxNumberX, NO_BOXES_X, BoxNumber, BoxNumberXY, \
    EIR, SUNQ, xEPS, Depth, NO_BOXES_XY, NO_BOXES = None

    HOURS_PER_DAY, E2W = None
    p_PAR = None

    EIR[0] = EIR[0] * p_PAR / E2W
    for BoxNumberZ in range(1, NO_BOXES_Z):
        box_no = BoxNumberZ
        EIR[box_no] = EIR[box_no - 1] * np.exp(- 1.0 * xEPS[box_no - 1] * Depth[box_no - 1])

    return

# EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
