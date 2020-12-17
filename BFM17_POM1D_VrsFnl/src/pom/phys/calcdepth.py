# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: Calcdepth
#
# DESCRIPTION
#
#   This subroutine establishes the vertical resolution log distributions
#   at the top and bottom, and a linear distribution between KL1 and KL2.
#   Default values: KL1 = .3*KB AND KL2 = KB-2.
#   Yields a log distribution at the top and none at the bottom.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
import numpy as np
from decimal import *


def CALCDEPTH(Z,ZZ,DZ,DZZ,KB,KL1,KL2):

    getcontext().prec = 12

    # DZ,DZZ,Z,ZZ = np.empty(KB,dtype=float)

    BB = float(KL2-KL1) + 4.
    CC = float(KL1) - 2.
    DEL1 = 2./BB/np.exp(.693147*float(KL1-2))
    DEL2 = 2./BB/np.exp(.693147*float(KB-KL2-1))
    Z[0] = 0.
    ZZ[0] = -DEL1/2.

    for K in range(1,KL1 - 2):
        Z[K] = -DEL1 * np.exp(.693147 * float(K - 2))
        ZZ[K] = -DEL1 * np.exp(.693147 * (float(K) - 1.5))

    for K in range(KL1 - 2,KL2 + 1):
        Z[K] = - (float(K) - CC) / BB
        ZZ[K] = - (float(K) - CC + 0.5) / BB

    for K in range(0,KB - 1):
        DZ[K] = Z[K] - Z[K+1]
        DZZ[K] = ZZ[K] - ZZ[K+1]

    return

# EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
