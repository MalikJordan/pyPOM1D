from cppdefs import *
import numpy as np

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: ANALYTICAL_PROFILE
#
# DESCRIPTION:  This routine creates a vertical profile {\tt prof} with value
#               {\tt v1} in a surface layer down to depth {\tt z1} and a bottom
#               layer of value {\tt v2} reaching from depth {\tt z2} down to the bottom.
#               Both layers are connected by an intermediate layer reaching from {\tt z1}
#               to {\tt z2} with values linearly varying from {\tt v1} to {\tt v2}.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def analytical_profile(nlev,z,z1,v1,z2,v2,prof):

    from pom.modules import RLEN, zero, bfm_lwp, LOGUNIT

    # INPUT PARAMETERS
    # nlev = int()
    # z = np.empty(nlev,dtype=float)
    # z1 = float()
    # v1 = float()
    # z2 = float()
    # v2 = float()

    # OUTPUT PARAMETERS
    prof = np.empty(nlev,dtype=float)

    # LOCAL VARIABLES
    i = int()
    alpha = float()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if z2 - z1 > -1E-15:
        alpha = (v2 - v1) / (z2 - z1 + 2E-15)

    else:
        STDERR('*************************************************')
        STDERR('*  Error detected by analytical_profile.F90:    *')
        STDERR('*  anz2 should be larger than anz1.             *')
        STDERR('*  Please edit BFM_General.nml or bio_bfm.nml.  *')
        STDERR('*************************************************')

    for i in range(nlev-1,0,-1):
        if z[i] < z1:
            prof[i] = v1

        if alpha < 1E15:
            if z1 < z[i] < z2:
                prof[i] = v1 + alpha*(z[i] - z1)

        if z[i] > z2:
            prof[i] = v2

    return prof

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   Copyright 2013 BFM System Team (bfm_st@lists.cmcc.it)
#   Copyright by the GOTM-team under the GNU Public License - www.gnu.org
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
