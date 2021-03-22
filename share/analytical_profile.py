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

def analytical_profile(nlev, vertical_layers, depth1, surface_layer_value, depth2, bottom_layer_value):

    """
    Description: Creates a vertical profile with a value (surface_layer_value) in a
                 surface layer down to depth (depth1) and a bottom layer with a value
                 (bottom_layer_value) from depth (depth2) down to the bottom.
                 Both layers are connected by an intermediate layer reaching from
                 depth1 to depth2 with values varying linearly from surface_layer_value
                 to bottom_layer_value.

    :param nlev:
    :param vertical_layers:
    :param depth1:
    :param surface_layer_value: value at the surface
    :param depth2:
    :param bottom_layer_value: value at bottom
    :return: vertical profile
    """

    # OUTPUT PARAMETERS
    vertical_profile = np.empty(nlev,dtype=float)

    # LOCAL VARIABLES
    alpha = float()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if (depth2 - depth1) > -1E-15:
        alpha = (bottom_layer_value - surface_layer_value) / (depth2 - depth1 + 2E-15)

    else:
        STDERR('*************************************************')
        STDERR('*  Error detected by analytical_profile.py:     *')
        STDERR('*  anz2 should be larger than anz1.             *')
        STDERR('*  Please edit BFM_General.nml or bio_bfm.nml.  *')
        STDERR('*************************************************')

    for i in range(nlev-1,0,-1):
        if vertical_layers[i] < depth1:
            vertical_profile[i] = surface_layer_value

        if alpha < 1E15:
            if depth1 < vertical_layers[i] < depth2:
                vertical_profile[i] = surface_layer_value + alpha * (vertical_layers[i] - depth1)

        if vertical_layers[i] > depth2:
            vertical_profile[i] = bottom_layer_value

    return vertical_profile

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   Copyright 2013 BFM System Team (bfm_st@lists.cmcc.it)
#   Copyright by the GOTM-team under the GNU Public License - www.gnu.org
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
