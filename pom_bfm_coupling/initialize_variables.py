import numpy as np
from include import IMPFLUX, INCLUDE_BEN
from pom.constants import current_path, seconds_per_day, twice_the_timestep, vertical_layers
from inputs import params_POMBFM
from bfm.constants import num_d3_box_states #, num_d2_box_states_ben
from bfm.variable_info import set_variable_info_bfm
from os import path
from inputs.namelist_input_data import o2o0, n1p0, n3n0, n4n0, n5s0, n6r0, o3c0, o3h0, o4n0, p1c0, p2c0, p3c0, p4c0, \
    z3c0, z4c0, z5c0, z6c0, b1c0, r1c0, r2c0, r3c0, r6c0, \
    y1c0, y2c0, y3c0, y4c0, y5c0, h1c0, h2c0, k1p0, k11p0, k21p0, k4n0, k14n0, k24n0, k3n0, k5s0, k6r0, \
    d1m0, d2m0, d6m0, d7m0, d8m0, d9m0, q6c0, q6n0, q6p0, q6s0, q1c0, q11c0, g2o0, g3c0, g13c0, g23c0, g3h0, g13h0, g23h0, \
    calcphytoplankton, calcmesozooplankton, calcmicrozooplankton, calcpelbacteria, calcbenorganisms, calcbenbacteria

if INCLUDE_BEN:
    from bfm.constants import num_d2_box_states_ben

def initialize_bfm_in_pom(vertical_grid):

    ocean_wet_points = 1
    ocean_surface_points = 1
    ocean_bottom_points = 1

    # BFM dimensions (1, 1, vertical_layers-1)
    num_boxes = vertical_layers - 1
    num_boxes_x = 1
    num_boxes_y = 1
    num_boxes_z = num_boxes
    num_boxes_xy = 1

    SEAmask = np.full((num_boxes_x,num_boxes_y,num_boxes_z),True)
    BOTmask = np.full((num_boxes_x,num_boxes_y,num_boxes_z),True)
    SRFmask = np.full((num_boxes_x,num_boxes_y,num_boxes_z),True)

    zeros = np.zeros((num_boxes_x,num_boxes_y,num_boxes_z))

    # try:
    #     INCLUDE_BEN
    # except NameError:
    #     INCLUDE_BEN = False
    # else:
    #     INCLUDE_BEN = True
    if INCLUDE_BEN:
        num_states = num_d3_box_states * num_boxes + num_d2_box_states_ben
        num_boxes_z_ben = 0
        num_boxes_ben = num_boxes_xy * num_boxes_z_ben
        num_states_ben = num_boxes_ben * num_d2_box_states_ben
    else:
        num_states = num_d3_box_states * num_boxes

    # Compresses coordinates
    lonitude_length = num_boxes_x
    latitude_length = num_boxes_y

    # Array containing the indices of the elements in pelagic BFM 1D arrays that have a benthic layer
    bottom_indicies = (vertical_layers-1) * np.ones(num_boxes_xy)

    # Array containing the indices of the elements in pelagic BFM 1D arrays that have a surface layer
    surface_indices = np.ones(num_boxes_xy)

    # Initialize ancillary arrays for output
    # init_bfm()

    # Initialize state variable names and diagnostics
    variable_abbrev, variable_names, variable_units, index, pelagic_index, seaice_index, benthic_index = set_variable_info_bfm()

    # Allocate memory and give initial values (the argument list is mandatory with BFM)
    # init_var_bfm()

    # Set leading restart flag
    bfm_init = params_POMBFM.ihotst

    # Set the thickness of the vertical_layers-1 layers
    depth = np.zeros(vertical_layers)
    depth[:] = vertical_grid.vertical_spacing[:] * params_POMBFM.h

    # Allocate and initialize additional integration arrays
    d3state = np.zeros((num_boxes,num_d3_box_states))
    # try:
    #     INCLUDE_BEN
    # except NameError:
    #     INCLUDE_BEN = False
    # else:
    #     INCLUDE_BEN = True
    if INCLUDE_BEN:
        d2stateb_ben = np.zeros((num_boxes_xy,num_d2_box_states_ben))

    # Define initial conditions
    bfm_variables = set_initial_conditions()

    # init_save_bfm()

    # Initialize prior time step for leapfrog
    d3state = bfm_variables[0:num_boxes,:]
    d3stateb = d3state

    return d3state, d3stateb


def read_bfm_input():

    phyto_input = current_path + '/inputs/BFM17_BERM_INIT/init_prof_Pc_150m_bermuda_killworth.da'
    zoop_input  = current_path + '/inputs/BFM17_BERM_INIT/init_prof_Zc_150m_bermuda_killworth.da'
    poc_input   = current_path + '/inputs/BFM17_BERM_INIT/init_prof_POC_150m_bermuda_killworth.da'
    doc_input   = current_path + '/inputs/BFM17_BERM_INIT/init_prof_DOC_150m_bermuda_killworth.da'
    phos_input  = current_path + '/inputs/BFM17_BERM_INIT/init_prof_P_150m_bermuda_killworth.da'
    nit_input   = current_path + '/inputs/BFM17_BERM_INIT/init_prof_N_150m_bermuda_killworth.da'
    am_input    = current_path + '/inputs/BFM17_BERM_INIT/init_prof_Am_150m_bermuda_killworth.da'
    oxy_input   = current_path + '/inputs/BFM17_BERM_INIT/init_prof_Oxy_150m_bermuda_killworth.da'

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   PHYTOPLANKTON CARBON IC
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(phyto_input):
        p2c = np.fromfile(phyto_input)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   ZOOPLANKTON CARBON IC
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(zoop_input):
        z5c = np.fromfile(zoop_input)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   PARTICULATE ORGANIC CARBON IC
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(poc_input):
        r6c = np.fromfile(poc_input)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   DISOLVED ORGANIC CARBON IC
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(doc_input):
        r1c = np.fromfile(doc_input)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   PHOSPHATE IC
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(phos_input):
        n1p = np.fromfile(phos_input)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   NITRATE IC
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(nit_input):
        n3n = np.fromfile(nit_input)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   AMMONIUM IC
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(am_input):
        n4n = np.fromfile(am_input)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   OXYGEN IC
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(oxy_input):
        o2o = np.fromfile(oxy_input)

    return p2c, z5c, r6c, r1c, n1p, n3n, n4n, o2o


def set_initial_conditions():

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Local Variables
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    p_nRc = 0.0126
    p_pRc = 0.7862e-3
    p_sRc = 0.0118
    p_iRc = 1./25.

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Definition of BGC global variables
    #   IrrOPT in the equation of Steele and light
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    eir = np.zeros(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Definition of general pelagic state variables: Pelagic Fases
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    p2c, z5c, r6c, r1c, n1p, n3n, n4n, o2o = read_bfm_input()
    n4n[:] = 0. # From get_IC.f90 --> N4n(k) = 0.0

    o3c = o3c0 * np.ones(vertical_layers)
    o3h = o3h0 * np.ones(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Pelagic nutrients (mMol / m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    n5s = n5s0 * np.ones(vertical_layers)
    o4n = o4n0 * np.ones(vertical_layers)
    n6r = n6r0 * np.ones(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Pelagic detritus (respectively mg C/m3 mMol N/m3 mMol P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    r6n = r6c * p_nRc
    r6p = r6c * p_pRc
    r6s = r6c * p_sRc
    r2c = r2c0 * np.ones(vertical_layers)
    r3c = r3c0 * np.ones(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Dissolved organic matter
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    r1n = r1c * p_nRc * 0.5
    r1p = r1c * p_pRc * 0.5

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   State variables for phytoplankton model
    #   pelagic diatoms (respectively mg C/m3 mMol N/m3 mMol P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    d1cc = p1c0
    d1cn = d1cc * p_nRc
    d1cp = d1cc * p_pRc
    d1cs = d1cc * p_sRc
    d1ci = d1cc * p_iRc

    if calcphytoplankton[0]:
        p1c = d1cc * np.ones(vertical_layers)
        p1n = d1cn * np.ones(vertical_layers)
        p1p = d1cp * np.ones(vertical_layers)
        p1s = d1cs * np.ones(vertical_layers)
        p1l = d1ci * np.ones(vertical_layers)
    else:
        p1c = np.zeros(vertical_layers)
        p1n = np.zeros(vertical_layers)
        p1p = np.zeros(vertical_layers)
        p1s = np.zeros(vertical_layers)
        p1l = np.zeros(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Pelagic flagellates (respectively mg C/m3 mMol N/m3 mMol P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    d1cc = p2c0
    d1cn = d1cc * p_nRc
    d1cp = d1cc * p_pRc
    d1ci = d1cc * p_iRc * 0.5

    if calcphytoplankton[1]:
        d2cc = p2c
        p2n = d2cc * p_nRc
        p2p = d2cc * p_pRc
        p2l = d2cc * p_iRc * 0.5
    else:
        p2c = np.zeros(vertical_layers)
        p2n = np.zeros(vertical_layers)
        p2p = np.zeros(vertical_layers)
        p2l = np.zeros(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Picophytoplankton (respectively mg C/m3 mMol N/m3 mMol P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    d1cc = p3c0
    d1cn = d1cc * p_nRc
    d1cp = d1cc * p_pRc
    d1ci = d1cc * p_iRc * 0.5

    if calcphytoplankton[2]:
        p3c = d1cc * np.ones(vertical_layers)
        p3n = d1cn * np.ones(vertical_layers)
        p3p = d1cp * np.ones(vertical_layers)
        p3l = d1ci * np.ones(vertical_layers)
    else:
        p3c = np.zeros(vertical_layers)
        p3n = np.zeros(vertical_layers)
        p3p = np.zeros(vertical_layers)
        p3l = np.zeros(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Large phytoplankton (respectively mg C/m3 mMol N/m3 mMol P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    d1cc = p4c0
    d1cn = d1cc * p_nRc
    d1cp = d1cc * p_pRc
    d1ci = d1cc * p_iRc * 0.5

    if calcphytoplankton[3]:
        p4c = d1cc * np.ones(vertical_layers)
        p4n = d1cn * np.ones(vertical_layers)
        p4p = d1cp * np.ones(vertical_layers)
        p4l = d1ci * np.ones(vertical_layers)
    else:
        p4c = np.zeros(vertical_layers)
        p4n = np.zeros(vertical_layers)
        p4p = np.zeros(vertical_layers)
        p4l = np.zeros(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   State variables for mesozooplankton model
    #   Carnivorous mesozooplankton ( mg C/m3 )
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if calcmesozooplankton[0]:
        z3c = z3c0 * np.ones(vertical_layers)
        z3n = z3c * p_nRc
        z3p = z3c * p_pRc
    else:
        z3c = np.zeros(vertical_layers)
        z3n = np.zeros(vertical_layers)
        z3p = np.zeros(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Omnivorous mesozooplankton ( mg C/m3 )
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if calcmesozooplankton[1]:
        z4c = z4c0 * np.ones(vertical_layers)
        z4n = z4c * p_nRc
        z4p = z4c * p_pRc
    else:
        z4c = np.zeros(vertical_layers)
        z4n = np.zeros(vertical_layers)
        z4p = np.zeros(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   State variables for microzooplankton model
    #   Pelagic microzooplankton  (respectively mg C/m3 mMol N/m3 mMol P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if calcmicrozooplankton[0]:
        # z5c = z5c0 * np.ones(vertical_layers)
        z5n = z5c * p_nRc
        z5p = z5c * p_pRc
    else:
        z5c = np.zeros(vertical_layers)
        z5n = np.zeros(vertical_layers)
        z5p = np.zeros(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Heterotrophic flagellates (respectively mg C/m3 mMol N/m3 mMol P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if calcmicrozooplankton[1]:
        z6c = z6c0 * np.ones(vertical_layers)
        z6n = z6c * p_nRc
        z6p = z6c * p_pRc
    else:
        z6c = np.zeros(vertical_layers)
        z6n = np.zeros(vertical_layers)
        z6p = np.zeros(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   State variables for pelagic bacteria model B1
    #   Pelagic bacteria (respectively mg C/m3 mMol N/m3 mMol P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if calcpelbacteria:
        b1c = b1c0 * np.ones(vertical_layers)
        b1n = b1c * p_nRc
        b1p = b1c * p_pRc
    else:
        b1c = np.zeros(vertical_layers)
        b1n = np.zeros(vertical_layers)
        b1p = np.zeros(vertical_layers)


    # try:
    #     INCLUDE_BEN
    # except NameError:
    #     INCLUDE_BEN = False
    # else:
    #     INCLUDE_BEN = True
    if INCLUDE_BEN:

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   State variables for the benthic modules
        #   Zoobenthos
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        if calcbenorganisms[0]:
            y1c = y1c0 * np.ones(vertical_layers)
        else:
            y1c = np.zeros(vertical_layers)

        if calcbenorganisms[1]:
            y2c = y2c0 * np.ones(vertical_layers)
        else:
            y2c = np.zeros(vertical_layers)

        if calcbenorganisms[2]:
            y3c = y3c0 * np.ones(vertical_layers)
        else:
            y3c = np.zeros(vertical_layers)

        if calcbenorganisms[3]:
            y4c = y4c0 * np.ones(vertical_layers)
        else:
            y4c = np.zeros(vertical_layers)

        if calcbenorganisms[4]:
            y5c = y5c0 * np.ones(vertical_layers)
        else:
            y5c = np.zeros(vertical_layers)

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   Bacteria
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        if calcbenbacteria[0]:
            h1c = h1c0 * np.ones(vertical_layers)
        else:
            h1c = np.zeros(vertical_layers)
        if calcbenbacteria[1]:
            h2c = h2c0 * np.ones(vertical_layers)
        else:
            h2c = np.zeros(vertical_layers)

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   Benthic nutrients
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        k5s = 20.75 * np.ones(vertical_layers)
        k6r = k6r0 * np.ones(vertical_layers)
        k4n = k4n0 * np.ones(vertical_layers)
        k14n = k14n0 * np.ones(vertical_layers)
        k24n = k24n0 * np.ones(vertical_layers)
        k1p = k1p0 * np.ones(vertical_layers)
        k11p = k11p0 * np.ones(vertical_layers)
        k21p = k21p0 * np.ones(vertical_layers)
        k3n = k3n0 * np.ones(vertical_layers)

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   Benthic detritus (respectively mg C/m3 mMol N/m3 mMol P/m3)
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        q1c = q1c0 * np.ones(vertical_layers)
        q11c = q11c0 * np.ones(vertical_layers)

        # try:
        #     IMPFLUX
        # except NameError:
        #     IMPFLUX = False
        # else:
        #     IMPFLUX = True
        if IMPFLUX:
            q6c = 1.E9 * np.ones(vertical_layers)
            q6n = 1.E9 * np.ones(vertical_layers)
            q6p = 1.E9 * np.ones(vertical_layers)
            q6s = 1.E9 * np.ones(vertical_layers)
        else:
            q6c = q6c0 * np.ones(vertical_layers)
            q6n = q6n0 * np.ones(vertical_layers)
            q6p = q6p0 * np.ones(vertical_layers)
            q6s = q6s0 * np.ones(vertical_layers)

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   Gases
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        g2o = g2o0 * np.ones(vertical_layers)
        g3c = g3c0 * np.ones(vertical_layers)
        g4n = 37. * np.ones(vertical_layers)

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   Layers
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        d1m = d1m0 * np.ones(vertical_layers)
        d2m = d2m0 * np.ones(vertical_layers)
        d6m = d6m0 * np.ones(vertical_layers)
        d7m = d7m0 * np.ones(vertical_layers)
        d8m = d8m0 * np.ones(vertical_layers)
        d9m = d9m0 * np.ones(vertical_layers)

    # Fill matrix with bfm variable data
    bfm_variables = np.zeros((vertical_layers,50))

    bfm_variables[:,0] = o2o[:]
    bfm_variables[:,1] = n1p[:]
    bfm_variables[:,2] = n3n[:]
    bfm_variables[:,3] = n4n[:]
    bfm_variables[:,4] = o4n[:]
    bfm_variables[:,5] = n5s[:]
    bfm_variables[:,6] = n6r[:]

    bfm_variables[:,7] = b1c[:]
    bfm_variables[:,8] = b1n[:]
    bfm_variables[:,9] = b1p[:]
 
    bfm_variables[:,10] = p1c[:]
    bfm_variables[:,11] = p1n[:]
    bfm_variables[:,12] = p1p[:]
    bfm_variables[:,13] = p1l[:]
    bfm_variables[:,14] = p1s[:]
 
    bfm_variables[:,15] = p2c[:]
    bfm_variables[:,16] = p2n[:]
    bfm_variables[:,17] = p2p[:]
    bfm_variables[:,18] = p2l[:]

    bfm_variables[:,19] = p3c[:]
    bfm_variables[:,20] = p3n[:]
    bfm_variables[:,21] = p3p[:]
    bfm_variables[:,22] = p3l[:]

    bfm_variables[:,23] = p4c[:]
    bfm_variables[:,24] = p4n[:]
    bfm_variables[:,25] = p4p[:]
    bfm_variables[:,26] = p4l[:]

    bfm_variables[:,27] = z3c[:]
    bfm_variables[:,28] = z3n[:]
    bfm_variables[:,29] = z3p[:]

    bfm_variables[:,30] = z4c[:]
    bfm_variables[:,31] = z4n[:]
    bfm_variables[:,32] = z4p[:]
 
    bfm_variables[:,33] = z5c[:]
    bfm_variables[:,34] = z5n[:]
    bfm_variables[:,35] = z5p[:]
 
    bfm_variables[:,36] = z6c[:]
    bfm_variables[:,37] = z6n[:]
    bfm_variables[:,38] = z6p[:]

    bfm_variables[:,39] = r1c[:]
    bfm_variables[:,40] = r1n[:]
    bfm_variables[:,41] = r1p[:]
 
    bfm_variables[:,42] = r2c[:]
    bfm_variables[:,43] = r3c[:]

    bfm_variables[:,44] = r6c[:]
    bfm_variables[:,45] = r6n[:]
    bfm_variables[:,46] = r6p[:]
    bfm_variables[:,47] = r6s[:]

    bfm_variables[:,48] = o3c[:]
    bfm_variables[:,49] = o3h[:]

    return bfm_variables


# UPDATED VARIABLE NAMES    -    NAME_TYPE_CONSTITUENT
#
# n1p - n1p
# n3n - n3n
# n4n - n4n
# n5s - n5s
# n6r - n6r
# n7f - n7f
# o2o - o2o
# o3c - o3c
# o3h - o3h
# o4n - o4n
# p1c - p1c
# p1n - p1n
# p1p - p1p
# p1l - p1l
# p1s - p1s
# p1f - p1f
# p2c - p2c
# p2n - p2n
# p2p - p2p
# p2l - p2l
# p2f - p2f
# p3c - p3c
# p3n - p3n
# p3p - p3p
# p3l - p3l
# p3f - p3f
# p4c - p4c
# p4n - p4n
# p4p - p4p
# p4l - p4l
# p4f - p4f
# b1c - b1c
# b1n - b1n
# b1p - b1p
# z3c - z3c
# z3n - z3n
# z3p - z3p
# z4c - z4c
# z4n - z4n
# z4p - z4p
# z5c - z5c
# z5n - z5n
# z5p - z5p
# z6c - z6c
# z6n - z6n
# z6p - z6p
# r1c - r1c
# r1n - r1n
# r1p - r1p
# r1f - r1f
# r2c - r2c
# r3c - r3c
# r6c - r6c
# r6n - r6n
# r6p - r6p
# r6s - r6s
# r6f - r6f











