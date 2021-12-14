import numpy as np
from pom.constants import current_path, seconds_per_day, twice_the_timestep, vertical_layers
from inputs import params_POMBFM
from bfm.constants import num_d3_box_states, num_d2_box_states_ben
from bfm.variable_info import set_variable_info_bfm
from os import path
from inputs.namelist_input_data import disOxygen_IO_O0, phospate_IO_P0, nitrate_IO_N0, ammonium_IO_N0, silicate_IO_Si0, reductEquiv_IO_R0, disInorgCarbon_IO_C0, o3h0, o4n0, diatoms_LO_C0, nanoflagellates_LO_C0, picophyto_LO_C0, largephyto_LO_C0, \
    carnivMesozoo_LO_C0, omnivMesozoo_LO_C0, microzoo_LO_C0, z6c0, pelBacteria_LO_C0, labileDOM_NO_C0, semilabileDOC_NO_C0, semirefractDOC_NO_C0, particOrganDetritus_NO_C0, \
    y1c0, y2c0, y3c0, y4c0, y5c0, h1c0, h2c0, k1p0, k11p0, k21p0, k4n0, k14n0, k24n0, k3n0, k5s0, k6r0, \
    d1m0, d2m0, d6m0, d7m0, d8m0, d9m0, q6c0, q6n0, q6p0, q6s0, q1c0, q11c0, g2o0, g3c0, g13c0, g23c0, g3h0, g13h0, g23h0, \
    calcphytoplankton, calcmesozooplankton, calcmicrozooplankton, calcpelbacteria, calcbenorganisms, calcbenbacteria

# Unresolved reference: INCLUDE_BEN

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

    try:
        INCLUDE_BEN
    except NameError:
        INCLUDE_BEN = False
    else:
        INCLUDE_BEN = True
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
    d3stateb = np.zeros((num_d3_box_states,num_boxes))
    try:
        INCLUDE_BEN
    except NameError:
        INCLUDE_BEN = False
    else:
        INCLUDE_BEN = True
    if INCLUDE_BEN:
        d2stateb_ben = np.zeros((num_d2_box_states_ben,num_boxes_xy))

    # Define initial conditions
    [disOxygen_IO_O, phospate_IO_P, nitrate_IO_N, ammonium_IO_N, o4n, silicate_IO_Si, reductEquiv_IO_R, pelBacteria_LO_C, pelBacteria_LO_N, pelBacteria_LO_P,
     diatoms_LO_C, diatoms_LO_N, diatoms_LO_P, diatoms_LO_Chl, diatoms_LO_Si, nanoflagellates_LO_C, nanoflagellates_LO_N, nanoflagellates_LO_P, nanoflagellates_LO_Chl, picophyto_LO_C, picophyto_LO_N, picophyto_LO_P, picophyto_LO_Chl, largephyto_LO_C, largephyto_LO_N, largephyto_LO_P, largephyto_LO_Chl,
     carnivMesozoo_LO_C, carnivMesozoo_LO_N, carnivMesozoo_LO_P, omnivMesozoo_LO_C, omnivMesozoo_LO_N, omnivMesozoo_LO_P, microzoo_LO_C, microzoo_LO_N, microzoo_LO_P, z6c, z6n, z6p,
     labileDOM_NO_C, labileDOM_NO_N, labileDOM_NO_P, semilabileDOC_NO_C, semirefractDOC_NO_C, particOrganDetritus_NO_C, particOrganDetritus_NO_N, particOrganDetritus_NO_P, particOrganDetritus_NO_Si, disInorgCarbon_IO_C, o3h] = set_initial_conditions()

    # init_save_bfm()

    # Initialize prior time step for leapfrog
    d3state = np.zeros((num_d3_box_states,num_boxes))   # this isn't the real matrix, just a placeholder until api_bfm is translated
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
        nanoflagellates_LO_C = np.fromfile(phyto_input)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   ZOOPLANKTON CARBON IC
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(zoop_input):
        microzoo_LO_C = np.fromfile(zoop_input)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   PARTICULATE ORGANIC CARBON IC
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(poc_input):
        particOrganDetritus_NO_C = np.fromfile(poc_input)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   DISOLVED ORGANIC CARBON IC
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(doc_input):
        labileDOM_NO_C = np.fromfile(doc_input)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   PHOSPHATE IC
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(phos_input):
        phospate_IO_P = np.fromfile(phos_input)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   NITRATE IC
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(nit_input):
        nitrate_IO_N = np.fromfile(nit_input)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   AMMONIUM IC
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(am_input):
        ammonium_IO_N = np.fromfile(am_input)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   OXYGEN IC
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if path.exists(oxy_input):
        disOxygen_IO_O = np.fromfile(oxy_input)

    return nanoflagellates_LO_C, microzoo_LO_C, particOrganDetritus_NO_C, labileDOM_NO_C, phospate_IO_P, nitrate_IO_N, ammonium_IO_N, disOxygen_IO_O


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
    nanoflagellates_LO_C, microzoo_LO_C, particOrganDetritus_NO_C, labileDOM_NO_C, phospate_IO_P, nitrate_IO_N, ammonium_IO_N, disOxygen_IO_O = read_bfm_input()

    disInorgCarbon_IO_C = disInorgCarbon_IO_C0 * np.ones(vertical_layers)
    o3h = o3h0 * np.ones(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Pelagic nutrients (mMol / m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    silicate_IO_Si = silicate_IO_Si0 * np.ones(vertical_layers)
    o4n = o4n0 * np.ones(vertical_layers)
    reductEquiv_IO_R = reductEquiv_IO_R0 * np.ones(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Pelagic detritus (respectively mg C/m3 mMol N/m3 mMol P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    particOrganDetritus_NO_N = particOrganDetritus_NO_C * p_nRc
    particOrganDetritus_NO_P = particOrganDetritus_NO_C * p_pRc
    particOrganDetritus_NO_Si = particOrganDetritus_NO_C * p_sRc
    semilabileDOC_NO_C = semilabileDOC_NO_C0 * np.ones(vertical_layers)
    semirefractDOC_NO_C = semirefractDOC_NO_C0 * np.ones(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Dissolved organic matter
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    labileDOM_NO_N = labileDOM_NO_C * p_nRc * 0.5
    labileDOM_NO_P = labileDOM_NO_C * p_pRc * 0.5

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   State variables for phytoplankton model
    #   pelagic diatoms (respectively mg C/m3 mMol N/m3 mMol P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    d1cc = diatoms_LO_C0
    d1cn = d1cc * p_nRc
    d1cp = d1cc * p_pRc
    d1cs = d1cc * p_sRc
    d1ci = d1cc * p_iRc

    if calcphytoplankton[0]:
        diatoms_LO_C = d1cc * np.ones(vertical_layers)
        diatoms_LO_N = d1cn * np.ones(vertical_layers)
        diatoms_LO_P = d1cp * np.ones(vertical_layers)
        diatoms_LO_Si = d1cs * np.ones(vertical_layers)
        diatoms_LO_Chl = d1ci * np.ones(vertical_layers)
    else:
        diatoms_LO_C = np.zeros(vertical_layers)
        diatoms_LO_N = np.zeros(vertical_layers)
        diatoms_LO_P = np.zeros(vertical_layers)
        diatoms_LO_Si = np.zeros(vertical_layers)
        diatoms_LO_Chl = np.zeros(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Pelagic flagellates (respectively mg C/m3 mMol N/m3 mMol P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    d1cc = nanoflagellates_LO_C0
    d1cn = d1cc * p_nRc
    d1cp = d1cc * p_pRc
    d1ci = d1cc * p_iRc * 0.5

    if calcphytoplankton[1]:
        d2cc = nanoflagellates_LO_C
        nanoflagellates_LO_N = d2cc * p_nRc
        nanoflagellates_LO_P = d2cc * p_pRc
        nanoflagellates_LO_Chl = d2cc * p_iRc * 0.5
    else:
        nanoflagellates_LO_C = np.zeros(vertical_layers)
        nanoflagellates_LO_N = np.zeros(vertical_layers)
        nanoflagellates_LO_P = np.zeros(vertical_layers)
        nanoflagellates_LO_Chl = np.zeros(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Picophytoplankton (respectively mg C/m3 mMol N/m3 mMol P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    d1cc = picophyto_LO_C0
    d1cn = d1cc * p_nRc
    d1cp = d1cc * p_pRc
    d1ci = d1cc * p_iRc * 0.5

    if calcphytoplankton[2]:
        picophyto_LO_C = d1cc * np.ones(vertical_layers)
        picophyto_LO_N = d1cn * np.ones(vertical_layers)
        picophyto_LO_P = d1cp * np.ones(vertical_layers)
        picophyto_LO_Chl = d1ci * np.ones(vertical_layers)
    else:
        picophyto_LO_C = np.zeros(vertical_layers)
        picophyto_LO_N = np.zeros(vertical_layers)
        picophyto_LO_P = np.zeros(vertical_layers)
        picophyto_LO_Chl = np.zeros(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Large phytoplankton (respectively mg C/m3 mMol N/m3 mMol P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    d1cc = largephyto_LO_C0
    d1cn = d1cc * p_nRc
    d1cp = d1cc * p_pRc
    d1ci = d1cc * p_iRc * 0.5

    if calcphytoplankton[3]:
        largephyto_LO_C = d1cc * np.ones(vertical_layers)
        largephyto_LO_N = d1cn * np.ones(vertical_layers)
        largephyto_LO_P = d1cp * np.ones(vertical_layers)
        largephyto_LO_Chl = d1ci * np.ones(vertical_layers)
    else:
        largephyto_LO_C = np.zeros(vertical_layers)
        largephyto_LO_N = np.zeros(vertical_layers)
        largephyto_LO_P = np.zeros(vertical_layers)
        largephyto_LO_Chl = np.zeros(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   State variables for mesozooplankton model
    #   Carnivorous mesozooplankton ( mg C/m3 )
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if calcmesozooplankton[0]:
        carnivMesozoo_LO_C = carnivMesozoo_LO_C0 * np.ones(vertical_layers)
        carnivMesozoo_LO_N = carnivMesozoo_LO_C * p_nRc
        carnivMesozoo_LO_P = carnivMesozoo_LO_C * p_pRc
    else:
        carnivMesozoo_LO_C = np.zeros(vertical_layers)
        carnivMesozoo_LO_N = np.zeros(vertical_layers)
        carnivMesozoo_LO_P = np.zeros(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   Omnivorous mesozooplankton ( mg C/m3 )
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if calcmesozooplankton[1]:
        omnivMesozoo_LO_C = omnivMesozoo_LO_C0 * np.ones(vertical_layers)
        omnivMesozoo_LO_N = omnivMesozoo_LO_C * p_nRc
        omnivMesozoo_LO_P = omnivMesozoo_LO_C * p_pRc
    else:
        omnivMesozoo_LO_C = np.zeros(vertical_layers)
        omnivMesozoo_LO_N = np.zeros(vertical_layers)
        omnivMesozoo_LO_P = np.zeros(vertical_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   State variables for microzooplankton model
    #   Pelagic microzooplankton  (respectively mg C/m3 mMol N/m3 mMol P/m3)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if calcmicrozooplankton[0]:
        microzoo_LO_C = microzoo_LO_C0 * np.ones(vertical_layers)
        microzoo_LO_N = microzoo_LO_C * p_nRc
        microzoo_LO_P = microzoo_LO_C * p_pRc
    else:
        microzoo_LO_C = np.zeros(vertical_layers)
        microzoo_LO_N = np.zeros(vertical_layers)
        microzoo_LO_P = np.zeros(vertical_layers)

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
    if calcpelbacteria[0]:
        pelBacteria_LO_C = pelBacteria_LO_C0 * np.ones(vertical_layers)
        pelBacteria_LO_N = pelBacteria_LO_C * p_nRc
        pelBacteria_LO_P = pelBacteria_LO_C * p_pRc
    else:
        pelBacteria_LO_C = np.zeros(vertical_layers)
        pelBacteria_LO_N = np.zeros(vertical_layers)
        pelBacteria_LO_P = np.zeros(vertical_layers)


    try:
        INCLUDE_BEN
    except NameError:
        INCLUDE_BEN = False
    else:
        INCLUDE_BEN = True
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

        try:
            IMPFLUX
        except NameError:
            IMPFLUX = False
        else:
            IMPFLUX = True
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




    return [disOxygen_IO_O, phospate_IO_P, nitrate_IO_N, ammonium_IO_N, o4n, silicate_IO_Si, reductEquiv_IO_R, pelBacteria_LO_C, pelBacteria_LO_N, pelBacteria_LO_P,
            diatoms_LO_C, diatoms_LO_N, diatoms_LO_P, diatoms_LO_Chl, diatoms_LO_Si, nanoflagellates_LO_C, nanoflagellates_LO_N, nanoflagellates_LO_P, nanoflagellates_LO_Chl, picophyto_LO_C, picophyto_LO_N, picophyto_LO_P, picophyto_LO_Chl, largephyto_LO_C, largephyto_LO_N, largephyto_LO_P, largephyto_LO_Chl,
            carnivMesozoo_LO_C, carnivMesozoo_LO_N, carnivMesozoo_LO_P, omnivMesozoo_LO_C, omnivMesozoo_LO_N, omnivMesozoo_LO_P, microzoo_LO_C, microzoo_LO_N, microzoo_LO_P, z6c, z6n, z6p,
            labileDOM_NO_C, labileDOM_NO_N, labileDOM_NO_P, semilabileDOC_NO_C, semirefractDOC_NO_C, particOrganDetritus_NO_C, particOrganDetritus_NO_N, particOrganDetritus_NO_P, particOrganDetritus_NO_Si, disInorgCarbon_IO_C, o3h]


# UPDATED VARIABLE NAMES    -    NAME_TYPE_CONSTITUENT
#
# n1p - phospate_IO_P
# n3n - nitrate_IO_N
# n4n - ammonium_IO_N
# n5s - silicate_IO_Si
# n6r - reductEquiv_IO_R
# n7f - reductEquiv_IO_Fe
# o2o - disOxygen_IO_O
# o3c - disInorgCarbon_IO_C
# o3h - totalAlkalinity_IO
# p1c - diatoms_LO_C
# p1n - diatoms_LO_N
# p1p - diatoms_LO_P
# p1l - diatoms_LO_Chl
# p1s - diatoms_LO_Si
# p1f - diatoms_LO_Fe
# p2c - nanoflagellates_LO_C
# p2n - nanoflagellates_LO_N
# p2p - nanoflagellates_LO_P
# p2l - nanoflagellates_LO_Chl
# p2f - nanoflagellates_LO_Fe
# p3c - picophyto_LO_C
# p3n - picophyto_LO_N
# p3p - picophyto_LO_P
# p3l - picophyto_LO_L
# p3f - picophyto_LO_Fe
# p4c - largephyto_LO_C
# p4n - largephyto_LO_N
# p4p - largephyto_LO_P
# p4l - largephyto_LO_Chl
# p4f - largephyto_LO_Fe
# b1c - pelBacteria_LO_C
# b1n - pelBacteria_LO_N
# b1p - pelBacteria_LO_P
# z3c - carnivMesozoo_LO_C
# z3n - carnivMesozoo_LO_N
# z3p - carnivMesozoo_LO_P
# z4c - omnivMesozoo_LO_C
# z4n - omnivMesozoo_LO_N
# z4p - omnivMesozoo_LO_P
# z5c - microzoo_LO_C
# z5n - microzoo_LO_N
# z5p - microzoo_LO_P
# z6c - heteroFlagellates_LO_C
# z6n - heteroFlagellates_LO_N
# z6p - heteroFlagellates_LO_P
# r1c - labileDOM_NO_C
# r1n - labileDOM_NO_N
# r1p - labileDOM_NO_P
# r1f - labileDOM_NO_Fe
# r2c - semilabileDOC_NO_C
# r3c - semirefractDOC_NO_C
# r6c - particOrganDetritus_NO_C
# r6n - particOrganDetritus_NO_N
# r6p - particOrganDetritus_NO_P
# r6s - particOrganDetritus_NO_Si
# r6f - particOrganDetritus_NO_Fe











