import numpy as np
from pom.constants import seconds_per_day, twice_the_timestep, vertical_layers
from inputs import params_POMBFM
from bfm.constants import num_d3_box_states, num_d2_box_states_ben
from bfm.variable_info import set_variable_info_bfm


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
    # set_initial_conditions()
    # init_save_bfm()

    # Initialize prior time step for leapfrog
    # d3stateb = d3state









