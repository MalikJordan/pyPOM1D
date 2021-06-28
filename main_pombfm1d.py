from cppdefs import *
import numpy as np
from pom.modules import vertical_layers, DAYI
from inputs import params_POMBFM
from pom.calculations import *
from pom.initialize_variables import *
from pom.create_profiles import *

# pyPOM1D DIRECTORY, USED FOR READING INPUTS (TO BE CHANGED BY USER)
current_path = '/Users/malikjordan/Desktop/pyPOM1D'

earth_angular_velocity = 7.29E-5                                                        # OMEGA

# VARIABLE NAMES (FORTRAN --> PYTHON)
# Z = vertical_coordinates
# ZZ = vertical_coordinates_staggered
# DZ = vertical_spacing
# DZZ = vertical_spacing_staggered
# DZR = vertical_spacing_reciprocal
# T = temperature
# TF = temperature_forward
# TB = temperature_backward
# TSTAR = interpolated_temperature
# S = salinity
# SF = salinity_forward
# SB = salinity_backward
# SSTAR = interpolated_salinity
# RHO = density_profile
# U = velocity_zonal
# UF = velocity_zonal_forward
# UB = velocity_zonal_backward
# V = velocity_meridional
# VF = velocity_meridional_forward
# VB = velocity_meridional_backward
# Q2 = kinetic_energy
# Q2F = kinetic_energy_forward
# Q2B = kinetic_energy_backward
# Q2L = kinetic_energy_times_length
# Q2LF = kinetic_energy_times_length_forward
# Q2LB = kinetic_energy_times_length_backward
# L = length_scale
# KM = diffusion_coefficient_momentum
# KH = diffusion_coefficient_tracers
# KQ = diffusion_coefficient_kinetic_energy
# WTADV = temperature_lateral_advection
# WSADV = salinity_lateral_advection
# WUSURF = wind_stress_zonal
# WVSURF = wind_stress_meridional
# WUBOT = bottom_stress_zonal
# WVBOT = bottom_stress_meridional
# WTSURF = surface_heat_flux
# SWRAD = shortwave_radiation
# WSSURF = surface_salinity_flux
# TSURF = surface_temperature
# SSURF = surface_salinity

# GENERAL INITIALIZATION
vertical_spacing_reciprocal = np.zeros(vertical_layers)                     # DZR
interpolated_temperature = np.zeros(vertical_layers)                        # TSTAR, from ModuleForcing.F90, line 78
interpolated_salinity = np.zeros(vertical_layers)                           # SSTAR, from ModuleForcing.F90, line 78
temperature_lateral_advection = np.zeros(vertical_layers)                   # WTADV
salinity_lateral_advection = np.zeros(vertical_layers)                      # WSADV


# DEFINE VERTICAL COORDINATE SYSTEM
vertical_coordinates, vertical_coordinates_staggered, \
    vertical_spacing, vertical_spacing_staggered        = create_vertical_coordinate_system(vertical_layers, params_POMBFM.kl1, params_POMBFM.kl2)

# VERTICAL SPACING RECIROCAL (DZR)
vertical_spacing[vertical_layers-1] = 1.E-6  # Small value to avoid division by zero for vertical_spacing_reciprocal
vertical_spacing_reciprocal[:] = 1. / vertical_spacing[:]
vertical_spacing[vertical_layers-1] = 0.

# CORIOLIS PARAMETER
coriolis_parameter = 2. * earth_angular_velocity * np.sin(alat * 2. * np.pi / 360.)        # COR

# TWICE THE TIMESTEP
twice_the_timestep = 2. * params_POMBFM.dti                                                # DT2

# ITERATIONS NEEDED TO CARRY OUT AN"IDAYS" SIMULATION
iterations_needed = idays * seconds_per_day / params_POMBFM.dti                            # iend

# READ  T&S INITIAL CONDITIONS (IHOTST=0) OR RESTART FILE (IHOTST=1)
if params_POMBFM.ihotst == 0:
    time0 = 0.
    temperature, temperature_backward, salinity, salinity_backward = get_temperature_and_salinity_initial_coditions()
elif params_POMBFM.ihotst == 1:
    # get_rst()
    pass

try:
    POM_only
except NameError:
    POM_only = False
else:
    POM_only = True
if not POM_only:
    # INITIALIZATION OF BFM
    # pom_ini_bfm_1d()
    pass

# # BEGIN THE TIME MARCH
# print('ICOUNT before time march loop = ')
# # print(icountf)


# BEGIN THE TIME MARCH
for i in range(0, int(iterations_needed)):
    time = time0 + (dti * i * DAYI)

    # TURBULENCE CLOSURE
    kinetic_energy_forward[:] = kinetic_energy_backward[:]
    kinetic_energy_times_length_forward[:] = kinetic_energy_times_length_backward[:]

    PROFQ()
    # DEFINE ALL FORCINGS
    forcing_manager(i)

    # T&S COMPUTATION
    if params_POMBFM.idiagn == 0:
        # PROGNOSTIC MODE
        # T&S FULLY COMPUTED BY MODEL

        if params_POMBFM.trt != 0:
            for j in range(0, vertical_layers):
                if (-vertical_coordinates_staggered[j] * params_POMBFM.h) >= params_POMBFM.upperh:
                    temperature_lateral_advection[j] = (interpolated_temperature[j] - temperature[j]) / (params_POMBFM.trt * seconds_per_day)

        if params_POMBFM.srt != 0:
            for j in range(0, vertical_layers):
                if (-vertical_coordinates_staggered[j] * params_POMBFM.h) >= params_POMBFM.upperh:
                    salinity_lateral_advection[j] = (interpolated_salinity[j] - salinity[j]) / (params_POMBFM.srt * seconds_per_day)

        # COMPUTE SURFACE SALINITY FLUX
        surface_salinity_flux = -(surface_salinity - salinity[0]) * params_POMBFM.srt / seconds_per_day

        # COMPUTE TEMPREATURE
        temperature_forward[:] = temperature_backward[:] + (temperature_lateral_advection[:] * twice_the_timestep)
        calculate_vertical_temperature_and_salinity_profiles()

        # CALCULATE SALINITY
        salinity_forward[:] = salinity_backward[:] + (salinity_lateral_advection[:] * twice_the_timestep)
        calculate_vertical_temperature_and_salinity_profiles()

        # MIXING THE TIMESTEP (ASSELIN)
        temperature[:] = temperature[:] + 0.5 * params_POMBFM.smoth * (temperature_forward[:] + temperature_backward[:] - 2. * temperature_forward[:])
        salinity[:] = salinity[:] + 0.5 * params_POMBFM.smoth * (salinity_forward[:] + salinity_backward[:] - 2. * salinity[:])

    # COMPUTE VELOCITY
    velocity_zonal_forward[:] = velocity_zonal_backward[:] + twice_the_timestep * coriolis_parameter * velocity_meridional[:]
    calculate_vertical_zonal_velocity_profile(vertical_spacing, vertical_spacing_staggered, wind_stress_zonal,
                                              diffusion_coefficient_momentum, velocity_zonal_forward)

    for j in range(0, vertical_layers):
        velocity_meridional_forward[j] = velocity_meridional_backward[j] - twice_the_timestep * coriolis_parameter * velocity_zonal[j]
    calculate_vertical_meridional_velocity_profile(vertical_spacing, vertical_spacing_staggered, wind_stress_meridional,
                                                   diffusion_coefficient_momentum, velocity_meridional_forward)

    # MIX TIME STEL (ASSELIN FILTER)
    for j in range(0, vertical_layers):
        kinetic_energy[j] = kinetic_energy[j] + 0.5 * params_POMBFM.smoth * (kinetic_energy_forward[j] + kinetic_energy_backward[j] - 2. * kinetic_energy[j])
        kinetic_energy_times_length[j] = kinetic_energy_times_length[j] * params_POMBFM.smoth * (kinetic_energy_times_length_forward[j] + kinetic_energy_times_length_backward[j] - 2. * kinetic_energy_times_length[j])

        velocity_zonal[j] = velocity_zonal[j] + 0.5 * params_POMBFM.smoth * (velocity_zonal_forward[j] + velocity_zonal_backward[j] - 2. * velocity_zonal[j])
        velocity_meridional[j] = velocity_meridional[j] + 0.5 * params_POMBFM.smoth * (velocity_meridional_forward[j] + velocity_meridional_backward[j] - 2. * velocity_meridional[j])

    # RESTORE TIME SEQUENCE
    kinetic_energy_backward[:] = kinetic_energy[:]
    kinetic_energy[:] = kinetic_energy_forward[:]
    kinetic_energy_times_length_backward[:] = kinetic_energy_times_length[:]
    kinetic_energy_times_length[:] = kinetic_energy_times_length_forward[:]

    velocity_zonal_backward[:] = velocity_zonal[:]
    velocity_zonal[:] = velocity_zonal_forward[:]
    velocity_meridional_backward[:] = velocity_meridional[:]
    velocity_meridional[:] = velocity_meridional_forward[:]

    temperature_backward[:] = temperature[:]
    temperature[:] = temperature_forward[:]
    salinity_backward[:] = salinity[:]
    salinity[:] = salinity_forward[:]

    # UPDATE DENSITY
    create_vertical_diffusivity_profile(temperature, salinity, vertical_spacing_staggered, dti, vertical_layers)

try:
    POM_only
except NameError:
    POM_only = False
else:
    POM_only = True
if not POM_only:
    # INITIALIZATION OF BFM
    # pom_ini_bfm_1d()
    pass

# WRITING OF RESTART
print(time)
print(velocity_zonal, velocity_zonal_backward, velocity_meridional, velocity_meridional_backward)
print(temperature, temperature_backward, salinity, salinity_backward)
print(kinetic_energy, kinetic_energy_backward, kinetic_energy_times_length, kinetic_energy_times_length_backward)
print(diffusion_coefficient_tracers, diffusion_coefficient_momentum, diffusion_coefficient_kinetic_energy)
print(length_scale)
print(bottom_stress_zonal, bottom_stress_meridional)
print(density_profile)

#BFM RESTART
try:
    POM_only
except NameError:
    POM_only = False
else:
    POM_only = True
if not POM_only:
    # INITIALIZATION OF BFM
    # restart_BFM_inPOM()
    pass

print('Main done')
