from cppdefs import *
import numpy as np
from pom.modules import forcing_manager
from inputs import params_POMBFM
from pom.calculations import calculate_vertical_density_profile, create_vertical_coordinate_system
from pom.initialize_variables import get_temperature_and_salinity_initial_coditions
from pom.create_profiles import create_kinetic_energy_profile, create_vertical_diffusivity_profile, \
    calculate_vertical_temperature_and_salinity_profiles, calculate_vertical_zonal_velocity_profile, calculate_vertical_meridional_velocity_profile
from pom.data_classes import DiffusionCoefficients, ForcingManagerCounters, LeapFrogTimeLevels, MonthlyForcingData, Stresses, TemperatureSalinityData, VelocityData
from pom.pom_constants import earth_angular_velocity, DAYI, water_specific_heat_times_density, vertical_layers, seconds_per_day

np.set_printoptions(precision=16)
# pyPOM1D DIRECTORY, USED FOR READING INPUTS (TO BE CHANGED BY USER)
# current_path = '/Users/malikjordan/Desktop/pyPOM1D'

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

# earth_angular_velocity = 7.29E-5                                                        # OMEGA
# vertical_layers = 151
# DAYI = 1. / seconds_per_day
# water_specific_heat_times_density = 4.187E6
# GENERAL INITIALIZATION
length_scale = np.ones(vertical_layers)
length_scale[0] = 0.
length_scale[vertical_layers-1] = 0.
diffusion = DiffusionCoefficients()
kinetic_energy = LeapFrogTimeLevels(1.E-07 * np.ones(vertical_layers),1.E-07 * np.ones(vertical_layers),1.E-07 * np.ones(vertical_layers))
kinetic_energy_times_length = LeapFrogTimeLevels(1.E-07 * np.ones(vertical_layers),1.E-07 * np.ones(vertical_layers),1.E-07 * np.ones(vertical_layers))
velocity = VelocityData()
wind_stress = Stresses()
bottom_stress = Stresses()
temperature = TemperatureSalinityData()
salinity = TemperatureSalinityData()

# DEFINE VERTICAL COORDINATE SYSTEM
vertical_grid = create_vertical_coordinate_system(params_POMBFM.kl1, params_POMBFM.kl2)
vertical_grid.length_scale = length_scale
# CORIOLIS PARAMETER
coriolis_parameter = 2. * earth_angular_velocity * np.sin(params_POMBFM.alat * 2. * np.pi / 360.)        # COR

# TWICE THE TIMESTEP
twice_the_timestep = 2. * params_POMBFM.dti                                                              # DT2

# ITERATIONS NEEDED TO CARRY OUT AN "IDAYS" SIMULATION
iterations_needed = params_POMBFM.idays * seconds_per_day / params_POMBFM.dti                            # iend

# READ  T&S INITIAL CONDITIONS (IHOTST=0) OR RESTART FILE (IHOTST=1)
if params_POMBFM.ihotst == 0:
    time0 = 0.
    temperature.current, temperature.backward, salinity.current, salinity.backward = get_temperature_and_salinity_initial_coditions()
    vertical_density_profile = calculate_vertical_density_profile(temperature, salinity, vertical_grid)
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
counters = ForcingManagerCounters()
month1_data = MonthlyForcingData()
month2_data = MonthlyForcingData()
# for i in range(0, int(iterations_needed)):
for i in range(0, 10000):
    time = time0 + (params_POMBFM.dti * i * DAYI)

    # TURBULENCE CLOSURE
    kinetic_energy.forward[:] = kinetic_energy.backward[:]
    kinetic_energy_times_length.forward[:] = kinetic_energy_times_length.backward[:]

    # print(kinetic_energy_times_length.backward[1])
    # print(kinetic_energy_times_length.current[1])
    # print(kinetic_energy_times_length.forward[1])
    # print()
    # print(velocity.zonal_current[1])
    # print()
    kinetic_energy, kinetic_energy_times_length, diffusion, vertical_grid = create_kinetic_energy_profile(vertical_grid, diffusion, temperature, salinity, vertical_density_profile, velocity,
                                                                                                          kinetic_energy, kinetic_energy_times_length, wind_stress, bottom_stress)
    # DEFINE ALL FORCINGS
    temperature.forward, temperature.interpolated, salinity.forward, salinity.interpolated, shortwave_radiation, \
        temperature.surface_flux, wind_stress, month1_data, month2_data, counters = forcing_manager(i,counters,month1_data,month2_data)

    # T&S COMPUTATION
    if params_POMBFM.idiagn == 0:
        # PROGNOSTIC MODE
        # T&S FULLY COMPUTED BY MODEL
        temperature.surface_value = temperature.forward[0]
        salinity.surface_value = salinity.forward[0]

        if params_POMBFM.trt != 0:
            for j in range(0, vertical_layers):
                if (-vertical_grid.vertical_coordinates_staggered[j] * params_POMBFM.h) >= params_POMBFM.upperh:
                    temperature.lateral_advection[j] = (temperature.interpolated[j] - temperature.current[j]) / (params_POMBFM.trt * seconds_per_day)

        if params_POMBFM.srt != 0:
            for j in range(0, vertical_layers):
                if (-vertical_grid.vertical_coordinates_staggered[j] * params_POMBFM.h) >= params_POMBFM.upperh:
                    salinity.lateral_advection[j] = (salinity.interpolated[j] - salinity.current[j]) / (params_POMBFM.srt * seconds_per_day)

        # COMPUTE SURFACE SALINITY FLUX
        salinity.surface_flux = -(salinity.surface_value - salinity.current[0]) * params_POMBFM.srt / seconds_per_day

        # COMPUTE TEMPREATURE
        temperature.forward[:] = temperature.backward[:] + (temperature.lateral_advection[:] * twice_the_timestep)
        calculate_vertical_temperature_and_salinity_profiles(vertical_grid, diffusion, temperature, shortwave_radiation, params_POMBFM.nbct)

        # CALCULATE SALINITY
        salinity.forward[:] = salinity.backward[:] + (salinity.lateral_advection[:] * twice_the_timestep)
        calculate_vertical_temperature_and_salinity_profiles(vertical_grid, diffusion, salinity, shortwave_radiation, params_POMBFM.nbcs)

        # MIXING THE TIMESTEP (ASSELIN)
        temperature.current[:] = temperature.current[:] + 0.5 * params_POMBFM.smoth * (temperature.forward[:] + temperature.backward[:] - 2. * temperature.current[:])
        salinity.current[:] = salinity.current[:] + 0.5 * params_POMBFM.smoth * (salinity.forward[:] + salinity.backward[:] - 2. * salinity.current[:])

    # print(diffusion.momentum[1])
    # COMPUTE VELOCITY
    # print(velocity.zonal_forward[1])
    # print()
    velocity.zonal_forward[:] = velocity.zonal_backward[:] + twice_the_timestep * coriolis_parameter * velocity.meridional_current[:]
    # print(velocity.zonal_forward[1])
    # print()
    velocity, bottom_stress = calculate_vertical_zonal_velocity_profile(vertical_grid, wind_stress, bottom_stress, diffusion, velocity)
    # print(velocity.zonal_forward[1])
    # print('-------------------------------------------------')

    velocity.meridional_forward[:] = velocity.meridional_backward[:] - twice_the_timestep * coriolis_parameter * velocity.zonal_current[:]
    velocity, bottom_stress = calculate_vertical_meridional_velocity_profile(vertical_grid, wind_stress, bottom_stress, diffusion, velocity)

    # MIX TIME STEL (ASSELIN FILTER)
    # print(kinetic_energy_times_length.backward[1])
    # print(kinetic_energy_times_length.current[1])
    # print(kinetic_energy_times_length.forward[1])
    # print('-------------------------------------------------')
    kinetic_energy.current[:] = kinetic_energy.current[:] + 0.5 * params_POMBFM.smoth * (kinetic_energy.forward[:] + kinetic_energy.backward[:] - 2. * kinetic_energy.current[:])
    kinetic_energy_times_length.current[:] = kinetic_energy_times_length.current[:] + 0.5 * params_POMBFM.smoth * (kinetic_energy_times_length.forward[:] + kinetic_energy_times_length.backward[:] - 2. * kinetic_energy_times_length.current[:])

    velocity.zonal_current[:] = velocity.zonal_current[:] + 0.5 * params_POMBFM.smoth * (velocity.zonal_forward[:] + velocity.zonal_backward[:] - 2. * velocity.zonal_current[:])
    velocity.meridional_current[:] = velocity.meridional_current[:] + 0.5 * params_POMBFM.smoth * (velocity.meridional_forward[:] + velocity.meridional_backward[:] - 2. * velocity.meridional_current[:])

    # RESTORE TIME SEQUENCE
    kinetic_energy.backward[:] = kinetic_energy.current[:]
    kinetic_energy.current[:] = kinetic_energy.forward[:]
    kinetic_energy_times_length.backward[:] = kinetic_energy_times_length.current[:]
    kinetic_energy_times_length.current[:] = kinetic_energy_times_length.forward[:]

    # print(kinetic_energy.backward[2])
    velocity.zonal_backward[:] = velocity.zonal_current[:]
    velocity.zonal_current[:] = velocity.zonal_forward[:]
    velocity.meridional_backward[:] = velocity.meridional_current[:]
    velocity.meridional_current[:] = velocity.meridional_forward[:]
    # print(velocity.zonal_current[1])
    # print('-------------------------------------------------')

    temperature.backward[:] = temperature.current[:]
    temperature.current[:] = temperature.forward[:]
    salinity.backward[:] = salinity.current[:]
    salinity.current[:] = salinity.forward[:]

    # UPDATE DENSITY
    vertical_density_profile = calculate_vertical_density_profile(temperature,salinity,vertical_grid)

    try:
        POM_only
    except NameError:
        POM_only = False
    else:
        POM_only = True
    if not POM_only:
        # pom_bfm_1d()
        pass

# WRITING OF RESTART
print(time,'\n')
print('-------------------------------------------------\n')
print(velocity.zonal_current,'\n', velocity.zonal_backward,'\n', velocity.meridional_current,'\n', velocity.meridional_backward,'\n')
print('-------------------------------------------------\n')
print(temperature.current,'\n', temperature.backward,'\n', salinity.current,'\n', salinity.backward,'\n')
print('-------------------------------------------------\n')
print(kinetic_energy.current,'\n', kinetic_energy.backward,'\n', kinetic_energy_times_length.current,'\n', kinetic_energy_times_length.backward,'\n')
print('-------------------------------------------------\n')
print(diffusion.tracers,'\n', diffusion.momentum,'\n', diffusion.kinetic_energy,'\n')
print('-------------------------------------------------\n')
print(vertical_grid.length_scale,'\n')
print('-------------------------------------------------\n')
print(bottom_stress.zonal,'\n', bottom_stress.meridional,'\n')
print('-------------------------------------------------\n')
print(vertical_density_profile)

# BFM RESTART
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
