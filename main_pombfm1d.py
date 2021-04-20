from cppdefs import *
import numpy as np
from pom.modules import zero, one, pi
from pom.modules import *
from inputs.params_POMBFM import *
from pom.calculations import *
from inputs.params_POMBFM import *  # Parameters read from params_POMBFM.nml, from main_pombfm1d.F90 lines 147-149
from pom.initialize_variables import *
from pom.create_profiles import *

# pyPOM1D DIRECTORY, USED FOR READING INPUTS (TO BE CHANGED BY USER)
current_path = '/Users/malikjordan/Desktop/pyPOM1D'



earth_angular_velocity = 7.29E-5                                                        # OMEGA

# GENERAL INITIALIZATION
vertical_coordinates = np.zeros(vertical_layers, dtype=float)                           # Z
vertical_coordinates_staggered = np.zeros(vertical_layers, dtype=float)                 # ZZ
vertical_spacing = np.zeros(vertical_layers,dtype=float)                                # DZ
vertical_spacing_staggered = np.zeros(vertical_layers,dtype=float)                      # DZZ
temperature_forward = np.zeros(vertical_layers,dtype=float)                             # TF
temperature = np.zeros(vertical_layers,dtype=float)                                     # T
temperature_backward = np.zeros(vertical_layers,dtype=float)                            # TB
salinity_forward = np.zeros(vertical_layers,dtype=float)                                # SF
salinity = np.zeros(vertical_layers,dtype=float)                                        # S
salinity_backward = np.zeros(vertical_layers,dtype=float)                               # SB
density_profile = np.zeros(vertical_layers,dtype=float)                                 # RHO
velocity_zonal_forward = np.zeros(vertical_layers,dtype=float)                          # UF
velocity_zonal = np.zeros(vertical_layers,dtype=float)                                  # U
velocity_zonal_backward = np.zeros(vertical_layers,dtype=float)                         # UB
velocity_meridional_forward = np.zeros(vertical_layers,dtype=float)                     # VF
velocity_meridional = np.zeros(vertical_layers,dtype=float)                             # V
velocity_meridional_backward = np.zeros(vertical_layers,dtype=float)                    # VB
kinetic_energy_forward = np.ones(vertical_layers,dtype=float) * 1.E-07                  # Q2F
kinetic_energy = np.ones(vertical_layers,dtype=float) * 1.E-07                          # Q2
kinetic_energy_backward = np.ones(vertical_layers,dtype=float) * 1.E-07                 # Q2B
kinetic_energy_times_length_forward = np.ones(vertical_layers,dtype=float) * 1.E-07     # Q2LF
kinetic_energy_times_length = np.ones(vertical_layers,dtype=float) * 1.E-07             # Q2L
kinetic_energy_times_length_backward = np.ones(vertical_layers,dtype=float) * 1.E-07    # Q2LB
length_scale = np.ones(vertical_layers,dtype=float)                                     # L
length_scale[0] = zero                                                                  # L(1)
length_scale[vertical_layers-1] = zero                                                  # L(KB)
diffusion_coefficient_momentum = np.zeros(vertical_layers,dtype=float)                  # KM
diffusion_coefficient_tracers = np.zeros(vertical_layers,dtype=float)                   # KH
diffusion_coefficient_kinetic_energy = np.zeros(vertical_layers,dtype=float)            # KQ
GM = np.zeros(vertical_layers,dtype=float)
GH = np.zeros(vertical_layers, dtype=float)
SM = np.zeros(vertical_layers,dtype=float)
SH = np.zeros(vertical_layers,dtype=float)
KN = np.zeros(vertical_layers,dtype=float)
SPROD = np.zeros(vertical_layers,dtype=float)
BPROD = np.zeros(vertical_layers,dtype=float)
PROD = np.zeros(vertical_layers,dtype=float)
VH = np.zeros(vertical_layers,dtype=float)
VHP = np.zeros(vertical_layers,dtype=float)
DTEF = np.zeros(vertical_layers,dtype=float)
D = np.zeros(vertical_layers,dtype=float)
DT = np.zeros(vertical_layers,dtype=float)
A = np.zeros(vertical_layers,dtype=float)
C = np.zeros(vertical_layers,dtype=float)
temperature_lateral_advection = np.zeros(vertical_layers, dtype=float)                  # WTADV
salinity_lateral_advection = np.zeros(vertical_layers, dtype=float)                     # WSADV
wind_stress_zonal = zero                                                                # WUSURF
wind_stress_meridional = zero                                                           # WVSURF
bottom_stress_zonal = zero                                                              # WUBOT
bottom_stress_meridional = zero                                                         # WVBOT
surface_heat_flux = zero                                                                # WTSURF
shortwave_radiation = zero                                                              # SWRAD
surface_salinity_flux = zero                                                            # WSSURF
surfave_temperature = zero                                                              # TSURF
surface_salinity = zero                                                                 # SSURF

# DEFINE VERTICAL COORDINATE SYSTEM
vertical_coordinates, vertical_coordinates_staggered, \
    vertical_spacing, vertical_spacing_staggered        = create_vertical_coordinate_system(vertical_layers, kl1, kl2)

vertical_spacing[vertical_layers-1] = 1.E-6  # Small value to avoid division by zero for vertical_spacing_reciprocal
vertical_spacing_reciprocal = np.zeros(vertical_layers,dtype=float)                     # DZR
for i in range(0,vertical_layers):
    vertical_spacing_reciprocal[i] = one / vertical_spacing[i]

vertical_spacing[vertical_layers] = zero

# CORIOLIS PARAMETER
coriolis_parameter = 2. * earth_angular_velocity * np.sin(alat * 2. * pi / 360.)        # COR

# TWICE THE TIMESTEP
twice_the_timestep = 2. * dti                                                           # DT2

# ITERATIONS NEEDED TO CARRY OUT AN"IDAYS" SIMULATION
iterations_needed = idays * seconds_per_day / dti                                       # iend

# READ  T&S INITIAL CONDITIONS (IHOTST=0) OR RESTART FILE (IHOTST=1)
if ihotst == 0:
    time0 = zero
    get_temperature_and_salinity_initial_coditions()
elif ihotst == 1:
    # get_rst()
    pass

try:
    import POM_only
    POM_only = True
except FileNotFoundError:
    POM_only = False
if not POM_only:
    # INITIALIZATION OF BFM
    # pom_ini_bfm_1d()
    pass

# BEGIN THE TIME MARCH
print('ICOUNT befor time march loop = ')
# print(icountf)


interpolated_temperature = np.zeros(vertical_layers,dtype=float)                        # TSTAR, from ModuleForcing.F90, line 78
interpolated_salinity = np.zeros(vertical_layers,dtype=float)                           # SSTAR, from ModuleForcing.F90, line 78
# BEGIN THE TIME MARCH

for i in range(0, iterations_needed):
    time = time0 + (dti * i * DAYI)

    # TURBULENCE CLOSURE
    for j in range(0, vertical_layers):
        kinetic_energy_forward[j] = kinetic_energy_backward[j]
        kinetic_energy_times_length_forward[j] = kinetic_energy_times_length_backward[j]

        # DEFINE ALL FORCINGS
        forcing_manager()

    # T&S COMPUTATION
    if idiagn == zero:
        # PROGNOSTIC MODE
        # T&S FULLY COMPUTED BY MODEL

        # COMPUTE LATERAL ADVECTION TERM FOR T&S
        # for j in range(0, vertical_layers):
        #     temperature_lateral_advection[j] = zero
        #     salinity_lateral_advection[j] = zero

        if trt == zero:
            for j in range(0, vertical_layers):
                if (-vertical_coordinates_staggered * bottom_depth) >= upperh:
                    temperature_lateral_advection[j] = (interpolated_temperature[j] - temperature[j]) / (trt * seconds_per_day)

        if srt == zero:
            for j in range(0, vertical_layers):
                if (-vertical_coordinates_staggered * bottom_depth) >= upperh:
                    salinity_lateral_advection[j] = (interpolated_salinity[j] - salinity[j]) / (srt * seconds_per_day)

        # COMPUTE SURFACE SALINITY FLUX
        surface_salinity_flux = -(surface_salinity - salinity[0]) * srt / seconds_per_day

        # COMPUTE TEMPREATURE
        for j in range(0, vertical_layers):
            temperature_forward[j] = temperature_backward[j] + (temperature_lateral_advection[j] * twice_the_timestep)
            calculate_vertical_temperature_and_salinity_profiles()

        # CALCULATE SALINITY
        for j in range(0, vertical_layers):
            salinity_forward[j] = salinity_backward[j] + (salinity_lateral_advection[j] * twice_the_timestep)
            calculate_vertical_temperature_and_salinity_profiles()

        # MIXING THE TIMESTEP (ASSELIN)
        for j in range(0, vertical_layers):
            temperature[j] = temperature[j] + 0.5 * smoth * (temperature_forward[j] + temperature_backward[j] - 2. * temperature_forward[j])
            salinity[j] = salinity[j] + 0.5 * smoth * (salinity_forward[j] + salinity_backward[j] - 2. * salinity[j])

    # COMPUTE VELOCITY
    for j in range(0, vertical_layers):
        velocity_zonal_forward[j] = velocity_zonal_backward[j] + twice_the_timestep * coriolis_parameter * velocity_meridional[j]
    calculate_vertical_zonal_velocity_profile(twice_the_timestep, bottom_depth, diffusion_coefficient_momentum, background_diffusion_momentum,
                                              vertical_spacing, vertical_spacing_staggered, wind_stress_zonal)

    for j in range(0, vertical_layers):
        velocity_meridional_forward[j] = velocity_meridional_backward[j] - twice_the_timestep * coriolis_parameter * velocity_zonal[j]
    calculate_vertical_meridional_velocity_profile(twice_the_timestep, bottom_depth, diffusion_coefficient_momentum, background_diffusion_momentum,
                                                   vertical_spacing, vertical_spacing_staggered, wind_stress_meridional)

    # MIX TIME STEL (ASSELIN FILTER)
    for j in range(0, vertical_layers):
        kinetic_energy[j] = kinetic_energy[j] + 0.5 * smoth * (kinetic_energy_forward[j] + kinetic_energy_backward[j] - 2. * kinetic_energy[j])
        kinetic_energy_times_length[j] = kinetic_energy_times_length[j] * smoth * (kinetic_energy_times_length_forward[j] + kinetic_energy_times_length_backward[j] - 2. * kinetic_energy_times_length[j])

        velocity_zonal[j] = velocity_zonal[j] + 0.5 * smoth * (velocity_zonal_forward[j] + velocity_zonal_backward[j] - 2. * velocity_zonal[j])
        velocity_meridional[j] = velocity_meridional[j] + 0.5 * smoth * (velocity_meridional_forward[j] + velocity_meridional_backward[j] - 2. * velocity_meridional[j])

    # RESTORE TIME SEQUENCE
    for j in range(0, vertical_layers):
        kinetic_energy_backward[j] = kinetic_energy[j]
        kinetic_energy[j] = kinetic_energy_forward[j]
        kinetic_energy_times_length_backward[j] = kinetic_energy_times_length[j]
        kinetic_energy_times_length[j] = kinetic_energy_times_length_forward[j]

    for j in range(0, vertical_layers):
        velocity_zonal_backward[j] = velocity_zonal[j]
        velocity_zonal[j] = velocity_zonal_forward[j]
        velocity_meridional_backward[j] = velocity_meridional[j]
        velocity_meridional[j] = velocity_meridional_forward[j]

    for j in range(0, vertical_layers):
        temperature_backward[j] = temperature[j]
        temperature[j] = temperature_forward[j]
        salinity_backward[j] = salinity[j]
        salinity[j] = salinity_forward[j]

    # UPDATE DENSITY
    create_vertical_diffusivity_profile(temperature, salinity, vertical_spacing_staggered, dti, vertical_layers)

try:
    import POM_only
    POM_only = True
except FileNotFoundError:
    POM_only = False
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
    import POM_only
    POM_only = True
except FileNotFoundError:
    POM_only = False
if not POM_only:
    # INITIALIZATION OF BFM
    # restart_BFM_inPOM()
    pass

print('Main done')
