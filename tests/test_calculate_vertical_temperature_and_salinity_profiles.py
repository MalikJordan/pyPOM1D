import numpy as np
from pom.create_profiles import calculate_vertical_temperature_and_salinity_profiles

def powspace(start, stop, power, num):
    start = np.power(start, 1/float(power))
    stop = np.power(stop, 1/float(power))
    return np.power(np.linspace(start, stop, num=num), power)


# ======================================================================================================================
# Test 1 - Division by zero
# vertical spacing and vertical spacing staggered = 0 (should cause a division by zero error in first loop)

print("--------------------------------------------------------------")
print("Test 1")
print()

vertical_layers = 151
vertical_spacing = np.zeros(vertical_layers)
vertical_spacing_staggered = np.zeros(vertical_layers)
vertical_coordinates = np.zeros(vertical_layers)
vertical_coordinates_staggered = np.zeros(vertical_layers)
diffusion_coefficient_tracers = 0.001 * np.ones(vertical_layers)

property_forward_profile = np.zeros(vertical_layers)  # FF
property_surface_flux = 0.  # WFSURF
property_bottom_flux = 0.  # WFBOT
shortwave_radiation = 0.  # SWRAD
property_surface_value = 0.  # FSURF
nbc = 1

# property_forward_profile = calculate_vertical_temperature_and_salinity_profiles(vertical_spacing, vertical_spacing_staggered, vertical_coordinates, vertical_coordinates_staggered,
#                                                                                 diffusion_coefficient_tracers, property_forward_profile, property_surface_flux, property_bottom_flux,
#                                                                                 shortwave_radiation, property_surface_value, nbc)

# Result: division by zero error
# ======================================================================================================================

# ======================================================================================================================
# Test 2 - Constant Profile
# constant vertical spacing and vertical spacing staggered --> [:] = -0.1
# constant diffusion coefficient --> [:]= 0.001
# profile initialized to zero
# no flux or shortwave radiation

print("--------------------------------------------------------------")
print("Test 2")
print()

vertical_layers = 151
vertical_spacing = -0.1 * np.ones(vertical_layers)
vertical_spacing_staggered = -0.1 * np.ones(vertical_layers)
vertical_coordinates = np.linspace(0,-15,vertical_layers)
vertical_coordinates_staggered = np.linspace(-0.05,-15.05,vertical_layers)
diffusion_coefficient_tracers = 0.001 * np.ones(vertical_layers)

property_forward_profile = np.zeros(vertical_layers)
property_surface_flux = 0
property_bottom_flux = 0
shortwave_radiation = 0
property_surface_value = 0
nbc = 1

property_forward_profile = calculate_vertical_temperature_and_salinity_profiles(vertical_spacing, vertical_spacing_staggered, vertical_coordinates, vertical_coordinates_staggered,
                                                                                diffusion_coefficient_tracers, property_forward_profile, property_surface_flux, property_bottom_flux,
                                                                                shortwave_radiation, property_surface_value, nbc)

print("Test 2: Property Foward Profile = \n",property_forward_profile)
print()

# Result: constant profile --> [:] = 0
# ======================================================================================================================

# ======================================================================================================================
# Test 3 - Flux depandancy
# constant vertical spacing and vertical spacing staggered --> [:] = -0.1
# constant diffusion coefficient --> [:]= 0.001
# profile initialized to zero
# flux and radiation = 1 vs 10

print("--------------------------------------------------------------")
print("Test 3")
print()

vertical_layers = 151
vertical_spacing = -0.1 * np.ones(vertical_layers)
vertical_spacing_staggered = -0.1 * np.ones(vertical_layers)
vertical_coordinates = np.linspace(0,-15,vertical_layers)
vertical_coordinates_staggered = np.linspace(-0.05,-15.05,vertical_layers)
diffusion_coefficient_tracers = 0.001 * np.ones(vertical_layers)

property_forward_profile = np.zeros(vertical_layers)
property_surface_flux = 1
property_bottom_flux = 1
shortwave_radiation = 1
property_surface_value = 1
nbc = 1

property_forward_profile = calculate_vertical_temperature_and_salinity_profiles(vertical_spacing, vertical_spacing_staggered, vertical_coordinates, vertical_coordinates_staggered,
                                                                                diffusion_coefficient_tracers, property_forward_profile, property_surface_flux, property_bottom_flux,
                                                                                shortwave_radiation, property_surface_value, nbc)

print("Test 3: Property Foward Profile (flux = 1) = \n",property_forward_profile)
print()

property_forward_profile = np.zeros(vertical_layers)
property_surface_flux = 10
property_bottom_flux = 10
shortwave_radiation = 10
property_surface_value = 10
nbc = 1

property_forward_profile = calculate_vertical_temperature_and_salinity_profiles(vertical_spacing, vertical_spacing_staggered, vertical_coordinates, vertical_coordinates_staggered,
                                                                                diffusion_coefficient_tracers, property_forward_profile, property_surface_flux, property_bottom_flux,
                                                                                shortwave_radiation, property_surface_value, nbc)

print("Test 3: Property Foward Profile (flux = 10) = \n",property_forward_profile)
print()

# Result: similar values, differ by a factor of 10
#         larger flux results in larger profile values
# ======================================================================================================================

# ======================================================================================================================
# Test 4 - Property profile depandancy
# constant vertical spacing and vertical spacing staggered --> [:] = -0.1
# constant diffusion coefficient --> [:]= 0.001
# profile constant nonzero value, linear, quatratic
# flux and radiation = 1

print("--------------------------------------------------------------")
print("Test 4")
print()

vertical_layers = 151
vertical_spacing = -0.1 * np.ones(vertical_layers)
vertical_spacing_staggered = -0.1 * np.ones(vertical_layers)
vertical_coordinates = np.linspace(0,-15,vertical_layers)
vertical_coordinates_staggered = np.linspace(-0.05,-15.05,vertical_layers)
diffusion_coefficient_tracers = 0.001 * np.ones(vertical_layers)

property_forward_profile = 20 * np.ones(vertical_layers)
property_surface_flux = 1
property_bottom_flux = 1
shortwave_radiation = 1
property_surface_value = 1
nbc = 1

property_forward_profile = calculate_vertical_temperature_and_salinity_profiles(vertical_spacing, vertical_spacing_staggered, vertical_coordinates, vertical_coordinates_staggered,
                                                                                diffusion_coefficient_tracers, property_forward_profile, property_surface_flux, property_bottom_flux,
                                                                                shortwave_radiation, property_surface_value, nbc)

print("Test 4: Property Foward Profile (constant nonzero initial profile) = \n",property_forward_profile)
print()

property_forward_profile = np.linspace(50,20,vertical_layers)
property_surface_flux = 1
property_bottom_flux = 1
shortwave_radiation = 1
property_surface_value = 1
nbc = 1

property_forward_profile = calculate_vertical_temperature_and_salinity_profiles(vertical_spacing, vertical_spacing_staggered, vertical_coordinates, vertical_coordinates_staggered,
                                                                                diffusion_coefficient_tracers, property_forward_profile, property_surface_flux, property_bottom_flux,
                                                                                shortwave_radiation, property_surface_value, nbc)

print("Test 4: Property Foward Profile (linear initial profile) = \n",property_forward_profile)
print()

property_forward_profile = powspace(50,20,2,vertical_layers)
property_surface_flux = 1
property_bottom_flux = 1
shortwave_radiation = 1
property_surface_value = 1
nbc = 1

property_forward_profile = calculate_vertical_temperature_and_salinity_profiles(vertical_spacing, vertical_spacing_staggered, vertical_coordinates, vertical_coordinates_staggered,
                                                                                diffusion_coefficient_tracers, property_forward_profile, property_surface_flux, property_bottom_flux,
                                                                                shortwave_radiation, property_surface_value, nbc)

print("Test 4: Property Foward Profile (quadratic initial profile) = \n",property_forward_profile)
print()

# Result: function nearly returned initial profile
#         first few indices are larger than the initial profile, likely due to flux
#         second to last indice is zero for each case, error in code
# ======================================================================================================================

# ======================================================================================================================
# Test 5 - Diffusion depandancy
# constant vertical spacing and vertical spacing staggered --> [:] = -0.1
# constant diffusion coefficient ([:]= 0.001 vs 1), linear, quadratic
# profile initialized to zero
# flux and radiation = 1

print("--------------------------------------------------------------")
print("Test 5")
print()

vertical_layers = 151
vertical_spacing = -0.1 * np.ones(vertical_layers)
vertical_spacing_staggered = -0.1 * np.ones(vertical_layers)
vertical_coordinates = np.linspace(0,-15,vertical_layers)
vertical_coordinates_staggered = np.linspace(-0.05,-15.05,vertical_layers)
diffusion_coefficient_tracers = 0.001 * np.ones(vertical_layers)

property_forward_profile = np.zeros(vertical_layers)
property_surface_flux = 1
property_bottom_flux = 1
shortwave_radiation = 1
property_surface_value = 1
nbc = 1

property_forward_profile = calculate_vertical_temperature_and_salinity_profiles(vertical_spacing, vertical_spacing_staggered, vertical_coordinates, vertical_coordinates_staggered,
                                                                                diffusion_coefficient_tracers, property_forward_profile, property_surface_flux, property_bottom_flux,
                                                                                shortwave_radiation, property_surface_value, nbc)

print("Test 5: Property Foward Profile (constant diffusion profile= 0.001) = \n",property_forward_profile)
print()

diffusion_coefficient_tracers = np.ones(vertical_layers)
property_forward_profile = np.zeros(vertical_layers)

property_forward_profile = calculate_vertical_temperature_and_salinity_profiles(vertical_spacing, vertical_spacing_staggered, vertical_coordinates, vertical_coordinates_staggered,
                                                                                diffusion_coefficient_tracers, property_forward_profile, property_surface_flux, property_bottom_flux,
                                                                                shortwave_radiation, property_surface_value, nbc)

print("Test 5: Property Foward Profile (constant diffusion profile = 1) = \n",property_forward_profile)
print()

diffusion_coefficient_tracers = np.linspace(0,15,vertical_layers)
property_forward_profile = np.zeros(vertical_layers)

property_forward_profile = calculate_vertical_temperature_and_salinity_profiles(vertical_spacing, vertical_spacing_staggered, vertical_coordinates, vertical_coordinates_staggered,
                                                                                diffusion_coefficient_tracers, property_forward_profile, property_surface_flux, property_bottom_flux,
                                                                                shortwave_radiation, property_surface_value, nbc)

print("Test 5: Property Foward Profile (linear diffusion profile) = \n",property_forward_profile)
print()

diffusion_coefficient_tracers = powspace(0,15,2,vertical_layers)
property_forward_profile = np.zeros(vertical_layers)

property_forward_profile = calculate_vertical_temperature_and_salinity_profiles(vertical_spacing, vertical_spacing_staggered, vertical_coordinates, vertical_coordinates_staggered,
                                                                                diffusion_coefficient_tracers, property_forward_profile, property_surface_flux, property_bottom_flux,
                                                                                shortwave_radiation, property_surface_value, nbc)

print("Test 5: Property Foward Profile (quadratic diffusion profile) = \n",property_forward_profile)
print()

# Result: for constant profile, values are larger for larger diffusion coefficient
#         for changing profile, values decrease at a faster rate for quadratic vs linear
# ======================================================================================================================

# ======================================================================================================================
# Test 6 - nbc depandancy
# constant vertical spacing and vertical spacing staggered --> [:] = -0.1
# constant diffusion coefficient --> [:]= 0.001
# profile initialized to zero
# flux and radiation = 1

print("--------------------------------------------------------------")
print("Test 6")
print()

vertical_layers = 151
vertical_spacing = -0.1 * np.ones(vertical_layers)
vertical_spacing_staggered = -0.1 * np.ones(vertical_layers)
vertical_coordinates = np.linspace(0,-15,vertical_layers)
vertical_coordinates_staggered = np.linspace(-0.05,-15.05,vertical_layers)
diffusion_coefficient_tracers = 0.001 * np.ones(vertical_layers)

property_forward_profile = np.zeros(vertical_layers)
property_surface_flux = 1
property_bottom_flux = 1
shortwave_radiation = 1
property_surface_value = 1
nbc = 1

property_forward_profile = calculate_vertical_temperature_and_salinity_profiles(vertical_spacing, vertical_spacing_staggered, vertical_coordinates, vertical_coordinates_staggered,
                                                                                diffusion_coefficient_tracers, property_forward_profile, property_surface_flux, property_bottom_flux,
                                                                                shortwave_radiation, property_surface_value, nbc)

print("Test 6: Property Foward Profile (nbc = 1) = \n",property_forward_profile)
print()

property_forward_profile = np.zeros(vertical_layers)
nbc = 2

property_forward_profile = calculate_vertical_temperature_and_salinity_profiles(vertical_spacing, vertical_spacing_staggered, vertical_coordinates, vertical_coordinates_staggered,
                                                                                diffusion_coefficient_tracers, property_forward_profile, property_surface_flux, property_bottom_flux,
                                                                                shortwave_radiation, property_surface_value, nbc)

print("Test 6: Property Foward Profile (nbc = 2) = \n",property_forward_profile)
print()

property_forward_profile = np.zeros(vertical_layers)
nbc = 3

property_forward_profile = calculate_vertical_temperature_and_salinity_profiles(vertical_spacing, vertical_spacing_staggered, vertical_coordinates, vertical_coordinates_staggered,
                                                                                diffusion_coefficient_tracers, property_forward_profile, property_surface_flux, property_bottom_flux,
                                                                                shortwave_radiation, property_surface_value, nbc)

print("Test 6: Property Foward Profile (nbc = 3) = \n",property_forward_profile)
print()

property_forward_profile = np.zeros(vertical_layers)
nbc = 4

property_forward_profile = calculate_vertical_temperature_and_salinity_profiles(vertical_spacing, vertical_spacing_staggered, vertical_coordinates, vertical_coordinates_staggered,
                                                                                diffusion_coefficient_tracers, property_forward_profile, property_surface_flux, property_bottom_flux,
                                                                                shortwave_radiation, property_surface_value, nbc)

print("Test 6: Property Foward Profile (nbc = 4) = \n",property_forward_profile)
print()

# Result: different profiles for different nbc values
#         nbc = 2 and nbc = 4 return nearly identical profiles
# ======================================================================================================================




