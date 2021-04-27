import numpy as np
from pom.create_profiles import calculate_vertical_zonal_velocity_profile

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
wind_stress_zonal = 0.
diffusion_coefficient_momentum = 0.001 * np.ones(vertical_layers)
velocity_zonal_forward = np.zeros(vertical_layers)

# velocity_zonal_forward = calculate_vertical_zonal_velocity_profile(vertical_spacing,vertical_spacing_staggered,wind_stress_zonal,
#                                                                    diffusion_coefficient_momentum,velocity_zonal_forward)

# print("Test 1: Velocity Zonal Forward = \n",velocity_zonal_forward)
# print()

# Result: division by zero error
# ======================================================================================================================

# ======================================================================================================================
# Test 2 - Constant Profile
# constant vertical spacing and vertical spacing staggered --> [:] = -0.1
# constant diffusion coefficient --> [:]= 0.001
# velocity initialized to zero
# no wind stress

print("--------------------------------------------------------------")
print("Test 2")
print()

vertical_layers = 151
vertical_spacing = -0.1 * np.ones(vertical_layers)
vertical_spacing_staggered = -0.1 * np.ones(vertical_layers)
wind_stress_zonal = 0.
diffusion_coefficient_momentum = 0.001 * np.ones(vertical_layers)
velocity_zonal_forward = np.zeros(vertical_layers)

velocity_zonal_forward = calculate_vertical_zonal_velocity_profile(vertical_spacing,vertical_spacing_staggered,wind_stress_zonal,
                                                                   diffusion_coefficient_momentum,velocity_zonal_forward)

print("Test 2: Velocity Zonal Forward = \n",velocity_zonal_forward)
print()

# Result: constant velocity profile [:] = 0
# ======================================================================================================================

# ======================================================================================================================
# Test 3 - Wind Stress dependancy
# constant profiles for vertical spacing and vertical spacing staggered --> [:] = -0.1
# constant diffusion coefficient --> [:]= 0.001
# velocity initialized to zero
# wind stress = 0.1 vs 10

print("--------------------------------------------------------------")
print("Test 3")
print()

vertical_layers = 151
vertical_spacing = -0.1 * np.ones(vertical_layers)
vertical_spacing_staggered = -0.1 * np.ones(vertical_layers)
wind_stress_zonal = 0.1
diffusion_coefficient_momentum = 0.001 * np.ones(vertical_layers)
velocity_zonal_forward = np.zeros(vertical_layers)

velocity_zonal_forward = calculate_vertical_zonal_velocity_profile(vertical_spacing,vertical_spacing_staggered,wind_stress_zonal,
                                                                   diffusion_coefficient_momentum,velocity_zonal_forward)

print("Test 3: Velocity Zonal Forward (wind stress = 0.1) = \n",velocity_zonal_forward)
print()

wind_stress_zonal = 10
velocity_zonal_forward = np.zeros(vertical_layers)
velocity_zonal_forward = calculate_vertical_zonal_velocity_profile(vertical_spacing,vertical_spacing_staggered,wind_stress_zonal,
                                                                   diffusion_coefficient_momentum,velocity_zonal_forward)

print("Test 3: Velocity Zonal Forward (wind stress = 10) = \n",velocity_zonal_forward)
print()

# Result: velocity values are similar but differ by a factor of 100
#         larger wind stress values result in larger velocity values
# ======================================================================================================================

# ======================================================================================================================
# Test 3 - Velocity dependancy
# constant profiles for vertical spacing and vertical spacing staggered --> [:] = -0.1
# constant diffusion coefficient --> [:]= 0.001
# velocity profile constant (nonzero), linear, quadratic
# no wind stress

print("--------------------------------------------------------------")
print("Test 4")
print()

vertical_layers = 151
vertical_spacing = -0.1 * np.ones(vertical_layers)
vertical_spacing_staggered = -0.1 * np.ones(vertical_layers)
wind_stress_zonal = 0
diffusion_coefficient_momentum = 0.001 * np.ones(vertical_layers)
velocity_zonal_forward = 10 * np.ones(vertical_layers)

velocity_zonal_forward = calculate_vertical_zonal_velocity_profile(vertical_spacing,vertical_spacing_staggered,wind_stress_zonal,
                                                                   diffusion_coefficient_momentum,velocity_zonal_forward)

print("Test 4: Velocity Zonal Forward (constant initial velocity) = \n",velocity_zonal_forward)
print()

velocity_zonal_forward = np.linspace(10,40,vertical_layers)
velocity_zonal_forward = calculate_vertical_zonal_velocity_profile(vertical_spacing,vertical_spacing_staggered,wind_stress_zonal,
                                                                   diffusion_coefficient_momentum,velocity_zonal_forward)

print("Test 4: Velocity Zonal Forward (linear initial velocity) = \n",velocity_zonal_forward)
print()

velocity_zonal_forward = powspace(10,40,2,vertical_layers)
velocity_zonal_forward = calculate_vertical_zonal_velocity_profile(vertical_spacing,vertical_spacing_staggered,wind_stress_zonal,
                                                                   diffusion_coefficient_momentum,velocity_zonal_forward)

print("Test 4: Velocity Zonal Forward (quadratic initial velocity) = \n",velocity_zonal_forward)
print()

# Result: output velocity nearly returns the initial velocity profile
#         second to last entry on each return is zero, error in code
# ======================================================================================================================
