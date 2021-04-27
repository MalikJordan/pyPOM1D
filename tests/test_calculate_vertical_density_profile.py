import numpy as np
from pom.calculations import calculate_vertical_density_profile

def powspace(start, stop, power, num):
    start = np.power(start, 1/float(power))
    stop = np.power(stop, 1/float(power))
    return np.power(np.linspace(start, stop, num=num), power)


# ======================================================================================================================
# Test 1 - Constant profile
# Linear spacing with 61 points --> dz = 0.1
# One hour --> timestep = 60 s
# Constant temperature --> T = 20 C
# Constant salinity --> S = 35 psu

print("--------------------------------------------------------------")
print("Test 1")
print()

vertical_layers = 61
timestep = 60
vertical_spacing_staggered = np.linspace(0,6,vertical_layers)
temperature = 20 * np.ones(vertical_layers)
salinity = 35 * np.ones(vertical_layers)

vertical_density_profile = calculate_vertical_density_profile(temperature,salinity,vertical_spacing_staggered,timestep,vertical_layers)
density_profile_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    density_profile_spacing[i] = vertical_density_profile[i+1] - vertical_density_profile[i]
print("Test 1: vertical density profile = \n", vertical_density_profile)
print()
print("Test 1: density profile spacing = \n", density_profile_spacing)
print()

# Result: constant density profile
# ======================================================================================================================


# ======================================================================================================================
# Test 2 - Temperature dependancy
# Linear spacing with 61 points --> dz = 0.1
# One hour --> timestep = 60 s
# Linear increase in temperature --> 10 C < T < 20 C
# Constant salinity --> S = 35 psu

print("--------------------------------------------------------------")
print("Test 2")
print()

vertical_layers = 61
timestep = 60
vertical_spacing_staggered = np.linspace(0,6,vertical_layers)
temperature = np.linspace(10,20,vertical_layers)
temperature_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    temperature_spacing[i] = temperature[i+1] - temperature[i]
print("Test 2: temperature = \n", temperature)
print()
print("Test 2: temperature spacing = \n", temperature_spacing)
print()
salinity = 35 * np.ones(vertical_layers)

vertical_density_profile = calculate_vertical_density_profile(temperature,salinity,vertical_spacing_staggered,timestep,vertical_layers)
density_profile_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    density_profile_spacing[i] = vertical_density_profile[i+1] - vertical_density_profile[i]
print("Test 2: vertical density profile = \n", vertical_density_profile)
print()
print("Test 2: density profile spacing = \n", density_profile_spacing)
print()

# Result: linearly decreasing density with linearly increasing temperature
# Further testing: - lower values of temperature (potentially test negative to positive)
#                  - decreasing temperature profile (ocean temperature decreases with depth)
# ======================================================================================================================


# ======================================================================================================================
# Test 3 - Temperature dependancy
# Linear spacing with 61 points --> dz = 0.1
# One hour --> timestep = 60 s
# Logarithmic increase in temperature --> 10 C < T < 20 C
# Constant salinity --> S = 35 psu

print("--------------------------------------------------------------")
print("Test 3")
print()

vertical_layers = 61
timestep = 60
vertical_spacing_staggered = np.linspace(0,6,vertical_layers)
temperature = np.logspace(1,np.log(20)/np.log(10),vertical_layers)
temperature_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    temperature_spacing[i] = temperature[i+1] - temperature[i]
print("Test 3: temperature = \n", temperature)
print()
print("Test 3: temperature spacing = \n", temperature_spacing)
print()
salinity = 35 * np.ones(vertical_layers)

vertical_density_profile = calculate_vertical_density_profile(temperature,salinity,vertical_spacing_staggered,timestep,vertical_layers)
density_profile_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    density_profile_spacing[i] = vertical_density_profile[i+1] - vertical_density_profile[i]
print("Test 3: vertical density profile = \n", vertical_density_profile)
print()
print("Test 3: density profile spacing = \n", density_profile_spacing)
print()

# Result: decreasing density with increasing temperature
# ======================================================================================================================

# ======================================================================================================================
# Test 4 - Temperature dependancy
# Linear spacing with 61 points --> dz = 0.1
# One hour --> timestep = 60 s
# 10th power increase in temperature --> 10 C < T < 20 C
# Constant salinity --> S = 35 psu

print("--------------------------------------------------------------")
print("Test 4")
print()

vertical_layers = 61
timestep = 60
vertical_spacing_staggered = np.linspace(0,6,vertical_layers)
temperature = powspace(10,20,10,vertical_layers)
temperature_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    temperature_spacing[i] = temperature[i+1] - temperature[i]
print("Test 4: temperature = \n", temperature)
print()
print("Test 4: temperature spacing = \n", temperature_spacing)
print()
salinity = 35 * np.ones(vertical_layers)

vertical_density_profile = calculate_vertical_density_profile(temperature,salinity,vertical_spacing_staggered,timestep,vertical_layers)
density_profile_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    density_profile_spacing[i] = vertical_density_profile[i+1] - vertical_density_profile[i]
print("Test 4: vertical density profile = \n", vertical_density_profile)
print()
print("Test 4: density profile spacing = \n", density_profile_spacing)
print()

# Result: decreasing density with increasing temperature
# ======================================================================================================================

# ======================================================================================================================
# Test 5 - Temperature dependancy
# Linear spacing with 61 points --> dz = 0.1
# One hour --> timestep = 60 s
# Linear Decrease in temperature --> 10 C < T < 20 C
# Constant salinity --> S = 35 psu

print("--------------------------------------------------------------")
print("Test 5")
print()

vertical_layers = 61
timestep = 60
vertical_spacing_staggered = np.linspace(0,6,vertical_layers)
temperature = np.linspace(20,10,vertical_layers)
temperature_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    temperature_spacing[i] = temperature[i+1] - temperature[i]
print("Test 5: temperature = \n", temperature)
print()
print("Test 5: temperature spacing = \n", temperature_spacing)
print()
salinity = 35 * np.ones(vertical_layers)

vertical_density_profile = calculate_vertical_density_profile(temperature,salinity,vertical_spacing_staggered,timestep,vertical_layers)
density_profile_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    density_profile_spacing[i] = vertical_density_profile[i+1] - vertical_density_profile[i]
print("Test 5: vertical density profile = \n", vertical_density_profile)
print()
print("Test 5: density profile spacing = \n", density_profile_spacing)
print()

# Result: increasing density with decreasing temperature
# ======================================================================================================================

# ======================================================================================================================
# Test 6 - Salinity dependancy
# Linear spacing with 61 points --> dz = 0.1
# One hour --> timestep = 60 s
# Constant temperature --> T = 20 C
# Linearly increasing salinity --> 33 psu < S < 37 psu

print("--------------------------------------------------------------")
print("Test 6")
print()

vertical_layers = 61
timestep = 60
vertical_spacing_staggered = np.linspace(0,6,vertical_layers)
temperature = 20 * np.ones(vertical_layers)
salinity = np.linspace(33,37,vertical_layers)
salinity_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    salinity_spacing[i] = salinity[i+1] - salinity[i]
print("Test 6: salinity = \n", salinity)
print()
print("Test 6: salinity spacing = \n", salinity_spacing)
print()

vertical_density_profile = calculate_vertical_density_profile(temperature,salinity,vertical_spacing_staggered,timestep,vertical_layers)
density_profile_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    density_profile_spacing[i] = vertical_density_profile[i+1] - vertical_density_profile[i]
print("Test 6: vertical density profile = \n", vertical_density_profile)
print()
print("Test 6: density profile spacing = \n", density_profile_spacing)
print()

# Result: increasing density with increasing salinity
# ======================================================================================================================


# ======================================================================================================================
# Test 7 - Salinity dependancy
# Linear spacing with 61 points --> dz = 0.1
# One hour --> timestep = 60 s
# Constant temperature --> T = 20 C
# 10th power increase salinity --> 33 psu < S < 37 psu

print("--------------------------------------------------------------")
print("Test 7")
print()

vertical_layers = 61
timestep = 60
vertical_spacing_staggered = np.linspace(0,6,vertical_layers)
temperature = 20 * np.ones(vertical_layers)
salinity = powspace(33,37,10,vertical_layers)
salinity_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    salinity_spacing[i] = salinity[i+1] - salinity[i]
print("Test 7: salinity = \n", salinity)
print()
print("Test 7: salinity spacing = \n", salinity_spacing)
print()

vertical_density_profile = calculate_vertical_density_profile(temperature,salinity,vertical_spacing_staggered,timestep,vertical_layers)
density_profile_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    density_profile_spacing[i] = vertical_density_profile[i+1] - vertical_density_profile[i]
print("Test 7: vertical density profile = \n", vertical_density_profile)
print()
print("Test 7: density profile spacing = \n", density_profile_spacing)
print()

# Result: increasing density with increasing salinity
# ======================================================================================================================

