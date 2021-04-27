import numpy as np
from pom.calculations import calculate_mixed_layer_depth

# ======================================================================================================================
# Test 1 - Constant profile
# Constant vertical coordinate profile --> [:] = 0 (surface)
# Constant Temperature --> T = 20 C

print("--------------------------------------------------------------")
print("Test 1")
print()

vertical_layers = 61
vertical_coordinates_staggered = np.zeros(vertical_layers)
temperature = 20 * np.ones(vertical_layers)

mixed_layer_depth = calculate_mixed_layer_depth(vertical_coordinates_staggered,temperature,vertical_layers)
print("Test 1: Mixed Layer Depth = \n", mixed_layer_depth)
print()

# Result: constant profile, mixed layer depth = zero (ocean surface)
# ======================================================================================================================

# ======================================================================================================================
# Test 2 - Temperature dependancy
# Linearly spaced vertical coordinates with 61 points --> dz = 0.1
# Decreasing vertical coordinates (0 at surface)
# Constant Temperature --> T = 20 C

print("--------------------------------------------------------------")
print("Test 2")
print()

vertical_layers = 61
vertical_coordinates_staggered = np.linspace(0,-6,vertical_layers)
temperature = 20 * np.ones(vertical_layers)

mixed_layer_depth = calculate_mixed_layer_depth(vertical_coordinates_staggered,temperature,vertical_layers)
print("Test 2: Mixed Layer Depth = \n", mixed_layer_depth)
print()

# Result: very negative mixed layer depth, decreasing by dz
# ======================================================================================================================

# ======================================================================================================================
# Test 3 - Temperature dependancy
# Linearly spaced vertical coordinates with 61 points --> dz = 0.1
# Decreasing vertical coordinates (0 at surface)
# Linearly decreasing Temperature with small dTemp --> 20 C > T > 19.8 C

print("--------------------------------------------------------------")
print("Test 3")
print()

vertical_layers = 61
vertical_coordinates_staggered = np.linspace(0,-6,vertical_layers)
temperature = np.linspace(20,19.8,vertical_layers)
temperature_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    temperature_spacing[i] = temperature[i+1] - temperature[i]
print("Test 3: Temperature = \n", temperature)
print()
print("Test 3: Temperature spacing = \n", temperature_spacing)
print()

mixed_layer_depth = calculate_mixed_layer_depth(vertical_coordinates_staggered,temperature,vertical_layers)
mixed_layer_depth_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    mixed_layer_depth_spacing[i] = mixed_layer_depth[i+1] - mixed_layer_depth[i]
print("Test 3: Mixed Layer Depth = \n", mixed_layer_depth)
print()
print("Test 3: Mixed Layer Depth Spacing = \n", mixed_layer_depth_spacing)
print()

# Result: decreasing mixed layer depth with decreasing temperature
# ======================================================================================================================

# ======================================================================================================================
# Test 4 - Temperature dependancy
# Linearly spaced vertical coordinates with 61 points --> dz = 0.1
# Decreasing vertical coordinates (0 at surface)
# Linearly decreasing Temperature with large dTemp --> 20 C > T > 10 C

print("--------------------------------------------------------------")
print("Test 4")
print()

vertical_layers = 61
vertical_coordinates_staggered = np.linspace(0,-6,vertical_layers)
temperature = np.linspace(20,10,vertical_layers)
temperature_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    temperature_spacing[i] = temperature[i+1] - temperature[i]
print("Test 4: Temperature = \n", temperature)
print()
print("Test 4: Temperature spacing = \n", temperature_spacing)
print()

mixed_layer_depth = calculate_mixed_layer_depth(vertical_coordinates_staggered,temperature,vertical_layers)
mixed_layer_depth_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    mixed_layer_depth_spacing[i] = mixed_layer_depth[i+1] - mixed_layer_depth[i]
print("Test 4: Mixed Layer Depth = \n", mixed_layer_depth)
print()
print("Test 4: Mixed Layer Depth Spacing = \n", mixed_layer_depth_spacing)
print()

# Result: loop broke after, confirmed break condition
# ======================================================================================================================

# ======================================================================================================================
# Test 5 - Vertical Coordinate dependancy
# Linearly spaced vertical coordinates with 61 points --> dz = 0.1
# Vertical coordinates starting above surface --> [:] > 0
# Linearly decreasing Temperature with small dTemp --> 20 C > T > 19.8 C

print("--------------------------------------------------------------")
print("Test 5")
print()

vertical_layers = 61
vertical_coordinates_staggered = np.linspace(7,1,vertical_layers)
temperature = np.linspace(20,19.8,vertical_layers)
temperature_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    temperature_spacing[i] = temperature[i+1] - temperature[i]
print("Test 5: Temperature = \n", temperature)
print()
print("Test 5: Temperature spacing = \n", temperature_spacing)
print()

mixed_layer_depth = calculate_mixed_layer_depth(vertical_coordinates_staggered,temperature,vertical_layers)
mixed_layer_depth_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    mixed_layer_depth_spacing[i] = mixed_layer_depth[i+1] - mixed_layer_depth[i]
print("Test 5: Mixed Layer Depth = \n", mixed_layer_depth)
print()
print("Test 5: Mixed Layer Depth Spacing = \n", mixed_layer_depth_spacing)
print()

# Result: positive mixed layer depth (above surface), decreasing
# ======================================================================================================================

# ======================================================================================================================
# Test 6 - Temperature dependancy
# Linearly spaced vertical coordinates with 61 points --> dz = 0.1
# Increasing values for vertical coordinates
# Linearly decreasing Temperature with small dTemp --> 20 C > T > 19.8 C

print("--------------------------------------------------------------")
print("Test 6")
print()

vertical_layers = 61
vertical_coordinates_staggered = np.linspace(-3,3,vertical_layers)
temperature = np.linspace(20,19.8,vertical_layers)
temperature_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    temperature_spacing[i] = temperature[i+1] - temperature[i]
print("Test 6: Temperature = \n", temperature)
print()
print("Test 6: Temperature spacing = \n", temperature_spacing)
print()

mixed_layer_depth = calculate_mixed_layer_depth(vertical_coordinates_staggered,temperature,vertical_layers)
mixed_layer_depth_spacing = np.zeros(vertical_layers-1)
for i in range(0,vertical_layers-1):
    mixed_layer_depth_spacing[i] = mixed_layer_depth[i+1] - mixed_layer_depth[i]
print("Test 6: Mixed Layer Depth = \n", mixed_layer_depth)
print()
print("Test 6: Mixed Layer Depth Spacing = \n", mixed_layer_depth_spacing)
print()

# Result: all positive values for mixed layer depth, increasing
# ======================================================================================================================
