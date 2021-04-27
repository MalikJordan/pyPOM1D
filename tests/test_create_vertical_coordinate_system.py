import numpy as np
from pom.calculations import create_vertical_coordinate_system

# ======================================================================================================================
# Test 1 - Division by zero
# bottom layers with log distriu=bution - surface layers with log distribution + 4 = 0 --> kl2 - kl1 + 4 = 0
# Expected Result: division by zero in second loop will break code

vertical_layers = 61
surface_layers_with_log_distribution = 6
bottom_layers_with_log_distribution = 2

# vertical_coordinates, vertical_coordinates_staggered, vertical_spacing, vertical_spacing_staggered \
#     = create_vertical_coordinate_system(vertical_layers,surface_layers_with_log_distribution,bottom_layers_with_log_distribution)

# Result: division by zero
# ======================================================================================================================

# ======================================================================================================================
# Test 2 - kl1 dependancy
# kl1 = 1 vs kl1 = 5
# hold kl2 = 2

print("--------------------------------------------------------------")
print("Test 2")
print()

vertical_layers = 61
surface_layers_with_log_distribution = 1
bottom_layers_with_log_distribution = 2

vertical_coordinates, vertical_coordinates_staggered, vertical_spacing, vertical_spacing_staggered \
    = create_vertical_coordinate_system(vertical_layers,surface_layers_with_log_distribution,bottom_layers_with_log_distribution)

print("Test 2: Vertical coorniates (kl1 = 1, kl2 = 2) = \n", vertical_coordinates)
print()
print("Test 2: Vertical coordinates staggered (kl1 = 1, kl2 = 2) = \n", vertical_coordinates_staggered)
print()

vertical_layers = 61
surface_layers_with_log_distribution = 5
bottom_layers_with_log_distribution = 2

vertical_coordinates, vertical_coordinates_staggered, vertical_spacing, vertical_spacing_staggered \
    = create_vertical_coordinate_system(vertical_layers,surface_layers_with_log_distribution,bottom_layers_with_log_distribution)

print("Test 2: Vertical coorniates (kl1 = 5, kl2 = 2) = \n", vertical_coordinates)
print()
print("Test 2: Vertical coordinates staggered (kl1 = 5, kl2 = 2) = \n", vertical_coordinates_staggered)
print()

# Result: kl1 = 1 --> first vertical coordinate is -0.2 when it should  be 0
#                     linear spacing for between all points for vertical coordinates and staggered vertical coordinates
#         kl1 = 5 --> first vertical coordinate is 0 as expected
#                     logarithmic spacing between first 4 coordinates
#                     first staggerred vertical coordinate is > 0 (above ocean surface) because kl2 - kl1 + 4 > 0
# ======================================================================================================================

# ======================================================================================================================
# Test 3 - kl2 dependancy
# kl2 = 1, kl2 = 5, kl2 = 15
# hold kl1 = 10

print("--------------------------------------------------------------")
print("Test 3")
print()

vertical_layers = 61
surface_layers_with_log_distribution = 10
bottom_layers_with_log_distribution = 1

vertical_coordinates, vertical_coordinates_staggered, vertical_spacing, vertical_spacing_staggered \
    = create_vertical_coordinate_system(vertical_layers,surface_layers_with_log_distribution,bottom_layers_with_log_distribution)

print("Test 3: Vertical coorniates (kl1 = 10, kl2 = 1) = \n", vertical_coordinates)
print()
print("Test 3: Vertical coordinates staggered (kl1 = 10, kl2 = 1) = \n", vertical_coordinates_staggered)
print()

vertical_layers = 61
surface_layers_with_log_distribution = 10
bottom_layers_with_log_distribution = 5

vertical_coordinates, vertical_coordinates_staggered, vertical_spacing, vertical_spacing_staggered \
    = create_vertical_coordinate_system(vertical_layers,surface_layers_with_log_distribution,bottom_layers_with_log_distribution)

print("Test 3: Vertical coorniates (kl1 = 10, kl2 = 5) = \n", vertical_coordinates)
print()
print("Test 3: Vertical coordinates staggered (kl1 = 10, kl2 = 5) = \n", vertical_coordinates_staggered)
print()

vertical_layers = 61
surface_layers_with_log_distribution = 10
bottom_layers_with_log_distribution = 15

vertical_coordinates, vertical_coordinates_staggered, vertical_spacing, vertical_spacing_staggered \
    = create_vertical_coordinate_system(vertical_layers,surface_layers_with_log_distribution,bottom_layers_with_log_distribution)

print("Test 3: Vertical coorniates (kl1 = 10, kl2 = 15) = \n", vertical_coordinates)
print()
print("Test 3: Vertical coordinates staggered (kl1 = 10, kl2 = 15) = \n", vertical_coordinates_staggered)
print()

# Result: profile increases for kl2 = 1 and kl2 = 5 (expected to increase whenever kl2 - kl1 + 4 < 0)
#         logarithmic distribution between first 9 coordinates, as expected
#         first staggerred vertical coordinate for kl2 = 15 is > 0 (above ocean surface) because kl2 - kl1 + 4 > 0
# ======================================================================================================================













