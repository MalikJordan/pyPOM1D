# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   TESTS FOR POM CALCULATIONS
#       - calcdepth
#       - dens
#       - mldpth
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from pom.calculations import CALCDEPTH, DENS, MLDPTH
import numpy as np
import matplotlib.pyplot as plt

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# CALCDEPTH
#
# Default values provided my documentation: KL1 = .3*KB, KL2 = KB-2.
# Tests will examine the calculation of the vertical distributions Z, ZZ, DZ, and DZZ
# by changing the input value of KB. KL1 and KL2 will be calculated as above.
#
# KB  = # vertical layers (currently set to 31 in namelist)
# KL1 = # surface vertical leyers with log dist
# KL2 = # bottom vertical leyers with log dist

KB = [14, 21, 31, 61]
plt.figure(1)

for i in range(0, len(KB)):
    KL1 = np.floor(0.3 * KB[i])
    KL2 = KB[i] - 2
    x = np.ones(KB[i])

    Z = np.empty(KB[i])
    ZZ = np.empty(KB[i])
    DZ = np.empty(KB[i])
    DZZ = np.empty(KB[i])
    Z, ZZ, DZ, DZZ = CALCDEPTH(Z, ZZ, DZ, DZZ, KB[i], KL1, KL2)

    # print('Z = ',Z,'\n')
    # print('ZZ = ',ZZ,'\n')
    # print('DZ = ',DZ,'\n')
    # print('DZZ = ',DZZ,'\n')

    if i == 0:
        plt.subplot(1, 2, 1)
        plt.plot(KB[i] * x, Z, 'b.', label='Z')
        plt.plot(KB[i] * x, ZZ, 'r.', label='ZZ')
        plt.legend(loc='best')
        plt.subplot(1, 2, 2)
        plt.plot(KB[i] * x, DZ, 'b.', label='DZ')
        plt.plot(KB[i] * x, DZZ, 'r.', label='DZZ')
        plt.legend(loc='best')
    else:
        plt.subplot(1, 2, 1)
        plt.plot(KB[i] * x, Z, 'b.', KB[i] * x, ZZ, 'r.')
        plt.subplot(1, 2, 2)
        plt.plot(KB[i] * x, DZ, 'b.', KB[i] * x, DZZ, 'r.')

plt.subplot(1, 2, 1)
plt.title('Vertical Layers')
plt.xlabel('# Layers')
plt.ylabel('$\sigma$')

plt.subplot(1, 2, 2)
plt.title('Spacing')
plt.xlabel('# Layers')
plt.ylabel('$\Delta$')

plt.show()

# Comments:
#   - Profiles are consistent for all KB values tested
#   - Second Z layer and first ZZ layer overlap rather than alternate
#   - Z should go to -1 but appears to stop one point short (potential solution is to add another
#       point to fill out profile, but Z & ZZ vectors will be of different lencths)
#   - All DZ and DZZ values are at least zero, indicating the current layer is always below the previous
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# DENS
#
# The first set of tests will examine the density profile along ZZ by changing the input values for
# KB and ZZ while holding T, S, and DT constant. The second set will change the T input while holding
# all other inputs constant. The third set will change the S input while holding all other inputs constant.
# Notes: All tests assume temperature and salinity decrease with depth.
#        Linear temperature and salinity profiles where used for all tests.
#        For tests 1 and 3, Tmax = 22, Tmin = 17. (C)
#        For tests 1 and 3, Smax = 37, Smin = 32. (psu)
#        Hold Case - KB = 31, Tmax = 22, Tmin = 17, Smax = 37, Smin = 32

# TEST 1
KB = [14, 21, 31, 61]
plt.figure(2)

for i in range(0, len(KB)):
    KL1 = np.floor(0.3 * KB[i])
    KL2 = KB[i] - 2
    x = np.ones(KB[i])

    Z = np.empty(KB[i])
    ZZ = np.empty(KB[i])
    DZ = np.empty(KB[i])
    DZZ = np.empty(KB[i])
    Z, ZZ, DZ, DZZ = CALCDEPTH(Z, ZZ, DZ, DZZ, KB[i], KL1, KL2)

    RHO = np.empty(KB[i])
    T = np.linspace(293,293-2*(KB[i]-1),KB[i])
    S = np.linspace(293,293-2*(KB[i]-1),KB[i])
    DT = 1
    RHO = DENS(T, S, ZZ, DT, RHO, KB[i])

    # print(RHO)
    plt.subplot(1,3,1)
    if i == 0:
        plt.plot(KB[i] * x, RHO, '.', label='KB=14')
    elif i == 1:
        plt.plot(KB[i] * x, RHO, '.', label='KB=21')
    elif i == 2:
        plt.plot(KB[i] * x, RHO, 'r.', label='KB=31')
    elif i == 3:
        plt.plot(KB[i] * x, RHO, '.', label='KB=61')
plt.subplot(1,3,1)
plt.title('Density Distribution, Varied KB & ZZ')
plt.xlabel('# Layers')
plt.ylabel('$\\rho$')
plt.legend(loc='best')


# TEST 2
KB = 31

KL1 = np.floor(0.3 * KB)
KL2 = KB - 2
x = np.ones(KB)

Z = np.empty(KB)
ZZ = np.empty(KB)
DZ = np.empty(KB)
DZZ = np.empty(KB)
Z, ZZ, DZ, DZZ = CALCDEPTH(Z, ZZ, DZ, DZZ, KB, KL1, KL2)

RHO = np.empty(KB)
T = np.linspace(22, 17, KB)
DT = 1

max_val = [31, 34, 37, 40]

for i in range(0,len(max_val)):
    S = np.linspace(max_val[i],max_val[i]-5, KB)
    RHO = DENS(T, S, ZZ, DT, RHO, KB)

    plt.subplot(1, 3, 2)
    if i == 0:
        plt.plot(max_val[i] * x, RHO, '.', label='S_max=31')
    elif i == 1:
        plt.plot(max_val[i] * x, RHO, '.', label='S_max=34')
    elif i == 2:
        plt.plot(max_val[i] * x, RHO, 'r.', label='S_max=37')
    elif i == 3:
        plt.plot(max_val[i] * x, RHO, '.', label='S_max=40')

plt.subplot(1, 3, 2)
plt.title('Density Distribution, Varied S')
plt.xlabel('S_max')
plt.ylabel('$\\rho$')
plt.legend(loc='best')

# TEST 3
KB = 31

KL1 = np.floor(0.3 * KB)
KL2 = KB - 2
x = np.ones(KB)

Z = np.empty(KB)
ZZ = np.empty(KB)
DZ = np.empty(KB)
DZZ = np.empty(KB)
Z, ZZ, DZ, DZZ = CALCDEPTH(Z, ZZ, DZ, DZZ, KB, KL1, KL2)

RHO = np.empty(KB)
S = np.linspace(37, 32, KB)
DT = 1

max_val = [16, 19, 22, 25]

for i in range(0,len(max_val)):
    T = np.linspace(max_val[i],max_val[i]-5, KB)
    RHO = DENS(T, S, ZZ, DT, RHO, KB)

    plt.subplot(1, 3, 3)
    if i == 0:
        plt.plot(max_val[i] * x, RHO, '.', label='T_max=16')
    elif i == 1:
        plt.plot(max_val[i] * x, RHO, '.', label='T_max=19')
    elif i == 2:
        plt.plot(max_val[i] * x, RHO, 'r.', label='T_max=22')
    elif i == 3:
        plt.plot(max_val[i] * x, RHO, '.', label='T_max=25')

plt.subplot(1,3,3)
plt.title('Density Distribution, Varied T & S')
plt.xlabel('T_max')
plt.ylabel('$\\rho$')
plt.legend(loc='best')

plt.show()

# Comments:
#   - Increased KB decreases the minimum density value (at the bottom), but the maximum value (at the surface)
#   - Indirect relationship between salinity and density magnitude. Spacing of density between layers unchanged.
#   - Direct relationship between temperature and density magnitude. Spacing of density between layers unchanged.
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MLDPTH
#
# The first set of tests will examine the mixed layer depth by changing the input values for KB and ZZ
# while holding T constant. The second set will change the input value of T while holding KB & ZZ constant.
# Note: Hold Case - KB = 31, Tmax = 22, Tmin = 17

# TEST 1
KB = [31, 46, 61, 76]

plt.figure(3)
for i in range(0, len(KB)):
    KL1 = np.floor(0.3 * KB[i])
    KL2 = KB[i] - 2
    x = np.ones(KB[i])

    Z = np.empty(KB[i])
    ZZ = np.empty(KB[i])
    DZ = np.empty(KB[i])
    DZZ = np.empty(KB[i])
    Z, ZZ, DZ, DZZ = CALCDEPTH(Z, ZZ, DZ, DZZ, KB[i], KL1, KL2)

    T = np.linspace(22,17,KB[i])
    ZZMLD = np.empty(KB[i])
    ZZMLD = MLDPTH(ZZ,T,KB[i],ZZMLD)
    print(ZZMLD)

    if i == 0:
        plt.subplot(1, 2, 1)
        plt.plot(KB[i] * x, ZZMLD, 'r.', label='KB=31')
    elif i == 1:
        plt.subplot(1, 2, 1)
        plt.plot(KB[i] * x, ZZMLD, '.', label='KB=46')
    elif i == 2:
        plt.subplot(1, 2, 1)
        plt.plot(KB[i] * x, ZZMLD, '.', label='KB=61')
    elif i == 3:
        plt.subplot(1, 2, 1)
        plt.plot(KB[i] * x, ZZMLD, '.', label='KB=76')


plt.subplot(1, 2, 1)
plt.title('Mixed Layer Depth, Varied KB&ZZ')
plt.xlabel('# Layers')
plt.ylabel('ZZMLD')
plt.legend(loc='best')

# TEST 2
KB = 31

KL1 = np.floor(0.3 * KB)
KL2 = KB - 2
x = np.ones(KB)

Z = np.empty(KB)
ZZ = np.empty(KB)
DZ = np.empty(KB)
DZZ = np.empty(KB)
Z, ZZ, DZ, DZZ = CALCDEPTH(Z, ZZ, DZ, DZZ, KB, KL1, KL2)

ZZMLD = np.empty(KB)
max_val = [16, 19, 22, 25]

for i in range(0,len(max_val)):
    T = np.linspace(max_val[i],max_val[i]-5, KB)
    ZZMLD = MLDPTH(ZZ,T,KB,ZZMLD)
    print(ZZMLD)

    plt.subplot(1, 2, 2)
    if i == 0:
        plt.plot(max_val[i] * x, ZZMLD, '.', label='T_max=16')
    elif i == 1:
        plt.plot(max_val[i] * x, ZZMLD, '.', label='T_max=19')
    elif i == 2:
        plt.plot(max_val[i] * x, ZZMLD, 'r.', label='T_max=22')
    elif i == 3:
        plt.plot(max_val[i] * x, ZZMLD, '.', label='T_max=25')


plt.subplot(1, 2, 2)
plt.title('Mixed Layer Depth, Varied T')
plt.xlabel('T_max')
plt.ylabel('ZZMLD')
plt.legend(loc='best')

plt.show()

# Comments:
#   - Function appears to break for KB = 46, 61, 76. Further tests will need to be run to verify wheteher
#       the issue is a result of KB or the range of temperatures used.
#   - Mixed Layer depth is independent of temperature.
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
