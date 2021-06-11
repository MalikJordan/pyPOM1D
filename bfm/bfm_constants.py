# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# MODULE: CONSTANTS
#
# DESCRIPTION: Full list of Fortran parameters  (comparable with Sesame constants)
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# THIS INITIALIZATION HAS TO BE DONE IN MODULEPARAM BECAUSE SOME COMPILERS
# DO NOT ALLOW THE INITIALIZATION OF CONSTANTS WITH INTRINSIC FUNCTIONS

MIN_VAL_EXPFUN = 0.
ECOLOGY = 1  # BASE TEMPERATURE FOR Q10
TRANSPORT = 2
zero_kelvin = -273.15
molar_gas_constant = 83.131  # GAS CONSTANT: bar mol^-1 deg-1
carbon_molecular_weight = 12.0  # MOLECULAR WEIGHT CARBON
nitrogen_molecular_weight = 14.0  # MOLECULAR WEIGHT NITROGEN
phosphorus_molecular_weight = 31.0  # MOLECULAR WEIGHT PHOSPHORUS
silica_molecular_weight = 28.0  # MOLECULAR WEIGHT SILICA
molecular_weight_conversion_factor = 0.217  # MOLECULAR WEIGHT CONVERSION FACTOR EINSTEIN->W
seconds_per_day = 86400.  # SECONDS IN DAY
days_per_second = 0.  # INVERSE OF SECONDS IN DAY
one_per_day = 1.  # RATE WHICH IS USED IN CASES WHERE IMPLICITLY ASSUMED

NO_BENTHOS = 0
BENTHIC_RETURN = 1
BENTHIC_BIO = 2
BENTHIC_FULL = 3
HOURS_PER_DAY = 24.  # HOURS IN DAY
SOLAR_RADIATION = 1368.  # SOLAR RADIATION CONSTANT
NODATA = 0
DAILYDATA = 1
CLOUDDATA = 3
ALLDATA = 4
ONLYDAILY = 5

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   INTEGER CONSTANTS USED IN THE BENTHIC NUTRIENT DYNAMICS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

DEFINE = -1000
DOUBLE_DEFINE = -1100
PARAMETER_DEFINE = -1200
EQUATION = 0
SDERIVATIVE = -2
DERIVATIVE = -1
INTEGRAL = 1
EXPONENTIAL_INTEGRAL = 11
SHIFT = 21
LINEAR_TERM = 1
CONSTANT_TERM = 0
QUADRATIC_TERM = 2
EXPONENTIAL_TERM = -1
ZERO_EXPONENTIAL_TERM = -2
BESSELI_EXP_TERM = -5
BESSELK_EXP_TERM = -6
ADD = 1000
INPUT_TERM = 6001
START_ADD_TERM = 6002
INPUT_ADD_TERM = 6000
RFLUX = 1
MASS = 2
AVERAGE = 3
PARAMETER = 4
STANDARD = 0
SET_CONTINUITY = 100
SET_LAYER_INTEGRAL = 200
SET_LAYER_INTEGRAL_UNTIL = 300
SET_BOUNDARY = 400
SET_DEPTH_INTEGRAL = 500
GET = 9000
LABDA_1 = 1
LABDA_2 = 2
COEFFICIENT = 3
COEFF2PARA = -1000
LAYERS = 4
DIFFUSION = 5
POROSITY = 6
ADSORPTION = 7
# INITIALIZE = 0
FLAG = 1
METHOD = 2
LAYER1 = 1
LAYER2 = 2
LAYER3 = 3
LAYER4 = 4
LAYER5 = 5
LAYER6 = 6
LAYER7 = 7
LAYER8 = 8
FOR_ALL_LAYERS = -1
NUMBER_OF_PROFILES = 12
NCOEFF = 22
NLAYER = 8
