from include import INCLUDE_BEN
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# MODULE: PARAM
#
# DESCRIPTION: List of global model parameters.
#              (global variables that can be changed during the model initialization
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ! Global Switches : turn on/off or choose model components
# ! NAME                          KIND    DESCRIPTION
# ! CalcPelagicFlag               logical Pelagic System
# ! CalcBenthicFlag               numeric Benthic system
# !                                       0 = No Benthic System
# !                                       The following are Not Yet Activated
# !                                       1 = Simple Benthic Return
# !                                       2 = Benthic organisms and intermediate
# !                                           complexity nutrient regeneration
# !                                       3 = Benthic organisms and full nutrient
# !                                           regeneration (early diagenesis)
# ! CalcTransportFlag             logical Compute Transport Term (when coupled
# !                                       with a OGCM)
# ! CalcConservationFlag          logical Mass Conservation Check
# ! CalcPhytoPlankton             logical Pelagic Phytoplankton (vector)
# ! CalcPelBacteria               logical Pelagic Bacteria (vector)
# ! CalcMesoZooPlankton           logical Mesozooplankton (vector)
# ! CalcMicroZooPlankton          logical Microzooplankton (vector)
# ! CalcPelChemistry              logical Pelagic Hydrochemical Processes
# ! AssignPelBenFluxesInBFMFlag   logical Benthic-pelagic fluxes are added to the
# !                                       time integration
# ! AssignAirPelFluxesInBFMFlag   logical Air-sea fluxes are added to the
# !                                       time integration
# ! ChlDynamicsFlag               numeric Choose the dynamics of Chl-a
# !                                       1 = diagnostic, optimal light property
# !                                           in phytoplankton
# !                                           (Ebenhoeh et al 1995, ERSEM-II)
# !                                       2 = state variable, constituent of
# !                                           phytoplankton
# ! LightPeriodFlag               numeric Choose the light averaging period
# !                                       1 = Instantanous irradiance
# !                                       2 = Daily average
# !                                       3 = Daylight average with explicit
# !                                           photoperiod
# ! LightLocationFlag             numeric Choose the parameterization of light
# !                                       location in the discrete grid
# !                                       1 = Light at the top of the cell
# !                                       2 = Light in the middle of the cell
# !                                       3 = Average Light in the cell
# ! check_fixed_quota             numeric Check whether zooplankton have fixed quota
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
CalcPelagicFlag = True
CalcBenthicFlag = 0  # Switch for Benthic system
CalcSeaiceFlag  = True  # Switch for Seaice system

CalcTransportFlag           = False
CalcConservationFlag        = True
CalcPelChemistry            = True
AssignPelBenFluxesInBFMFlag = True
AssignAirPelFluxesInBFMFlag = True

ChlDynamicsFlag   = 2
LightPeriodFlag   = 1
LightLocationFlag = 3
check_fixed_quota = 0

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ! Global Parameters : used throughout the model and not related
# !                     to a specific component
# ! NAME          UNIT          DESCRIPTION
# ! p_small      [-]           Smallest numeric value (the model "zero")
# ! slp0         [mbar]        Reference sea level pressure
# ! p_PAR        [-]           Fraction of Photosynthetically Available Radiation
# ! p_eps0       [1/m]         Background extinction coefficient
# ! p_epsESS     [m2/g]        Specific attenuation coefficient of
# !                            suspended sediments
# ! p_epsChla   [m2/mgChla]    Chla-specific extinction coefficient
# ! p_epsR6      [m2/mgC]      Specific attenuation coefficient of particulate
# !                            detritus
# ! p_pe_r1c     [-]           Fractional content of C in cytoplasm
# ! p_pe_r1n     [-]           Fractional content of N in cytoplasm
# ! p_pe_r1p     [-]           Fractional content of P in cytoplasm
# ! p_qro        [mmolHS-/     Stoichiometric coefficient for
# !               mmolO2]      anaerobic reactions
# ! p_qon_dentri [mmolO2/      Stoichiometric coefficient for
# !               mmolN]       denitrification
# ! p_qon_nitri  [mmolO2/      Stoichiometric coefficient for
# !               mmolN]       nitrification (3/2)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   PELAGIC MODEL PARAMETERS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
p_small      = 1.0E-20
slp0         = 1013.25
p_PAR        = 0.50
p_eps0       = 0.04
p_epsESS     = 0.04E-03
p_epsChla    = 0.03
p_epsR6      = 0.1E-03
p_pe_r1c     = 0.60
p_pe_r1n     = 0.72
p_pe_r1p     = 0.832
p_pe_R1s     = 0.06
p_qro        = 0.50
p_qon_dentri = 1.25
p_qon_nitri  = 1.5

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ! Benthic model parameters
# ! NAME          UNIT          DESCRIPTION
# ! p_sedlevels   [-]           Number of sigma levels for benthic nutrient
# ! p_poro0       [-]           Constant porosity for 0D and 1D runs
# ! p_InitSink    Logical       parameter to Initialize BenthicSInk var.
# ! p_q10diff     [-]           Temperature-dependency porewater diffusion
# ! p_clDxm       [m]           minimal value of D?.m for calculation of the alpha
# ! p_d_tot       [m]           Thickness of modelled benthic sediment layers
# ! p_clD1D2m     [m]           minimum distance between D1m and D2m
# ! p_d_tot_2     [m]           maximal Thickness of D2m
# ! p_sedsigma    [-]           Parameter for sigma level distribution
# ! p_poro        [-]           Sediment porosity
# ! p_p_ae        [-]           Adsorption coefficient
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# 0D-PARAMETERS
p_sedlevels = 20
p_sedsigma  = 2.0
p_d_tot     = 0.30
p_poro0     = 0.4

# 1D-PARAMETERS
# try:
#     INCLUDE_BEN
# except NameError:
#     INCLUDE_BEN = False
# else:
#     INCLUDE_BEN = True
if INCLUDE_BEN:
    p_InitSink = 100.0
    p_q10diff  = 1.49
    p_clDxm    = 0.001
    p_clD1D2m  = 0.01
    p_d_tot_2  = 0.35

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   SEAICE MODEL PARAMETERS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


