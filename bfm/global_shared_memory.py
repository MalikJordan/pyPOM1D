from include import INCLUDE_BEN, INCLUDE_SEAICE
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODULE: MEM
#
# DESCRIPTION: Definition of Global Shared Memory.
#
#              This module contains all the structural definitions of the BFM
#              and sets up the memory layout.
#              It is automatically generated from the prototype file
#              BFM/proto/ModuleMem.proto by including the information from
#              BFM/General/GlobalDefsBFM.model
#              Do not directly edit this code because changes will be lost at
#              any new compilation.
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   STATE VARIABLES INFO (PEL)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# !    3d name                                                  description           unit
# ! ---------- ------------------------------------------------------------ ---------------
# !        o2o                                                       Oxygen      mmol O2/m3
# !        n1p                                                    Phosphate       mmol P/m3
# !        n3n                                                      Nitrate       mmol N/m3
# !        n4n                                                     Ammonium       mmol N/m3
# !        o4n                                                 o4n       mmol N/m3
# !        n5s                                                     Silicate      mmol Si/m3
# !        n6r                                        Reduction Equivalents     mmol S--/m3
# !        b1c                               Aerobic and Anaerobic Bacteria         mg C/m3
# !        b1n                               Aerobic and Anaerobic Bacteria       mmol N/m3
# !        b1p                               Aerobic and Anaerobic Bacteria       mmol P/m3
# !        p1c                                                      Diatoms         mg C/m3
# !        p1n                                                      Diatoms       mmol N/m3
# !        p1p                                                      Diatoms       mmol P/m3
# !        p1l                                                      Diatoms       mg Chl/m3
# !        p1s                                                      Diatoms      mmol Si/m3
# !        p2c                                                  Flagellates         mg C/m3
# !        p2n                                                  Flagellates       mmol N/m3
# !        p2p                                                  Flagellates       mmol P/m3
# !        p2l                                                  Flagellates       mg Chl/m3
# !        p3c                                            PicoPhytoplankton         mg C/m3
# !        p3n                                            PicoPhytoplankton       mmol N/m3
# !        p3p                                            PicoPhytoplankton       mmol P/m3
# !        p3l                                            PicoPhytoplankton       mg Chl/m3
# !        p4c                                          Large Phytoplankton         mg C/m3
# !        p4n                                          Large Phytoplankton       mmol N/m3
# !        p4p                                          Large Phytoplankton       mmol P/m3
# !        p4l                                          Large Phytoplankton       mg Chl/m3
# !        z3c                                  Carnivorous Mesozooplankton         mg C/m3
# !        z3n                                  Carnivorous Mesozooplankton       mmol N/m3
# !        z3p                                  Carnivorous Mesozooplankton       mmol P/m3
# !        z4c                                   Omnivorous Mesozooplankton         mg C/m3
# !        z4n                                   Omnivorous Mesozooplankton       mmol N/m3
# !        z4p                                   Omnivorous Mesozooplankton       mmol P/m3
# !        z5c                                             Microzooplankton         mg C/m3
# !        z5n                                             Microzooplankton       mmol N/m3
# !        z5p                                             Microzooplankton       mmol P/m3
# !        z6c                         Heterotrophic Nanoflagellates (HNAN)         mg C/m3
# !        z6n                         Heterotrophic Nanoflagellates (HNAN)       mmol N/m3
# !        z6p                         Heterotrophic Nanoflagellates (HNAN)       mmol P/m3
# !        r1c                              Labile Dissolved Organic Matter         mg C/m3
# !        r1n                              Labile Dissolved Organic Matter       mmol N/m3
# !        r1p                              Labile Dissolved Organic Matter       mmol P/m3
# !        r2c                         Semi-labile Dissolved Organic Carbon         mg C/m3
# !        r3c                     Semi-refractory Dissolved Organic Carbon         mg C/m3
# !        r6c                                   Particulate Organic Matter         mg C/m3
# !        r6n                                   Particulate Organic Matter       mmol N/m3
# !        r6p                                   Particulate Organic Matter       mmol P/m3
# !        r6s                                   Particulate Organic Matter      mmol Si/m3
# !        o3c                                   Dissolved Inorganic Carbon         mg C/m3
# !        o3h                                   Dissolved Inorganic Carbon      mmol eq/m3

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   STATE VARIABLES INFO (ICE)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   DEFINITION OF ARRAYS WHICH WILL HOLD ALL STATE VARIABLES AND OTHER
#   GLOBAL VARIABLES USED FOR EXCHANGE BETWEEN SUBMODELS AND/OR OUTPUT
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL SYSTEM CONSTANTS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
iiPel   = 0
iiIce   = 700
iiBen   = 1000
iiReset = -1000

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL SYSTEM CONSTANTS (PEL)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
no_d3_box_states   = 50
no_d3_box_diagnoss = 91
no_d2_box_diagnoss = 162
no_d3_box_flux     = 30
#

# try:
#     INCLUDE_SEAICE
# except NameError:
#     INCLUDE_SEAICE = False
# else:
#     INCLUDE_SEAICE = True
if INCLUDE_SEAICE:
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   GLOBAL SYSTEM CONSTANTS (ICE)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    no_d2_box_states_ice   = 0
    no_d2_box_diagnoss_ice = 0
    no_d2_box_flux_ice     = 0

# try:
#     INCLUDE_BEN
# except NameError:
#     INCLUDE_BEN = False
# else:
#     INCLUDE_BEN = True
if INCLUDE_BEN:
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   GLOBAL SYSTEM CONSTANTS (BEN)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    no_d2_box_states_ben   = 0
    no_d2_box_diagnoss_ben = 0
    no_d2_box_flux_ben     = 0

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL DEFINITION OF PELAGIC (D3/D2) STATE VARIABLES
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
ppo2o = 1;  ppn1p = 2;  ppn3n = 3
ppn4n = 4;  ppo4n = 5;  ppn5s = 6;  ppn6r = 7;  ppb1c = 8;  ppb1n = 9;  ppb1p = 10; ppp1c = 11
ppp1n = 12; ppp1p = 13; ppp1l = 14; ppp1s = 15; ppp2c = 16; ppp2n = 17; ppp2p = 18
ppp2l = 19; ppP2s = 0;  ppp3c = 20; ppp3n = 21; ppp3p = 22; ppp3l = 23; ppP3s = 0
ppp4c = 24; ppp4n = 25; ppp4p = 26; ppp4l = 27; ppP4s = 0;  ppz3c = 28; ppz3n = 29
ppz3p = 30; ppz4c = 31; ppz4n = 32; ppz4p = 33; ppz5c = 34; ppz5n = 35; ppz5p = 36
ppz6c = 37; ppz6n = 38; ppz6p = 39; ppr1c = 40; ppr1n = 41; ppr1p = 42; ppR1s = 0
ppr2c = 43; ppR2n = 0;  ppR2p = 0;  ppR2s = 0;  ppr3c = 44; ppR3n = 0;  ppR3p = 0; ppR3s = 0
ppr6c = 45; ppr6n = 46; ppr6p = 47; ppr6s = 48; ppo3c = 49; ppo3h = 50

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL DEFINITION OF SEAICE (D2) STATE VARIABLES
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL DEFINITION OF SEABEN (D2) STATE VARIABLES
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   CONSTITUENT PARAMETERS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
iiC = 1; iiN = 2; iiP = 3
iiL = 4; iiS = 5; iiH = 6
iiLastElement = 6

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GROUP PARAMETERS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
iiPelBacteria      = 1; iiB1 = 1
iiPhytoPlankton    = 4; iiP1 = 1; iiP2 = 2; iiP3 = 3; iiP4 = 4
iiMesoZooPlankton  = 2; iiZ3 = 1; iiZ4 = 2
iiMicroZooPlankton = 2; iiZ5 = 1; iiZ6 = 2
iiPelDetritus      = 4; iiR1 = 1; iiR2 = 2; iiR3 = 3; iiR6 = 4
iiInorganic        = 1; iiO3 = 1

# CalcPelBacteria(iiPelBacteria)           = True
# CalcPhytoPlankton(iiPhytoPlankton)       = True
# CalcMesoZooPlankton(iiMesoZooPlankton)   = True
# CalcMicroZooPlankton(iiMicroZooPlankton) = True
# CalcPelDetritus(iiPelDetritus)           = True
# CalcInorganic(iiInorganic)               = True

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GROUP PARAMETERS (ICE)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GROUP PARAMETERS (BEN)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL VARIABLES (PEL)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL VARIABLES (ICE)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL VARIABLES (BEN)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL DEFINITION OF PELAGIC (D3/D2) STATE VARIABLES WHICH CAN BE OUTPUTTED IN NETCDF
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# !    3d name                                                  description           unit
# ! ---------- ------------------------------------------------------------ ---------------
# !        ETW                                                  temperature               C
# !        ESW                                                     Salinity               -
# !       ERHO                                             Seawater Density           kg/m3
# !        EIR                                                   Irradiance         uE/m2/s
# !        ESS                                          Suspended Sediments            g/m3
# !       exud                                                    exudation         mg C/m3
# !      Depth                                              Gridpoint Depth               m
# !     Volume                                             Gridpoint Volume              m3
# !       Area                                               Gridpoint Area              m2
# !        DIC                                   Dissolved Inorganic Carbon         umol/kg
# !        CO2                                                      CO2(aq)         umol/kg
# !       pCO2                                                 Oceanic pCO2            uatm
# !       HCO3                                                  Bicarbonate         umol/kg
# !        CO3                                                    Carbonate         umol/kg
# !        ALK                                                   Alkalinity      umol eq/kg
# !         pH                                                           pH               -
# !      OCalc                                  Saturation state of Calcite               -
# !      OArag                                Saturation state of Aragonite               -
# !        EPR                                               Water Pressure            dbar
# !    totpelc                                        Total Mass in Pelagic             g C
# !    totpeln                                        Total Mass in Pelagic             g N
# !    totpelp                                        Total Mass in Pelagic             g P
# !    totpels                                        Total Mass in Pelagic            g Si
# !      cxoO2                                            Oxygen Saturation      mmol O2/m3
# !     eO2mO2                                   Relative Oxygen saturation               -
# !       Chla                                                Chlorophyll-a       mg Chl/m3
# !    flPTn6r                        Pelagic Anaerobic Mineralization Rate    mmol O2/m3/d
# !    flN3o4n                                 Pelagic Denitrification Rate     mmol N/m3/d
# !    fln4n3n                                   Pelagic Nitrification Rate     mmol N/m3/d
# !     sediR2                                  Detritus sedimentation rate             m/d
# !     sediR6                                  Detritus sedimentation rate             m/d
# !       xEPS                                 Total Extinction Coefficient             1/m
# !   ABIO_eps                               Abiotic Extinction Coefficient             1/m
#
# ! qpcPPY(iiP1)                                                      Diatoms     mmol P/mg C
# ! qpcPPY(iiP2)                                                  Flagellates     mmol P/mg C
# ! qpcPPY(iiP3)                                            PicoPhytoplankton     mmol P/mg C
# ! qpcPPY(iiP4)                                          Large Phytoplankton     mmol P/mg C
# ! qncPPY(iiP1)                                                      Diatoms     mmol N/mg C
# ! qncPPY(iiP2)                                                  Flagellates     mmol N/mg C
# ! qncPPY(iiP3)                                            PicoPhytoplankton     mmol N/mg C
# ! qncPPY(iiP4)                                          Large Phytoplankton     mmol N/mg C
# ! qscPPY(iiP1)                                                      Diatoms    mmol Si/mg C
# ! qscPPY(iiP2)                                                  Flagellates    mmol Si/mg C
# ! qscPPY(iiP3)                                            PicoPhytoplankton    mmol Si/mg C
# ! qscPPY(iiP4)                                          Large Phytoplankton    mmol Si/mg C
# ! qlcPPY(iiP1)                                                      Diatoms    mg Chl /mg C
# ! qlcPPY(iiP2)                                                  Flagellates    mg Chl /mg C
# ! qlcPPY(iiP3)                                            PicoPhytoplankton    mg Chl /mg C
# ! qlcPPY(iiP4)                                          Large Phytoplankton    mg Chl /mg C
# ! qpcMEZ(iiZ3)                                  Carnivorous Mesozooplankton     mmol P/mg C
# ! qpcMEZ(iiZ4)                                   Omnivorous Mesozooplankton     mmol P/mg C
# ! qncMEZ(iiZ3)                                  Carnivorous Mesozooplankton     mmol N/mg C
# ! qncMEZ(iiZ4)                                   Omnivorous Mesozooplankton     mmol N/mg C
# ! qpcMIZ(iiZ5)                                             Microzooplankton     mmol P/mg C
# ! qpcMIZ(iiZ6)                         Heterotrophic Nanoflagellates (HNAN)     mmol P/mg C
# ! qncMIZ(iiZ5)                                             Microzooplankton     mmol N/mg C
# ! qncMIZ(iiZ6)                         Heterotrophic Nanoflagellates (HNAN)     mmol N/mg C
# ! qpcOMT(iiR1)                              Labile Dissolved Organic Matter     mmol N/mg C
# ! qpcOMT(iiR2)                         Semi-labile Dissolved Organic Carbon     mmol N/mg C
# ! qpcOMT(iiR3)                     Semi-refractory Dissolved Organic Carbon     mmol N/mg C
# ! qpcOMT(iiR6)                                   Particulate Organic Matter     mmol N/mg C
# ! qncOMT(iiR1)                              Labile Dissolved Organic Matter     mmol P/mg C
# ! qncOMT(iiR2)                         Semi-labile Dissolved Organic Carbon     mmol P/mg C
# ! qncOMT(iiR3)                     Semi-refractory Dissolved Organic Carbon     mmol P/mg C
# ! qncOMT(iiR6)                                   Particulate Organic Matter     mmol P/mg C
# ! qscOMT(iiR1)                              Labile Dissolved Organic Matter    mmol Si/mg C
# ! qscOMT(iiR2)                         Semi-labile Dissolved Organic Carbon    mmol Si/mg C
# ! qscOMT(iiR3)                     Semi-refractory Dissolved Organic Carbon    mmol Si/mg C
# ! qscOMT(iiR6)                                   Particulate Organic Matter    mmol Si/mg C
# ! qpcPBA(iiB1)                               Aerobic and Anaerobic Bacteria     mmol P/mg C
# ! qncPBA(iiB1)                               Aerobic and Anaerobic Bacteria     mmol N/mg C
# ! sediPPY(iiP1)                                                      Diatoms             m/d
# ! sediPPY(iiP2)                                                  Flagellates             m/d
# ! sediPPY(iiP3)                                            PicoPhytoplankton             m/d
# ! sediPPY(iiP4)                                          Large Phytoplankton             m/d
# ! sediMIZ(iiZ5)                                             Microzooplankton             m/d
# ! sediMIZ(iiZ6)                         Heterotrophic Nanoflagellates (HNAN)             m/d
# ! sediMEZ(iiZ3)                                  Carnivorous Mesozooplankton             m/d
# ! sediMEZ(iiZ4)                                   Omnivorous Mesozooplankton             m/d
# ! sunPPY(iiP1)                                                      Diatoms             1/d
# ! sunPPY(iiP2)                                                  Flagellates             1/d
# ! sunPPY(iiP3)                                            PicoPhytoplankton             1/d
# ! sunPPY(iiP4)                                          Large Phytoplankton             1/d
# ! eiPPY(iiP1)                                                      Diatoms               -
# ! eiPPY(iiP2)                                                  Flagellates               -
# ! eiPPY(iiP3)                                            PicoPhytoplankton               -
# ! eiPPY(iiP4)                                          Large Phytoplankton               -
# ! ELiPPY(iiP1)                                                      Diatoms            W/m2
# ! ELiPPY(iiP2)                                                  Flagellates            W/m2
# ! ELiPPY(iiP3)                                            PicoPhytoplankton            W/m2
# ! ELiPPY(iiP4)                                          Large Phytoplankton            W/m2

ppETW     = 1;  ppESW     = 2;  ppERHO    = 3;  ppEIR = 4
ppESS     = 5;  ppexud    = 6;  ppDepth   = 7;  ppVolume  = 8;  ppArea     = 9;  ppDIC   = 10; ppCO2   = 11
pppCO2    = 12; ppHCO3    = 13; ppCO3     = 14; ppALK     = 15; pppH       = 16; ppOCalc = 17; ppOArag = 18
ppEPR     = 19; pptotpelc = 20; pptotpeln = 21; pptotpelp = 22; pptotpels  = 23
ppcxoO2   = 24; ppeO2mO2  = 25; ppChla    = 26; ppflPTn6r = 27; ppflN3o4n  = 28
ppfln4n3n = 29; ppsediR2  = 30; ppsediR6  = 31; ppxEPS    = 32; ppABIO_eps = 33

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# !    2d name                                                  description           unit
# ! ---------- ------------------------------------------------------------ ---------------
# !      ETAUB                                                Bottom Stress           N m/s
# !   EPCO2air                             Atmospheric CO2 Partial Pressure            uatm
# ! CO2airflux                                             Sea-air CO2 Flux       mmol/m2/d
# !     Area2d                                           2-D Gridpoint Area              m2
# ! ThereIsLight                                   Switch for day/night cycle               -
# !       SUNQ                                           Daylength in hours               h
# !      EWIND                                                   Wind speed             m/s
# !    totsysc                                                   total mass             g C
# !    totsysn                                                   total mass             g N
# !    totsysp                                                   total mass             g P
# !    totsyss                                                   total mass            g Si
# !       EICE                                             Sea-ice fraction               -

ppETAUB   = 1; ppEPCO2air     = 2;  ppCO2airflux = 3
ppArea2d  = 4; ppThereIsLight = 5;  ppSUNQ       = 6;  ppEWIND = 7; pptotsysc = 8
pptotsysn = 9; pptotsysp      = 10; pptotsyss    = 11; ppEICE  = 12

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL DEFINITION OF SEAICE (D2) STATE VARIABLES WHICH CAN BE OUTPUTTED IN NETCDF
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL DEFINITION OF BENTHIC (D2) STATE VARIABLES WHICH CAN BE OUTPUTTED IN NETCDF
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   BOUNDARY FLUXES
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   OHTER 3D-GLOBAL VARIABLES
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   OHTER 2D-GLOBAL VARIABLES (PEL)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   OHTER 2D-GLOBAL VARIABLES (PEL)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   SHARED GLOBAL FUNCTIONS (PEL) (MUST BE BELOW CONTAINS)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# flux, flux_vector, Source, Source_D3_vector, \
#     fixed_quota_flux_vector
# ppPelBacteria, PelBacteria, ppPhytoPlankton, PhytoPlankton, \
#     ppMesoZooPlankton, MesoZooPlankton, ppMicroZooPlankton, MicroZooPlankton, \
#     ppPelDetritus, PelDetritus, ppInorganic, Inorganic



# try:
#     INCLUDE_SEAICE
# except NameError:
#     INCLUDE_SEAICE = False
# else:
#     INCLUDE_SEAICE = True
if INCLUDE_SEAICE:
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SHARED GLOBAL FUNCTIONS (ICE) (MUST BE BELOW CONTAINS)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Source_D2_vector_ice
    pass

# try:
#     INCLUDE_BEN
# except NameError:
#     INCLUDE_BEN = False
# else:
#     INCLUDE_BEN = True
if INCLUDE_BEN:
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SHARED GLOBAL FUNCTIONS (ICE) (MUST BE BELOW CONTAINS)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Source_D2_vector_ben
    pass

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GROUP PELAGIC (D3) STATE FUNCTIONS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# FUNCTIONS DEFINED IN MODULEMEM
def ppPelBacteria(n,constituent):

    pointers = [ppb1c, ppb1n, ppb1p, 0., 0., 0.]

    if n > 1 or n == 0:
        ppPelBacteria = 0
    elif constituent > 1 or constituent == 0:
        ppPelBacteria = 0
    else:
        ppPelBacteria = pointers[(n-1)*iiLastElement + constituent]

    return ppPelBacteria


def ppPhytoPlankton(n,constituent):

    pointers = [ppp1c, ppp1n, ppp1p, ppp1l, ppp1s, 0.,
                ppp2c, ppp2n, ppp2p, ppp2l, 0.   , 0.,
                ppp3c, ppp3n, ppp3p, ppp3l, 0.   , 0.,
                ppp4c, ppp4n, ppp4p, ppp4l, 0.   , 0.]

    if n > 4 or n == 0:
        ppPhytoPlankton = 0
    elif constituent > iiLastElement or constituent == 0:
        ppPhytoPlankton = 0
    else:
        ppPhytoPlankton = pointers[(n-1)*iiLastElement + constituent]

    return ppPhytoPlankton


def ppMesoZooPlankton(n,constituent):

    pointers = [ppz3c, ppz3n, ppz3p, 0., 0., 0.,
                ppz4c, ppz4n, ppz4p, 0., 0., 0.]

    if n > 2 or n == 0:
        ppMesoZooPlankton = 0
    elif constituent > iiLastElement or constituent == 0:
        ppMesoZooPlankton = 0
    else:
        ppMesoZooPlankton = pointers[(n-1)*iiLastElement + constituent]

    return ppMesoZooPlankton


def ppMicroZooPlankton(n,constituent):

    pointers = [ppz5c, ppz5n, ppz5p, 0., 0., 0.,
                ppz6c, ppz6n, ppz6p, 0., 0., 0.]

    if n > 2 or n == 0:
        ppMicroZooPlankton = 0
    elif constituent > iiLastElement or constituent == 0:
        ppMicroZooPlankton = 0
    else:
        ppMicroZooPlankton = pointers[(n-1)*iiLastElement + constituent]

    return ppMicroZooPlankton


def ppPelDetritus(n,constituent):

    pointers = [ppr1c, ppr1n, ppr1p, 0.   , 0.   , 0.,
                ppr2c, 0.   , 0.   , 0.   , 0.   , 0.,
                ppr3c, 0.   , 0.   , 0.   , 0.   , 0.,
                ppr6c, ppr6n, ppr6p, 0.   , ppr6s, 0.]

    if n > 4 or n == 0:
        ppPelDetritus = 0
    elif constituent > iiLastElement or constituent == 0:
        ppPelDetritus = 0
    else:
        ppPelDetritus = pointers[(n-1)*iiLastElement + constituent]

    return ppPelDetritus



