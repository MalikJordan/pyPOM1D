import numpy as np
from include import INCLUDE_BEN, INCLUDE_SEAICE

class PelagicIndex:
    def __init__(self, start, end, state_start, state_end, diag_start, diag_end, flux_start, flux_end, diag2d_start,
                 diag2d_end,surface_start, surface_end, bottom_start, bottom_end, river_start, river_end):
        start = start
        end = end
        state_start = state_start
        state_end = state_end
        diag_start = diag_start
        diag_end = diag_end
        flux_start = flux_start
        flux_end = flux_end
        diag2d_start = diag2d_start
        diag2d_end = diag2d_end
        surface_start = surface_start
        surface_end = surface_end
        bottom_start = bottom_start
        bottom_end = bottom_end
        river_start = river_start
        river_end = river_end

class SeaIceIndex:
    def __init__(self, start, end, state_start, state_end, diag2d_start, diag2d_end, flux2d_start, flux2d_end):
        start = start
        end = end
        state_start = state_start
        state_end = state_end
        diag2d_start = diag2d_start
        diag2d_end = diag2d_end
        flux2d_start = flux2d_start
        flux2d_end = flux2d_end

class BenthicIndex:
    def __init__(self, start, end, state_start, state_end, diag2d_start, diag2d_end, flux2d_start, flux2d_end):
        start = start
        end = end
        state_start = state_start
        state_end = state_end
        diag2d_start = diag2d_start
        diag2d_end = diag2d_end
        flux2d_start = flux2d_start
        flux2d_end = flux2d_end

class Index:
    def __init__(self, start, end):
        start = start
        end = end

def set_variable_info_bfm():
    variable_abbrev = np.empty(333, dtype=str)
    variable_names = np.empty(333, dtype=str)
    variable_units = np.empty(333, dtype=str)

    # 3D variables description:
    variable_abbrev[0] = "o2o"
    variable_names[0] = "Oxygen"
    variable_units[0] = "mmol O2/m3"

    variable_abbrev[1] = "n1p"
    variable_names[1] = "Phosphate"
    variable_units[1] = "mmol P/m3"

    variable_abbrev[2] = "n3n"
    variable_names[2] = "Nitrate"
    variable_units[2] = "mmol N/m3"

    variable_abbrev[3] = "n4n"
    variable_names[3] = "Ammonium"
    variable_units[3] = "mmol N/m3"

    variable_abbrev[4] = "o4n"
    variable_names[4] = "o4n"
    variable_units[4] = "mmol N/m3"

    variable_abbrev[5] = "n5s"
    variable_names[5] = "Silicate"
    variable_units[5] = "mmol Si/m3"

    variable_abbrev[6] = "n6r"
    variable_names[6] = "Reduction Equivalents"
    variable_units[6] = "mmol S--/m3"

    variable_abbrev[7] = "b1c"
    variable_names[7] = "Aerobic and Anaerobic Bacteria"
    variable_units[7] = "mg C/m3"

    variable_abbrev[8] = "b1n"
    variable_names[8] = "Aerobic and Anaerobic Bacteria"
    variable_units[8] = "mmol N/m3"

    variable_abbrev[9] = "b1p"
    variable_names[9] = "Aerobic and Anaerobic Bacteria"
    variable_units[9] = "mmol P/m3"

    variable_abbrev[10] = "p1c"
    variable_names[10] = "Diatoms"
    variable_units[10] = "mg C/m3"

    variable_abbrev[11] = "p1n"
    variable_names[11] = "Diatoms"
    variable_units[11] = "mmol N/m3"

    variable_abbrev[12] = "p1p"
    variable_names[12] = "Diatoms"
    variable_units[12] = "mmol P/m3"

    variable_abbrev[13] = "p1l"
    variable_names[13] = "Diatoms"
    variable_units[13] = "mg Chl/m3"

    variable_abbrev[14] = "p1s"
    variable_names[14] = "Diatoms"
    variable_units[14] = "mmol Si/m3"

    variable_abbrev[15] = "p2c"
    variable_names[15] = "Flagellates"
    variable_units[15] = "mg C/m3"

    variable_abbrev[16] = "p2n"
    variable_names[16] = "Flagellates"
    variable_units[16] = "mmol N/m3"

    variable_abbrev[17] = "p2p"
    variable_names[17] = "Flagellates"
    variable_units[17] = "mmol P/m3"

    variable_abbrev[18] = "p2l"
    variable_names[18] = "Flagellates"
    variable_units[18] = "mg Chl/m3"

    variable_abbrev[19] = "p3c"
    variable_names[19] = "PicoPhytoplankton"
    variable_units[19] = "mg C/m3"

    variable_abbrev[20] = "p3n"
    variable_names[20] = "PicoPhytoplankton"
    variable_units[20] = "mmol N/m3"

    variable_abbrev[21] = "p3p"
    variable_names[21] = "PicoPhytoplankton"
    variable_units[21] = "mmol P/m3"

    variable_abbrev[22] = "p3l"
    variable_names[22] = "PicoPhytoplankton"
    variable_units[22] = "mg Chl/m3"

    variable_abbrev[23] = "p4c"
    variable_names[23] = "Large Phytoplankton"
    variable_units[23] = "mg C/m3"

    variable_abbrev[24] = "p4n"
    variable_names[24] = "Large Phytoplankton"
    variable_units[24] = "mmol N/m3"

    variable_abbrev[25] = "p4p"
    variable_names[25] = "Large Phytoplankton"
    variable_units[25] = "mmol P/m3"

    variable_abbrev[26] = "p4l"
    variable_names[26] = "Large Phytoplankton"
    variable_units[26] = "mg Chl/m3"

    variable_abbrev[27] = "z3c"
    variable_names[27] = "Carnivorous Mesozooplankton"
    variable_units[27] = "mg C/m3"

    variable_abbrev[28] = "z3n"
    variable_names[28] = "Carnivorous Mesozooplankton"
    variable_units[28] = "mmol N/m3"

    variable_abbrev[29] = "z3p"
    variable_names[29] = "Carnivorous Mesozooplankton"
    variable_units[29] = "mmol P/m3"

    variable_abbrev[30] = "z4c"
    variable_names[30] = "Omnivorous Mesozooplankton"
    variable_units[30] = "mg C/m3"

    variable_abbrev[31] = "z4n"
    variable_names[31] = "Omnivorous Mesozooplankton"
    variable_units[31] = "mmol N/m3"

    variable_abbrev[32] = "z4p"
    variable_names[32] = "Omnivorous Mesozooplankton"
    variable_units[32] = "mmol P/m3"

    variable_abbrev[33] = "z5c"
    variable_names[33] = "Microzooplankton"
    variable_units[33] = "mg C/m3"

    variable_abbrev[34] = "z5n"
    variable_names[34] = "Microzooplankton"
    variable_units[34] = "mmol N/m3"

    variable_abbrev[35] = "z5p"
    variable_names[35] = "Microzooplankton"
    variable_units[35] = "mmol P/m3"

    variable_abbrev[36] = "z6c"
    variable_names[36] = "Heterotrophic Nanoflagellates (HNAN)"
    variable_units[36] = "mg C/m3"

    variable_abbrev[37] = "z6n"
    variable_names[37] = "Heterotrophic Nanoflagellates (HNAN)"
    variable_units[37] = "mmol N/m3"

    variable_abbrev[38] = "z6p"
    variable_names[38] = "Heterotrophic Nanoflagellates (HNAN)"
    variable_units[38] = "mmol P/m3"

    variable_abbrev[39] = "r1c"
    variable_names[39] = "Labile Dissolved Organic Matter"
    variable_units[39] = "mg C/m3"

    variable_abbrev[40] = "r1n"
    variable_names[40] = "Labile Dissolved Organic Matter"
    variable_units[40] = "mmol N/m3"

    variable_abbrev[41] = "r1p"
    variable_names[41] = "Labile Dissolved Organic Matter"
    variable_units[41] = "mmol P/m3"

    variable_abbrev[42] = "r2c"
    variable_names[42] = "Semi-labile Dissolved Organic Carbon"
    variable_units[42] = "mg C/m3"

    variable_abbrev[43] = "r3c"
    variable_names[43] = "Semi-refractory Dissolved Organic Carbon"
    variable_units[43] = "mg C/m3"

    variable_abbrev[44] = "r6c"
    variable_names[44] = "Particulate Organic Matter"
    variable_units[44] = "mg C/m3"

    variable_abbrev[45] = "r6n"
    variable_names[45] = "Particulate Organic Matter"
    variable_units[45] = "mmol N/m3"

    variable_abbrev[46] = "r6p"
    variable_names[46] = "Particulate Organic Matter"
    variable_units[46] = "mmol P/m3"

    variable_abbrev[47] = "r6s"
    variable_names[47] = "Particulate Organic Matter"
    variable_units[47] = "mmol Si/m3"

    variable_abbrev[48] = "o3c"
    variable_names[48] = "Dissolved Inorganic Carbon"
    variable_units[48] = "mg C/m3"

    variable_abbrev[49] = "o3h"
    variable_names[49] = "Dissolved Inorganic Carbon"
    variable_units[49] = "mmol eq/m3"

    variable_abbrev[50] = "ETW"
    variable_names[50] = "temperature"
    variable_units[50] = "C"

    variable_abbrev[51] = "ESW"
    variable_names[51] = "Salinity"
    variable_units[51] = "-"

    variable_abbrev[52] = "ERHO"
    variable_names[52] = "Seawater Density"
    variable_units[52] = "kg/m3"

    variable_abbrev[53] = "EIR"
    variable_names[53] = "Irradiance"
    variable_units[53] = "uE/m2/s"

    variable_abbrev[54] = "ESS"
    variable_names[54] = "Suspended Sediments"
    variable_units[54] = "g/m3"

    variable_abbrev[55] = "exud"
    variable_names[55] = "exudation"
    variable_units[55] = "mg C/m3"

    variable_abbrev[56] = "Depth"
    variable_names[56] = "Gridpoint Depth"
    variable_units[56] = "m"

    variable_abbrev[57] = "Volume"
    variable_names[57] = "Gridpoint Volume"
    variable_units[57] = "m3"

    variable_abbrev[58] = "Area"
    variable_names[58] = "Gridpoint Area"
    variable_units[58] = "m2"

    variable_abbrev[59] = "DIC"
    variable_names[59] = "Dissolved Inorganic Carbon"
    variable_units[59] = "umol/kg"

    variable_abbrev[60] = "CO2"
    variable_names[60] = "CO2(aq)"
    variable_units[60] = "umol/kg"

    variable_abbrev[61] = "pCO2"
    variable_names[61] = "Oceanic pCO2"
    variable_units[61] = "uatm"

    variable_abbrev[62] = "HCO3"
    variable_names[62] = "Bicarbonate"
    variable_units[62] = "umol/kg"

    variable_abbrev[63] = "CO3"
    variable_names[63] = "Carbonate"
    variable_units[63] = "umol/kg"

    variable_abbrev[64] = "ALK"
    variable_names[64] = "Alkalinity"
    variable_units[64] = "umol eq/kg"

    variable_abbrev[65] = "pH"
    variable_names[65] = "pH"
    variable_units[65] = "-"

    variable_abbrev[66] = "OCalc"
    variable_names[66] = "Saturation state of Calcite"
    variable_units[66] = "-"

    variable_abbrev[67] = "OArag"
    variable_names[67] = "Saturation state of Aragonite"
    variable_units[67] = "-"

    variable_abbrev[68] = "EPR"
    variable_names[68] = "Water Pressure"
    variable_units[68] = "dbar"

    variable_abbrev[69] = "totpelc"
    variable_names[69] = "Total Mass in Pelagic"
    variable_units[69] = "g C"

    variable_abbrev[70] = "totpeln"
    variable_names[70] = "Total Mass in Pelagic"
    variable_units[70] = "g N"

    variable_abbrev[71] = "totpelp"
    variable_names[71] = "Total Mass in Pelagic"
    variable_units[71] = "g P"

    variable_abbrev[72] = "totpels"
    variable_names[72] = "Total Mass in Pelagic"
    variable_units[72] = "g Si"

    variable_abbrev[73] = "cxoO2"
    variable_names[73] = "Oxygen Saturation"
    variable_units[73] = "mmol O2/m3"

    variable_abbrev[74] = "eO2mO2"
    variable_names[74] = "Relative Oxygen saturation"
    variable_units[74] = "-"

    variable_abbrev[75] = "Chla"
    variable_names[75] = "Chlorophyll-a"
    variable_units[75] = "mg Chl/m3"

    variable_abbrev[76] = "flPTn6r"
    variable_names[76] = "Pelagic Anaerobic Mineralization Rate"
    variable_units[76] = "mmol O2/m3/d"

    variable_abbrev[77] = "flN3o4n"
    variable_names[77] = "Pelagic Denitrification Rate"
    variable_units[77] = "mmol N/m3/d"

    variable_abbrev[78] = "fln4n3n"
    variable_names[78] = "Pelagic Nitrification Rate"
    variable_units[78] = "mmol N/m3/d"

    variable_abbrev[79] = "sediR2"
    variable_names[79] = "Detritus sedimentation rate"
    variable_units[79] = "m/d"

    variable_abbrev[80] = "sediR6"
    variable_names[80] = "Detritus sedimentation rate"
    variable_units[80] = "m/d"

    variable_abbrev[81] = "xEPS"
    variable_names[81] = "Total Extinction Coefficient"
    variable_units[81] = "1/m"

    variable_abbrev[82] = "ABIO_eps"
    variable_names[82] = "Abiotic Extinction Coefficient"
    variable_units[82] = "1/m"

    variable_abbrev[83] = "qpcPPY(iiP1)"
    variable_names[83] = "Quotum P/C in Phytoplankton"
    variable_units[83] = "mmol P/mg C"

    variable_abbrev[84] = "qpcPPY(iiP2)"
    variable_names[84] = "Quotum P/C in Phytoplankton"
    variable_units[84] = "mmol P/mg C"

    variable_abbrev[85] = "qpcPPY(iiP3)"
    variable_names[85] = "Quotum P/C in Phytoplankton"
    variable_units[85] = "mmol P/mg C"

    variable_abbrev[86] = "qpcPPY(iiP4)"
    variable_names[86] = "Quotum P/C in Phytoplankton"
    variable_units[86] = "mmol P/mg C"

    variable_abbrev[87] = "qncPPY(iiP1)"
    variable_names[87] = "Quotum N/C in Phytoplankton"
    variable_units[87] = "mmol N/mg C"

    variable_abbrev[88] = "qncPPY(iiP2)"
    variable_names[88] = "Quotum N/C in Phytoplankton"
    variable_units[88] = "mmol N/mg C"

    variable_abbrev[89] = "qncPPY(iiP3)"
    variable_names[89] = "Quotum N/C in Phytoplankton"
    variable_units[89] = "mmol N/mg C"

    variable_abbrev[90] = "qncPPY(iiP4)"
    variable_names[90] = "Quotum N/C in Phytoplankton"
    variable_units[90] = "mmol N/mg C"

    variable_abbrev[91] = "qscPPY(iiP1)"
    variable_names[91] = "Quotum Si/C in Phytoplankton"
    variable_units[91] = "mmol Si/mg C"

    variable_abbrev[92] = "qscPPY(iiP2)"
    variable_names[92] = "Quotum Si/C in Phytoplankton"
    variable_units[92] = "mmol Si/mg C"

    variable_abbrev[93] = "qscPPY(iiP3)"
    variable_names[93] = "Quotum Si/C in Phytoplankton"
    variable_units[93] = "mmol Si/mg C"

    variable_abbrev[94] = "qscPPY(iiP4)"
    variable_names[94] = "Quotum Si/C in Phytoplankton"
    variable_units[94] = "mmol Si/mg C"

    variable_abbrev[95] = "qlcPPY(iiP1)"
    variable_names[95] = "Quotum Chl/C in Phytoplankton"
    variable_units[95] = "mg Chl /mg C"

    variable_abbrev[96] = "qlcPPY(iiP2)"
    variable_names[96] = "Quotum Chl/C in Phytoplankton"
    variable_units[96] = "mg Chl /mg C"

    variable_abbrev[97] = "qlcPPY(iiP3)"
    variable_names[97] = "Quotum Chl/C in Phytoplankton"
    variable_units[97] = "mg Chl /mg C"

    variable_abbrev[98] = "qlcPPY(iiP4)"
    variable_names[98] = "Quotum Chl/C in Phytoplankton"
    variable_units[98] = "mg Chl /mg C"

    variable_abbrev[99] = "qpcMEZ(iiZ3)"
    variable_names[99] = "Quotum P/C in Mesozooplankton"
    variable_units[99] = "mmol P/mg C"

    variable_abbrev[100] = "qpcMEZ(iiZ4)"
    variable_names[100] = "Quotum P/C in Mesozooplankton"
    variable_units[100] = "mmol P/mg C"

    variable_abbrev[101] = "qncMEZ(iiZ3)"
    variable_names[101] = "Quotum N/C in Mesozooplankton"
    variable_units[101] = "mmol N/mg C"

    variable_abbrev[102] = "qncMEZ(iiZ4)"
    variable_names[102] = "Quotum N/C in Mesozooplankton"
    variable_units[102] = "mmol N/mg C"

    variable_abbrev[103] = "qpcMIZ(iiZ5)"
    variable_names[103] = "Quotum P/C in Z5(MicroZooPlankton)"
    variable_units[103] = "mmol P/mg C"

    variable_abbrev[104] = "qpcMIZ(iiZ6)"
    variable_names[104] = "Quotum P/C in Z6(MicroZooPlankton)"
    variable_units[104] = "mmol P/mg C"

    variable_abbrev[105] = "qncMIZ(iiZ5)"
    variable_names[105] = "Quotum N/C in Z5(MicroZooPlankton)"
    variable_units[105] = "mmol N/mg C"

    variable_abbrev[106] = "qncMIZ(iiZ6)"
    variable_names[106] = "Quotum N/C in Z6(MicroZooPlankton)"
    variable_units[106] = "mmol N/mg C"

    variable_abbrev[107] = "qpcOMT(iiR1)"
    variable_names[107] = "Quotum P/C in Organic Matter"
    variable_units[107] = "mmol N/mg C"

    variable_abbrev[108] = "qpcOMT(iiR2)"
    variable_names[108] = "Quotum P/C in Organic Matter"
    variable_units[108] = "mmol N/mg C"

    variable_abbrev[109] = "qpcOMT(iiR3)"
    variable_names[109] = "Quotum P/C in Organic Matter"
    variable_units[109] = "mmol N/mg C"

    variable_abbrev[110] = "qpcOMT(iiR6)"
    variable_names[110] = "Quotum P/C in Organic Matter"
    variable_units[110] = "mmol N/mg C"

    variable_abbrev[111] = "qncOMT(iiR1)"
    variable_names[111] = "Quotum N/C in Organic Matter"
    variable_units[111] = "mmol P/mg C"

    variable_abbrev[112] = "qncOMT(iiR2)"
    variable_names[112] = "Quotum N/C in Organic Matter"
    variable_units[112] = "mmol P/mg C"

    variable_abbrev[113] = "qncOMT(iiR3)"
    variable_names[113] = "Quotum N/C in Organic Matter"
    variable_units[113] = "mmol P/mg C"

    variable_abbrev[114] = "qncOMT(iiR6)"
    variable_names[114] = "Quotum N/C in Organic Matter"
    variable_units[114] = "mmol P/mg C"

    variable_abbrev[115] = "qscOMT(iiR1)"
    variable_names[115] = "Quotum Si/C in Organic Matter"
    variable_units[115] = "mmol Si/mg C"

    variable_abbrev[116] = "qscOMT(iiR2)"
    variable_names[116] = "Quotum Si/C in Organic Matter"
    variable_units[116] = "mmol Si/mg C"

    variable_abbrev[117] = "qscOMT(iiR3)"
    variable_names[117] = "Quotum Si/C in Organic Matter"
    variable_units[117] = "mmol Si/mg C"

    variable_abbrev[118] = "qscOMT(iiR6)"
    variable_names[118] = "Quotum Si/C in Organic Matter"
    variable_units[118] = "mmol Si/mg C"

    variable_abbrev[119] = "qpcPBA(iiB1)"
    variable_names[119] = "Quotum P/C in Pelagic Bacteria"
    variable_units[119] = "mmol P/mg C"

    variable_abbrev[120] = "qncPBA(iiB1)"
    variable_names[120] = "Quotum N/C in Pelagic Bacteria"
    variable_units[120] = "mmol N/mg C"

    variable_abbrev[121] = "sediPPY(iiP1)"
    variable_names[121] = "P1(PhytoPlankton) sedimentation rate"
    variable_units[121] = "m/d"

    variable_abbrev[122] = "sediPPY(iiP2)"
    variable_names[122] = "P2(PhytoPlankton) sedimentation rate"
    variable_units[122] = "m/d"

    variable_abbrev[123] = "sediPPY(iiP3)"
    variable_names[123] = "P3(PhytoPlankton) sedimentation rate"
    variable_units[123] = "m/d"

    variable_abbrev[124] = "sediPPY(iiP4)"
    variable_names[124] = "P4(PhytoPlankton) sedimentation rate"
    variable_units[124] = "m/d"

    variable_abbrev[125] = "sediMIZ(iiZ5)"
    variable_names[125] = "Z5(MicroZooPlankton) sedimentation rate"
    variable_units[125] = "m/d"

    variable_abbrev[126] = "sediMIZ(iiZ6)"
    variable_names[126] = "Z6(MicroZooPlankton) sedimentation rate"
    variable_units[126] = "m/d"

    variable_abbrev[127] = "sediMEZ(iiZ3)"
    variable_names[127] = "Z3(MesoZooPlankton) sedimentation rate"
    variable_units[127] = "m/d"

    variable_abbrev[128] = "sediMEZ(iiZ4)"
    variable_names[128] = "Z4(MesoZooPlankton) sedimentation rate"
    variable_units[128] = "m/d"

    variable_abbrev[129] = "sunPPY(iiP1)"
    variable_names[129] = "Specific Net Production of P1(PhytoPlankton)"
    variable_units[129] = "1/d"

    variable_abbrev[130] = "sunPPY(iiP2)"
    variable_names[130] = "Specific Net Production of P2(PhytoPlankton)"
    variable_units[130] = "1/d"

    variable_abbrev[131] = "sunPPY(iiP3)"
    variable_names[131] = "Specific Net Production of P3(PhytoPlankton)"
    variable_units[131] = "1/d"

    variable_abbrev[132] = "sunPPY(iiP4)"
    variable_names[132] = "Specific Net Production of P4(PhytoPlankton)"
    variable_units[132] = "1/d"

    variable_abbrev[133] = "eiPPY(iiP1)"
    variable_names[133] = "Regulating Factor for Light in P1(PhytoPlankton)"
    variable_units[133] = "-"

    variable_abbrev[134] = "eiPPY(iiP2)"
    variable_names[134] = "Regulating Factor for Light in P2(PhytoPlankton)"
    variable_units[134] = "-"

    variable_abbrev[135] = "eiPPY(iiP3)"
    variable_names[135] = "Regulating Factor for Light in P3(PhytoPlankton)"
    variable_units[135] = "-"

    variable_abbrev[136] = "eiPPY(iiP4)"
    variable_names[136] = "Regulating Factor for Light in P4(PhytoPlankton)"
    variable_units[136] = "-"

    variable_abbrev[137] = "ELiPPY(iiP1)"
    variable_names[137] = "Optimal light in P1(PhytoPlankton)"
    variable_units[137] = "W/m2"

    variable_abbrev[138] = "ELiPPY(iiP2)"
    variable_names[138] = "Optimal light in P2(PhytoPlankton)"
    variable_units[138] = "W/m2"

    variable_abbrev[139] = "ELiPPY(iiP3)"
    variable_names[139] = "Optimal light in P3(PhytoPlankton)"
    variable_units[139] = "W/m2"

    variable_abbrev[140] = "ELiPPY(iiP4)"
    variable_names[140] = "Optimal light in P4(PhytoPlankton)"
    variable_units[140] = "W/m2"

    variable_abbrev[141] = "fP1z4c"
    variable_names[141] = "diatom grazing by omniv.zooplankton"
    variable_units[141] = "mg C/m3/d"

    variable_abbrev[142] = "fP2z4c"
    variable_names[142] = "flagellates grazing by omniv.zooplankton"
    variable_units[142] = "mg C/m3/d"

    variable_abbrev[143] = "fP3z4c"
    variable_names[143] = "picophytoplankton grazing by omniv.zooplankton"
    variable_units[143] = "mg C/m3/d"

    variable_abbrev[144] = "fP4z4c"
    variable_names[144] = "large phytoplankton grazing by omniv.zooplankton"
    variable_units[144] = "mg C/m3/d"

    variable_abbrev[145] = "fP1z5c"
    variable_names[145] = "diatom grazing by microzooplankton"
    variable_units[145] = "mg C/m3/d"

    variable_abbrev[146] = "fP2z5c"
    variable_names[146] = "flagellates grazing by microzooplankton"
    variable_units[146] = "mg C/m3/d"

    variable_abbrev[147] = "fP3z5c"
    variable_names[147] = "picophytoplankton grazing by microzooplankton"
    variable_units[147] = "mg C/m3/d"

    variable_abbrev[148] = "fP4z5c"
    variable_names[148] = "large phytoplankton grazing by microzooplankton"
    variable_units[148] = "mg C/m3/d"

    variable_abbrev[149] = "fP1z6c"
    variable_names[149] = "diatom grazing by heterotrophic nanoflagellates"
    variable_units[149] = "mg C/m3/d"

    variable_abbrev[150] = "fP2z6c"
    variable_names[150] = "flagellates grazing by heterotrophic nanoflagellates"
    variable_units[150] = "mg C/m3/d"

    variable_abbrev[151] = "fP3z6c"
    variable_names[151] = "picophytoplankton grazing by heterotrophic nanoflagellates"
    variable_units[151] = "mg C/m3/d"

    variable_abbrev[152] = "fP4z6c"
    variable_names[152] = "large phytoplankton grazing by heterotrophic nanoflagellates"
    variable_units[152] = "mg C/m3/d"

    variable_abbrev[153] = "fB1z6c"
    variable_names[153] = "bacterial grazing by heterotrophic nanoflagellates"
    variable_units[153] = "mg C/m3/d"

    variable_abbrev[154] = "fB1z5c"
    variable_names[154] = "bacterial grazing by microzooplankton"
    variable_units[154] = "mg C/m3/d"

    variable_abbrev[155] = "ruPTc"
    variable_names[155] = "Gross Primary Production"
    variable_units[155] = "mg C/m3/d"

    variable_abbrev[156] = "resPP"
    variable_names[156] = "Respiration of phytoplankton"
    variable_units[156] = "mg C/m3/d"

    variable_abbrev[157] = "resZT"
    variable_names[157] = "Respiration of zooplankton"
    variable_units[157] = "mg C/m3/d"

    variable_abbrev[158] = "ruPTn"
    variable_names[158] = "net nutrient uptake"
    variable_units[158] = "mmoln /m3/d"

    variable_abbrev[159] = "ruPTp"
    variable_names[159] = "net phosphate uptake"
    variable_units[159] = "mmol P /m3/d"

    variable_abbrev[160] = "exPP"
    variable_names[160] = "C excretion from phytoplankton"
    variable_units[160] = "mg C/m3/d"

    variable_abbrev[161] = "ruZTc"
    variable_names[161] = "gross secondary production"
    variable_units[161] = "mg C/m3/d"

    variable_abbrev[162] = "netZTc"
    variable_names[162] = "net secondary production"
    variable_units[162] = "mg C/m3/d"

    variable_abbrev[163] = "rrPTo"
    variable_names[163] = "pelagic respiration"
    variable_units[163] = "mmol O2/m3/d"

    variable_abbrev[164] = "rePTn"
    variable_names[164] = "pelagic mineralization"
    variable_units[164] = "mmol N/m3/d"

    variable_abbrev[165] = "rePTp"
    variable_names[165] = "pelagic mineralization"
    variable_units[165] = "mmol P /m3/d"

    variable_abbrev[166] = "reBn"
    variable_names[166] = "bacterial mineralization"
    variable_units[166] = "mmol N/m3/d"

    variable_abbrev[167] = "ruBn"
    variable_names[167] = "bacterial uptake"
    variable_units[167] = "mmol N/m3/d"

    variable_abbrev[168] = "reBp"
    variable_names[168] = "bacterial mineralization"
    variable_units[168] = "mmol P/m3/d"

    variable_abbrev[169] = "ruBp"
    variable_names[169] = "bacterial uptake"
    variable_units[169] = "mmol P/m3/d"

    variable_abbrev[170] = "fR2b1c"
    variable_names[170] = "TEP uptake by bacteria"
    variable_units[170] = "mg C/m3/d"

    # 2D Pelagic variables description
    variable_abbrev[171] = "ETAUB"
    variable_names[171] = "Bottom Stress"
    variable_units[171] = "N m/s"

    variable_abbrev[172] = "EPCO2air"
    variable_names[172] = "Atmospheric CO2 Partial Pressure"
    variable_units[172] = "uatm"

    variable_abbrev[173] = "CO2airflux"
    variable_names[173] = "Sea-air CO2 Flux"
    variable_units[173] = "mmol/m2/d"

    variable_abbrev[174] = "Area2d"
    variable_names[174] = "2-D Gridpoint Area"
    variable_units[174] = "m2"

    variable_abbrev[175] = "ThereIsLight"
    variable_names[175] = "Switch for day/night cycle"
    variable_units[175] = "-"

    variable_abbrev[176] = "SUNQ"
    variable_names[176] = "Daylength in hours"
    variable_units[176] = "h"

    variable_abbrev[177] = "EWIND"
    variable_names[177] = "Wind speed"
    variable_units[177] = "m/s"

    variable_abbrev[178] = "totsysc"
    variable_names[178] = "total mass"
    variable_units[178] = "g C"

    variable_abbrev[179] = "totsysn"
    variable_names[179] = "total mass"
    variable_units[179] = "g N"

    variable_abbrev[180] = "totsysp"
    variable_names[180] = "total mass"
    variable_units[180] = "g P"

    variable_abbrev[181] = "totsyss"
    variable_names[181] = "total mass"
    variable_units[181] = "g Si"

    variable_abbrev[182] = "EICE"
    variable_names[182] = "Sea-ice fraction"
    variable_units[182] = "-"

    # Surface
    variable_abbrev[183] = "jsuro2o"
    variable_names[183] = "flux of Oxygen at SURFACE"
    variable_units[183] = "mmol O2/m2/day"

    variable_abbrev[184] = "jsurn1p"
    variable_names[184] = "flux of Phosphate at SURFACE"
    variable_units[184] = "mmol P/m2/day"

    variable_abbrev[185] = "jsurn3n"
    variable_names[185] = "flux of Nitrate at SURFACE"
    variable_units[185] = "mmol N/m2/day"

    variable_abbrev[186] = "jsurn4n"
    variable_names[186] = "flux of Ammonium at SURFACE"
    variable_units[186] = "mmol N/m2/day"

    variable_abbrev[187] = "jsuro4n"
    variable_names[187] = "flux of o4n at SURFACE"
    variable_units[187] = "mmol N/m2/day"

    variable_abbrev[188] = "jsurn5s"
    variable_names[188] = "flux of Silicate at SURFACE"
    variable_units[188] = "mmol Si/m2/day"

    variable_abbrev[189] = "jsurn6r"
    variable_names[189] = "flux of Reduction Equivalents at SURFACE"
    variable_units[189] = "mmol S--/m2/day"

    variable_abbrev[190] = "jsurb1c"
    variable_names[190] = "flux of Aerobic and Anaerobic Bacteria at SURFACE"
    variable_units[190] = "mg C/m2/day"

    variable_abbrev[191] = "jsurb1n"
    variable_names[191] = "flux of Aerobic and Anaerobic Bacteria at SURFACE"
    variable_units[191] = "mmol N/m2/day"

    variable_abbrev[192] = "jsurb1p"
    variable_names[192] = "flux of Aerobic and Anaerobic Bacteria at SURFACE"
    variable_units[192] = "mmol P/m2/day"

    variable_abbrev[193] = "jsurp1c"
    variable_names[193] = "flux of Diatoms at SURFACE"
    variable_units[193] = "mg C/m2/day"

    variable_abbrev[194] = "jsurp1n"
    variable_names[194] = "flux of Diatoms at SURFACE"
    variable_units[194] = "mmol N/m2/day"

    variable_abbrev[195] = "jsurp1p"
    variable_names[195] = "flux of Diatoms at SURFACE"
    variable_units[195] = "mmol P/m2/day"

    variable_abbrev[196] = "jsurp1l"
    variable_names[196] = "flux of Diatoms at SURFACE"
    variable_units[196] = "mg Chl/m2/day"

    variable_abbrev[197] = "jsurp1s"
    variable_names[197] = "flux of Diatoms at SURFACE"
    variable_units[197] = "mmol Si/m2/day"

    variable_abbrev[198] = "jsurp2c"
    variable_names[198] = "flux of Flagellates at SURFACE"
    variable_units[198] = "mg C/m2/day"

    variable_abbrev[199] = "jsurp2n"
    variable_names[199] = "flux of Flagellates at SURFACE"
    variable_units[199] = "mmol N/m2/day"

    variable_abbrev[200] = "jsurp2p"
    variable_names[200] = "flux of Flagellates at SURFACE"
    variable_units[200] = "mmol P/m2/day"

    variable_abbrev[201] = "jsurp2l"
    variable_names[201] = "flux of Flagellates at SURFACE"
    variable_units[201] = "mg Chl/m2/day"

    variable_abbrev[202] = "jsurp3c"
    variable_names[202] = "flux of PicoPhytoplankton at SURFACE"
    variable_units[202] = "mg C/m2/day"

    variable_abbrev[203] = "jsurp3n"
    variable_names[203] = "flux of PicoPhytoplankton at SURFACE"
    variable_units[203] = "mmol N/m2/day"

    variable_abbrev[204] = "jsurp3p"
    variable_names[204] = "flux of PicoPhytoplankton at SURFACE"
    variable_units[204] = "mmol P/m2/day"

    variable_abbrev[205] = "jsurp3l"
    variable_names[205] = "flux of PicoPhytoplankton at SURFACE"
    variable_units[205] = "mg Chl/m2/day"

    variable_abbrev[206] = "jsurp4c"
    variable_names[206] = "flux of Large Phytoplankton at SURFACE"
    variable_units[206] = "mg C/m2/day"

    variable_abbrev[207] = "jsurp4n"
    variable_names[207] = "flux of Large Phytoplankton at SURFACE"
    variable_units[207] = "mmol N/m2/day"

    variable_abbrev[208] = "jsurp4p"
    variable_names[208] = "flux of Large Phytoplankton at SURFACE"
    variable_units[208] = "mmol P/m2/day"

    variable_abbrev[209] = "jsurp4l"
    variable_names[209] = "flux of Large Phytoplankton at SURFACE"
    variable_units[209] = "mg Chl/m2/day"

    variable_abbrev[210] = "jsurz3c"
    variable_names[210] = "flux of Carnivorous Mesozooplankton at SURFACE"
    variable_units[210] = "mg C/m2/day"

    variable_abbrev[211] = "jsurz3n"
    variable_names[211] = "flux of Carnivorous Mesozooplankton at SURFACE"
    variable_units[211] = "mmol N/m2/day"

    variable_abbrev[212] = "jsurz3p"
    variable_names[212] = "flux of Carnivorous Mesozooplankton at SURFACE"
    variable_units[212] = "mmol P/m2/day"

    variable_abbrev[213] = "jsurz4c"
    variable_names[213] = "flux of Omnivorous Mesozooplankton at SURFACE"
    variable_units[213] = "mg C/m2/day"

    variable_abbrev[214] = "jsurz4n"
    variable_names[214] = "flux of Omnivorous Mesozooplankton at SURFACE"
    variable_units[214] = "mmol N/m2/day"

    variable_abbrev[215] = "jsurz4p"
    variable_names[215] = "flux of Omnivorous Mesozooplankton at SURFACE"
    variable_units[215] = "mmol P/m2/day"

    variable_abbrev[216] = "jsurz5c"
    variable_names[216] = "flux of Microzooplankton at SURFACE"
    variable_units[216] = "mg C/m2/day"

    variable_abbrev[217] = "jsurz5n"
    variable_names[217] = "flux of Microzooplankton at SURFACE"
    variable_units[217] = "mmol N/m2/day"

    variable_abbrev[218] = "jsurz5p"
    variable_names[218] = "flux of Microzooplankton at SURFACE"
    variable_units[218] = "mmol P/m2/day"

    variable_abbrev[219] = "jsurz6c"
    variable_names[219] = "flux of Heterotrophic Nanoflagellates (HNAN) at SURFACE"
    variable_units[219] = "mg C/m2/day"

    variable_abbrev[220] = "jsurz6n"
    variable_names[220] = "flux of Heterotrophic Nanoflagellates (HNAN) at SURFACE"
    variable_units[220] = "mmol N/m2/day"

    variable_abbrev[221] = "jsurz6p"
    variable_names[221] = "flux of Heterotrophic Nanoflagellates (HNAN) at SURFACE"
    variable_units[221] = "mmol P/m2/day"

    variable_abbrev[222] = "jsurr1c"
    variable_names[222] = "flux of Labile Dissolved Organic Matter at SURFACE"
    variable_units[222] = "mg C/m2/day"

    variable_abbrev[223] = "jsurr1n"
    variable_names[223] = "flux of Labile Dissolved Organic Matter at SURFACE"
    variable_units[223] = "mmol N/m2/day"

    variable_abbrev[224] = "jsurr1p"
    variable_names[224] = "flux of Labile Dissolved Organic Matter at SURFACE"
    variable_units[224] = "mmol P/m2/day"

    variable_abbrev[225] = "jsurr2c"
    variable_names[225] = "flux of Semi-labile Dissolved Organic Carbon at SURFACE"
    variable_units[225] = "mg C/m2/day"

    variable_abbrev[226] = "jsurr3c"
    variable_names[226] = "flux of Semi-refractory Dissolved Organic Carbon at SURFACE"
    variable_units[226] = "mg C/m2/day"

    variable_abbrev[227] = "jsurr6c"
    variable_names[227] = "flux of Particulate Organic Matter at SURFACE"
    variable_units[227] = "mg C/m2/day"

    variable_abbrev[228] = "jsurr6n"
    variable_names[228] = "flux of Particulate Organic Matter at SURFACE"
    variable_units[228] = "mmol N/m2/day"

    variable_abbrev[229] = "jsurr6p"
    variable_names[229] = "flux of Particulate Organic Matter at SURFACE"
    variable_units[229] = "mmol P/m2/day"

    variable_abbrev[230] = "jsurr6s"
    variable_names[230] = "flux of Particulate Organic Matter at SURFACE"
    variable_units[230] = "mmol Si/m2/day"

    variable_abbrev[231] = "jsuro3c"
    variable_names[231] = "flux of Dissolved Inorganic Carbon at SURFACE"
    variable_units[231] = "mg C/m2/day"

    variable_abbrev[232] = "jsuro3h"
    variable_names[232] = "flux of Dissolved Inorganic Carbon at SURFACE"
    variable_units[232] = "mmol eq/m2/day"

    # Bottom
    variable_abbrev[233] = "jboto2o"
    variable_names[233] = "flux of Oxygen at BOTTOM"
    variable_units[233] = "mmol O2/m2/day"

    variable_abbrev[234] = "jbotn1p"
    variable_names[234] = "flux of Phosphate at BOTTOM"
    variable_units[234] = "mmol P/m2/day"

    variable_abbrev[235] = "jbotn3n"
    variable_names[235] = "flux of Nitrate at BOTTOM"
    variable_units[235] = "mmol N/m2/day"

    variable_abbrev[236] = "jbotn4n"
    variable_names[236] = "flux of Ammonium at BOTTOM"
    variable_units[236] = "mmol N/m2/day"

    variable_abbrev[237] = "jboto4n"
    variable_names[237] = "flux of o4n at BOTTOM"
    variable_units[237] = "mmol N/m2/day"

    variable_abbrev[238] = "jbotn5s"
    variable_names[238] = "flux of Silicate at BOTTOM"
    variable_units[238] = "mmol Si/m2/day"

    variable_abbrev[239] = "jbotn6r"
    variable_names[239] = "flux of Reduction Equivalents at BOTTOM"
    variable_units[239] = "mmol S--/m2/day"

    variable_abbrev[240] = "jbotb1c"
    variable_names[240] = "flux of Aerobic and Anaerobic Bacteria at BOTTOM"
    variable_units[240] = "mg C/m2/day"

    variable_abbrev[241] = "jbotb1n"
    variable_names[241] = "flux of Aerobic and Anaerobic Bacteria at BOTTOM"
    variable_units[241] = "mmol N/m2/day"

    variable_abbrev[242] = "jbotb1p"
    variable_names[242] = "flux of Aerobic and Anaerobic Bacteria at BOTTOM"
    variable_units[242] = "mmol P/m2/day"

    variable_abbrev[243] = "jbotp1c"
    variable_names[243] = "flux of Diatoms at BOTTOM"
    variable_units[243] = "mg C/m2/day"

    variable_abbrev[244] = "jbotp1n"
    variable_names[244] = "flux of Diatoms at BOTTOM"
    variable_units[244] = "mmol N/m2/day"

    variable_abbrev[245] = "jbotp1p"
    variable_names[245] = "flux of Diatoms at BOTTOM"
    variable_units[245] = "mmol P/m2/day"

    variable_abbrev[246] = "jbotp1l"
    variable_names[246] = "flux of Diatoms at BOTTOM"
    variable_units[246] = "mg Chl/m2/day"

    variable_abbrev[247] = "jbotp1s"
    variable_names[247] = "flux of Diatoms at BOTTOM"
    variable_units[247] = "mmol Si/m2/day"

    variable_abbrev[248] = "jbotp2c"
    variable_names[248] = "flux of Flagellates at BOTTOM"
    variable_units[248] = "mg C/m2/day"

    variable_abbrev[249] = "jbotp2n"
    variable_names[249] = "flux of Flagellates at BOTTOM"
    variable_units[249] = "mmol N/m2/day"

    variable_abbrev[250] = "jbotp2p"
    variable_names[250] = "flux of Flagellates at BOTTOM"
    variable_units[250] = "mmol P/m2/day"

    variable_abbrev[251] = "jbotp2l"
    variable_names[251] = "flux of Flagellates at BOTTOM"
    variable_units[251] = "mg Chl/m2/day"

    variable_abbrev[252] = "jbotp3c"
    variable_names[252] = "flux of PicoPhytoplankton at BOTTOM"
    variable_units[252] = "mg C/m2/day"

    variable_abbrev[253] = "jbotp3n"
    variable_names[253] = "flux of PicoPhytoplankton at BOTTOM"
    variable_units[253] = "mmol N/m2/day"

    variable_abbrev[254] = "jbotp3p"
    variable_names[254] = "flux of PicoPhytoplankton at BOTTOM"
    variable_units[254] = "mmol P/m2/day"

    variable_abbrev[255] = "jbotp3l"
    variable_names[255] = "flux of PicoPhytoplankton at BOTTOM"
    variable_units[255] = "mg Chl/m2/day"

    variable_abbrev[256] = "jbotp4c"
    variable_names[256] = "flux of Large Phytoplankton at BOTTOM"
    variable_units[256] = "mg C/m2/day"

    variable_abbrev[257] = "jbotp4n"
    variable_names[257] = "flux of Large Phytoplankton at BOTTOM"
    variable_units[257] = "mmol N/m2/day"

    variable_abbrev[258] = "jbotp4p"
    variable_names[258] = "flux of Large Phytoplankton at BOTTOM"
    variable_units[258] = "mmol P/m2/day"

    variable_abbrev[259] = "jbotp4l"
    variable_names[259] = "flux of Large Phytoplankton at BOTTOM"
    variable_units[259] = "mg Chl/m2/day"

    variable_abbrev[260] = "jbotz3c"
    variable_names[260] = "flux of Carnivorous Mesozooplankton at BOTTOM"
    variable_units[260] = "mg C/m2/day"

    variable_abbrev[261] = "jbotz3n"
    variable_names[261] = "flux of Carnivorous Mesozooplankton at BOTTOM"
    variable_units[261] = "mmol N/m2/day"

    variable_abbrev[262] = "jbotz3p"
    variable_names[262] = "flux of Carnivorous Mesozooplankton at BOTTOM"
    variable_units[262] = "mmol P/m2/day"

    variable_abbrev[263] = "jbotz4c"
    variable_names[263] = "flux of Omnivorous Mesozooplankton at BOTTOM"
    variable_units[263] = "mg C/m2/day"

    variable_abbrev[264] = "jbotz4n"
    variable_names[264] = "flux of Omnivorous Mesozooplankton at BOTTOM"
    variable_units[264] = "mmol N/m2/day"

    variable_abbrev[265] = "jbotz4p"
    variable_names[265] = "flux of Omnivorous Mesozooplankton at BOTTOM"
    variable_units[265] = "mmol P/m2/day"

    variable_abbrev[266] = "jbotz5c"
    variable_names[266] = "flux of Microzooplankton at BOTTOM"
    variable_units[266] = "mg C/m2/day"

    variable_abbrev[267] = "jbotz5n"
    variable_names[267] = "flux of Microzooplankton at BOTTOM"
    variable_units[267] = "mmol N/m2/day"

    variable_abbrev[268] = "jbotz5p"
    variable_names[268] = "flux of Microzooplankton at BOTTOM"
    variable_units[268] = "mmol P/m2/day"

    variable_abbrev[269] = "jbotz6c"
    variable_names[269] = "flux of Heterotrophic Nanoflagellates (HNAN) at BOTTOM"
    variable_units[269] = "mg C/m2/day"

    variable_abbrev[270] = "jbotz6n"
    variable_names[270] = "flux of Heterotrophic Nanoflagellates (HNAN) at BOTTOM"
    variable_units[270] = "mmol N/m2/day"

    variable_abbrev[271] = "jbotz6p"
    variable_names[271] = "flux of Heterotrophic Nanoflagellates (HNAN) at BOTTOM"
    variable_units[271] = "mmol P/m2/day"

    variable_abbrev[272] = "jbotr1c"
    variable_names[272] = "flux of Labile Dissolved Organic Matter at BOTTOM"
    variable_units[272] = "mg C/m2/day"

    variable_abbrev[273] = "jbotr1n"
    variable_names[273] = "flux of Labile Dissolved Organic Matter at BOTTOM"
    variable_units[273] = "mmol N/m2/day"

    variable_abbrev[274] = "jbotr1p"
    variable_names[274] = "flux of Labile Dissolved Organic Matter at BOTTOM"
    variable_units[274] = "mmol P/m2/day"

    variable_abbrev[275] = "jbotr2c"
    variable_names[275] = "flux of Semi-labile Dissolved Organic Carbon at BOTTOM"
    variable_units[275] = "mg C/m2/day"

    variable_abbrev[276] = "jbotr3c"
    variable_names[276] = "flux of Semi-refractory Dissolved Organic Carbon at BOTTOM"
    variable_units[276] = "mg C/m2/day"

    variable_abbrev[277] = "jbotr6c"
    variable_names[277] = "flux of Particulate Organic Matter at BOTTOM"
    variable_units[277] = "mg C/m2/day"

    variable_abbrev[278] = "jbotr6n"
    variable_names[278] = "flux of Particulate Organic Matter at BOTTOM"
    variable_units[278] = "mmol N/m2/day"

    variable_abbrev[279] = "jbotr6p"
    variable_names[279] = "flux of Particulate Organic Matter at BOTTOM"
    variable_units[279] = "mmol P/m2/day"

    variable_abbrev[280] = "jbotr6s"
    variable_names[280] = "flux of Particulate Organic Matter at BOTTOM"
    variable_units[280] = "mmol Si/m2/day"

    variable_abbrev[281] = "jboto3c"
    variable_names[281] = "flux of Dissolved Inorganic Carbon at BOTTOM"
    variable_units[281] = "mg C/m2/day"

    variable_abbrev[282] = "jboto3h"
    variable_names[282] = "flux of Dissolved Inorganic Carbon at BOTTOM"
    variable_units[282] = "mmol eq/m2/day"

    # River
    variable_abbrev[283] = "jrivo2o"
    variable_names[283] = "flux of Oxygen at RIVER"
    variable_units[283] = "mmol O2/m2/day"

    variable_abbrev[284] = "jrivn1p"
    variable_names[284] = "flux of Phosphate at RIVER"
    variable_units[284] = "mmol P/m2/day"

    variable_abbrev[285] = "jrivn3n"
    variable_names[285] = "flux of Nitrate at RIVER"
    variable_units[285] = "mmol N/m2/day"

    variable_abbrev[286] = "jrivn4n"
    variable_names[286] = "flux of Ammonium at RIVER"
    variable_units[286] = "mmol N/m2/day"

    variable_abbrev[287] = "jrivo4n"
    variable_names[287] = "flux of o4n at RIVER"
    variable_units[287] = "flux of o4n at RIVER"

    variable_abbrev[288] = "jrivn5s"
    variable_names[288] = "flux of Silicate at RIVER"
    variable_units[288] = "mmol Si/m2/day"

    variable_abbrev[289] = "jrivn6r"
    variable_names[289] = "flux of Reduction Equivalents at RIVER"
    variable_units[289] = "mmol S--/m2/day"

    variable_abbrev[290] = "jrivb1c"
    variable_names[290] = "flux of Aerobic and Anaerobic Bacteria at RIVER"
    variable_units[290] = "mg C/m2/day"

    variable_abbrev[291] = "jrivb1n"
    variable_names[291] = "flux of Aerobic and Anaerobic Bacteria at RIVER"
    variable_units[291] = "mmol N/m2/day"

    variable_abbrev[292] = "jrivb1p"
    variable_names[292] = "flux of Aerobic and Anaerobic Bacteria at RIVER"
    variable_units[292] = "mmol P/m2/day"

    variable_abbrev[293] = "jrivp1c"
    variable_names[293] = "flux of Diatoms at RIVER"
    variable_units[293] = "mg C/m2/day"

    variable_abbrev[294] = "jrivp1n"
    variable_names[294] = "flux of Diatoms at RIVER"
    variable_units[294] = "mmol N/m2/day"

    variable_abbrev[295] = "jrivp1p"
    variable_names[295] = "flux of Diatoms at RIVER"
    variable_units[295] = "mmol P/m2/day"

    variable_abbrev[296] = "jrivp1l"
    variable_names[296] = "flux of Diatoms at RIVER"
    variable_units[296] = "mg Chl/m2/day"

    variable_abbrev[297] = "jrivp1s"
    variable_names[297] = "flux of Diatoms at RIVER"
    variable_units[297] = "mmol Si/m2/day"

    variable_abbrev[298] = "jrivp2c"
    variable_names[298] = "flux of Flagellates at RIVER"
    variable_units[298] = "mg C/m2/day"

    variable_abbrev[299] = "jrivp2n"
    variable_names[299] = "flux of Flagellates at RIVER"
    variable_units[299] = "mmol N/m2/day"

    variable_abbrev[300] = "jrivp2p"
    variable_names[300] = "flux of Flagellates at RIVER"
    variable_units[300] = "mmol P/m2/day"

    variable_abbrev[301] = "jrivp2l"
    variable_names[301] = "flux of Flagellates at RIVER"
    variable_units[301] = "flux of Flagellates at RIVER"

    variable_abbrev[302] = "jrivp3c"
    variable_names[302] = "flux of PicoPhytoplankton at RIVER"
    variable_units[302] = "mg C/m2/day"

    variable_abbrev[303] = "jrivp3n"
    variable_names[303] = "flux of PicoPhytoplankton at RIVER"
    variable_units[303] = "mmol N/m2/day"

    variable_abbrev[304] = "jrivp3p"
    variable_names[304] = "flux of PicoPhytoplankton at RIVER"
    variable_units[304] = "mmol P/m2/day"

    variable_abbrev[305] = "jrivp3l"
    variable_names[305] = "flux of PicoPhytoplankton at RIVER"
    variable_units[305] = "mg Chl/m2/day"

    variable_abbrev[306] = "jrivp4c"
    variable_names[306] = "flux of Large Phytoplankton at RIVER"
    variable_units[306] = "mg C/m2/day"

    variable_abbrev[307] = "jrivp4n"
    variable_names[307] = "flux of Large Phytoplankton at RIVER"
    variable_units[307] = "mmol N/m2/day"

    variable_abbrev[308] = "jrivp4p"
    variable_names[308] = "flux of Large Phytoplankton at RIVER"
    variable_units[308] = "mmol P/m2/day"

    variable_abbrev[309] = "jrivp4l"
    variable_names[309] = "flux of Large Phytoplankton at RIVER"
    variable_units[309] = "mg Chl/m2/day"

    variable_abbrev[310] = "jrivz3c"
    variable_names[310] = "flux of Carnivorous Mesozooplankton at RIVER"
    variable_units[310] = "mg C/m2/day"

    variable_abbrev[311] = "jrivz3n"
    variable_names[311] = "flux of Carnivorous Mesozooplankton at RIVER"
    variable_units[311] = "mmol N/m2/day"

    variable_abbrev[312] = "jrivz3p"
    variable_names[312] = "flux of Carnivorous Mesozooplankton at RIVER"
    variable_units[312] = "mmol P/m2/day"

    variable_abbrev[313] = "jrivz4c"
    variable_names[313] = "flux of Omnivorous Mesozooplankton at RIVER"
    variable_units[313] = "mg C/m2/day"

    variable_abbrev[314] = "jrivz4n"
    variable_names[314] = "flux of Omnivorous Mesozooplankton at RIVER"
    variable_units[314] = "mmol N/m2/day"

    variable_abbrev[315] = "jrivz4p"
    variable_names[315] = "flux of Omnivorous Mesozooplankton at RIVER"
    variable_units[315] = "mmol P/m2/day"

    variable_abbrev[316] = "jrivz5c"
    variable_names[316] = "flux of Microzooplankton at RIVER"
    variable_units[316] = "mg C/m2/day"

    variable_abbrev[317] = "jrivz5n"
    variable_names[317] = "flux of Microzooplankton at RIVER"
    variable_units[317] = "mmol N/m2/day"

    variable_abbrev[318] = "jrivz5p"
    variable_names[318] = "flux of Microzooplankton at RIVER"
    variable_units[318] = "mmol P/m2/day"

    variable_abbrev[319] = "jrivz6c"
    variable_names[319] = "flux of Heterotrophic Nanoflagellates (HNAN) at RIVER"
    variable_units[319] = "mg C/m2/day"

    variable_abbrev[320] = "jrivz6n"
    variable_names[320] = "flux of Heterotrophic Nanoflagellates (HNAN) at RIVER"
    variable_units[320] = "mmol N/m2/day"

    variable_abbrev[321] = "jrivz6p"
    variable_names[321] = "flux of Heterotrophic Nanoflagellates (HNAN) at RIVER"
    variable_units[321] = "mmol P/m2/day"

    variable_abbrev[322] = "jrivr1c"
    variable_names[322] = "flux of Labile Dissolved Organic Matter at RIVER"
    variable_units[322] = "mg C/m2/day"

    variable_abbrev[323] = "jrivr1n"
    variable_names[323] = "flux of Labile Dissolved Organic Matter at RIVER"
    variable_units[323] = "mmol N/m2/day"

    variable_abbrev[324] = "jrivr1p"
    variable_names[324] = "flux of Labile Dissolved Organic Matter at RIVER"
    variable_units[324] = "mmol P/m2/day"

    variable_abbrev[325] = "jrivr2c"
    variable_names[325] = "flux of Semi-labile Dissolved Organic Carbon at RIVER"
    variable_units[325] = "mg C/m2/day"

    variable_abbrev[326] = "jrivr3c"
    variable_names[326] = "flux of Semi-refractory Dissolved Organic Carbon at RIVER"
    variable_units[326] = "mg C/m2/day"

    variable_abbrev[327] = "jrivr6c"
    variable_names[327] = "flux of Particulate Organic Matter at RIVER"
    variable_units[327] = "mg C/m2/day"

    variable_abbrev[328] = "jrivr6n"
    variable_names[328] = "flux of Particulate Organic Matter at RIVER"
    variable_units[328] = "mmol N/m2/day"

    variable_abbrev[329] = "jrivr6p"
    variable_names[329] = "flux of Particulate Organic Matter at RIVER"
    variable_units[329] = "mmol P/m2/day"

    variable_abbrev[330] = "jrivr6s"
    variable_names[330] = "flux of Particulate Organic Matter at RIVER"
    variable_units[330] = "mmol Si/m2/day"

    variable_abbrev[331] = "jrivo3c"
    variable_names[331] = "flux of Dissolved Inorganic Carbon at RIVER"
    variable_units[331] = "mg C/m2/day"

    variable_abbrev[332] = "jrivo3h"
    variable_names[332] = "flux of Dissolved Inorganic Carbon at RIVER"
    variable_units[332] = "mmol eq/m2/day"


    # Indexes
    pelagic_index = PelagicIndex(0,332,0,49,50,140,141,170,171,182,183,232,233,282,283,331)
    index = Index(0,332)

    # try:
    #     INCLUDE_SEAICE
    # except NameError:
    #     INCLUDE_SEAICE = False
    # else:
    #     INCLUDE_SEAICE = True
    if INCLUDE_SEAICE:
        seaice_index = SeaIceIndex(333,332,333,332,333,332,333,332)
    else:
        seaice_index = None

    # try:
    #     INCLUDE_BEN
    # except NameError:
    #     INCLUDE_BEN = False
    # else:
    #     INCLUDE_BEN = True
    if INCLUDE_BEN:
        benthic_index = BenthicIndex(333,332,333,332,333,332,333,332)
    else:
        benthic_index = None


    return variable_abbrev, variable_names, variable_units, index, pelagic_index, seaice_index, benthic_index



























