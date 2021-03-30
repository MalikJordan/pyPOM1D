from cppdefs import *
import numpy as np
from main_pombfm1d import current_path

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# BFM_General.nml
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# bfm_nml
bio_calc = True
bfm_init = 0
bfm_rstctl = False
bio_setup = 1
out_fname = 'bfm17_pom1d'
out_dir = '.'
out_title = 'bfm17_pom1d'
out_delta = 1
parallel_log = False

# param_parameters
calcpelagicflag = True
calcbenthicflag = 0
calcconservationflag = False
calctransportflag = False
calcphytoplankton = [False,True,False,False]
calcpelbacteria = False
calcmicrozooplankton = [True,False]
calcmesozooplankton = [False,False]
calcpelchemistry = True
assignpelbenfluxesinbfmflag = False
assignairpelfluxesinbfmflag = True
chldynamicsflag = 2
lightperiodflag = 1
lightlocationflag = 3
check_fixed_quota = 0
p_small = 1e-20
slp0 = 1013.25
p_par = 0.4
p_eps0 = 0.0435
p_epsr6 = 0.0001
p_epsess = 0
p_pe_r1c = 0.6
p_pe_r1n = 0.72
p_pe_r1p = 0.832
p_qro = 0.5
p_qon_dentri = 1.25
p_qon_nitri = 2.0
p_poro0 = 0.75
p_d_tot = 0.3

# param_parameters_ben
calcbenorganisms = [False,False,False,False,False]
calcbenbacteria = [False,False]
p_initsink = 100.0
p_d_tot_2 = 0.35
p_cld1d2m = 0.01
p_cldxm = 0.001
p_q10diff = 1.49
calc_init_bennut_states = 0
p_qnqic = 0.11155
p_qpqic = 0.010255
p_qsqic = 0.0221

# bfm_init_nml
o2o0 = 219.0
n1p0 = 0.003
n3n0 = 0.04
n4n0 = 0.008
n5s0 = 0.0
n6r0 = 0.0
o3c0 = 0.0
o3h0 = 0.0
o4n0 = 0.0
p1c0 = 0.0
p2c0 = 11.5
p3c0 = 0.0
p4c0 = 0.0
z3c0 = 0.0
z4c0 = 0.0
z5c0 = 11.5
z6c0 = 0.0
b1c0 = 0.0
r1c0 = 12.4
r2c0 = 0.0
r3c0 = 0.0
r6c0 = 12.4

# bfm_ic_nml
phyto_input = '/inputs/BFM17_BERM_INIT/init_prof_Pc_150m_bermuda_killworth.da'
zoop_input = '/inputs/BFM17_BERM_INIT/init_prof_Zc_150m_bermuda_killworth.da'
poc_input = '/inputs/BFM17_BERM_INIT/init_prof_POC_150m_bermuda_killworth.da'
doc_input = '/inputs/BFM17_BERM_INIT/init_prof_DOC_150m_bermuda_killworth.da'
phos_input = '/inputs/BFM17_BERM_INIT/init_prof_P_150m_bermuda_killworth.da'
nit_input = '/inputs/BFM17_BERM_INIT/init_prof_N_150m_bermuda_killworth.da'
am_input = '/inputs/BFM17_BERM_INIT/init_prof_Am_150m_bermuda_killworth.da'
oxy_input = '/inputs/BFM17_BERM_INIT/init_prof_Oxy_150m_bermuda_killworth.da'

# bfm_init_nml_ben
y1c0 = 30.0
y2c0 = 610.0
y3c0 = 140.0
y4c0 = 220.0
y5c0 = 160.0
h1c0 = 120.0
h2c0 = 45.0
k1p0 = 0.1
k11p0 = 80.0
k21p0 = 19.11096
k4n0 = 0.0277856
k14n0 = 1.0838474
k24n0 = 100.0
k3n0 = 0.0252624
k5s0 = 8.463525
k6r0 = 100.0
d1m0 = 0.002
d2m0 = 0.025
d6m0 = 0.25
d7m0 = 0.25
d8m0 = 0.25
d9m0 = 0.25
q6c0 = 10250.0
q6n0 = 120.0
q6p0 = 10.0
q6s0 = 88.2
q1c0 = 10.4988
q11c0 = 10.4988
g2o0 = 0.67
g3c0 = 120.15
g13c0 = 440.8475
g23c0 = 11920.0
g3h0 = 10.35
g13h0 = 50.0
g23h0 = 1192.0


# bfm_save_nml
#     ave_save = 'ETW', 'O2o', 'DIC', 'EIR', 'ESW', 'ERHO', 'xEPS', 'Chla',
#                'N1p', 'N3n', 'N4n', 'P2c', 'P2n', 'P2p', 'P2l', 'Z5c', 'Z5n',
#                'Z5p', 'R1c', 'R1n', 'R1p', 'R6c', 'R6n', 'R6p', 'eiPPY(iiP1)',
#                'eiPPY(iiP2)', 'eiPPY(iiP3)', 'eiPPY(iiP4)', 'sunPPY(iiP1)',
#                'sunPPY(iiP2)', 'sunPPY(iiP3)', 'sunPPY(iiP4)', 'ruPTc',
#                'resPP', 'resZT', 'ruPTn', 'ruPTp', 'exPP', 'ruZTc', 'netZTc',
#                'rePTp', 'reBn', 'reBp', 'ruBn', 'ruBp', 'EPR'
