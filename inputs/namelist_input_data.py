from cppdefs import *
import numpy as np

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
p_pe_labileDOM_NO_C = 0.6
p_pe_labileDOM_NO_N = 0.72
p_pe_labileDOM_NO_P = 0.832
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
disOxygen_IO_O0 = 219.0
phospate_IO_P0 = 0.003
nitrate_IO_N0 = 0.04
ammonium_IO_N0 = 0.008
silicate_IO_Si0 = 0.0
reductEquiv_IO_R0 = 0.0
disInorgCarbon_IO_C0 = 0.0
o3h0 = 0.0
o4n0 = 0.0
diatoms_LO_C0 = 0.0
nanoflagellates_LO_C0 = 11.5
picophyto_LO_C0 = 0.0
largephyto_LO_C0 = 0.0
carnivMesozoo_LO_C0 = 0.0
omnivMesozoo_LO_C0 = 0.0
microzoo_LO_C0 = 11.5
z6c0 = 0.0
pelBacteria_LO_C0 = 0.0
labileDOM_NO_C0 = 12.4
semilabileDOC_NO_C0 = 0.0
semirefractDOC_NO_C0 = 0.0
particOrganDetritus_NO_C0 = 12.4

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
#     ave_save = 'ETW', 'disOxygen_IO_O', 'DIC', 'EIR', 'ESW', 'ERHO', 'xEPS', 'Chla',
#                'phospate_IO_P', 'nitrate_IO_N', 'ammonium_IO_N', 'nanoflagellates_LO_C', 'nanoflagellates_LO_N', 'nanoflagellates_LO_P', 'nanoflagellates_LO_Chl', 'microzoo_LO_C', 'microzoo_LO_N',
#                'microzoo_LO_P', 'labileDOM_NO_C', 'labileDOM_NO_N', 'labileDOM_NO_P', 'particOrganDetritus_NO_C', 'particOrganDetritus_NO_N', 'particOrganDetritus_NO_P', 'eiPPY(iiP1)',
#                'eiPPY(iiP2)', 'eiPPY(iiP3)', 'eiPPY(iiP4)', 'sunPPY(iiP1)',
#                'sunPPY(iiP2)', 'sunPPY(iiP3)', 'sunPPY(iiP4)', 'ruPTc',
#                'resPP', 'resZT', 'ruPTn', 'ruPTp', 'exPP', 'ruZTc', 'netZTc',
#                'rePTp', 'reBn', 'reBp', 'ruBn', 'ruBp', 'EPR'


# params_pombfm
h = 150.0           # bottom_depth
dti = 100.0         # time_step
alat = 45.0         # latitude
idiagn = 1          # prog_diag_switch
idays = 3600        # length_of_run
smoth = 0.1         # asselin_parameter
ihotst = 0          # hot_cold_switch
kl1 = 2             # surf_layers_log_dist
kl2 = 150           # bot_layers_log_dist
savef = 1
nrt_disOxygen_IO_O = 0.06      # relax_vel_oxygen
nrt_phospate_IO_P = 0.06      # relax_vel_phosphate
nrt_nitrate_IO_N = 0.06      # relax_vel_nitrate
nrt_ammonium_IO_N = 0.05      # relax_vel_ammonium
nbct = 2            # temp_bc_flag
nbcs = 1            # sal_bc_flag
nbcbfm = 1          # bfm_bc_flag
umol = 1e-06        # background_diffusion
umolt = 1e-07       # background_diffusion_temp
umols = 1.3e-07     # background_diffusion_sal
umolbfm = 0.0001    # background_diffusion_bfm
ntp = 2             # jerlov_flag
trt = 0             # relax_time_temp
srt = 1             # relax_time_sal
upperh = 5.0        #
ssrt = 5.68         # surf_sal_relaxation_vel






