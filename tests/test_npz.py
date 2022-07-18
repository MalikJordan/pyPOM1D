import numpy as np
from pom.constants import vertical_layers, seconds_per_day, twice_the_timestep
from inputs import params_POMBFM
import matplotlib.pyplot as plt
from pom.create_profiles import calculate_vertical_temperature_and_salinity_profiles
from pom_bfm_coupling.calculations import calculate_vertical_advection
from pom_bfm_coupling.calculate_average_field import calculate_daily_average_field
from pom_bfm_coupling.data_classes import BfmStateVariableData

max_growth_rate = 1. / seconds_per_day              # s^-1
nitrogen_half_saturation = 1.                       # μmol N l^-1
max_grazing_rate = 1. / seconds_per_day             # s^-1
zooplankton_death_rate = 0.2 / seconds_per_day      # s^-1
grazing = 0.2                                       # μmol N l^-1
phytoplankton_death_rate = 0.1 / seconds_per_day    # s^-1
zooplankton_assimiilated_nitrogen = 0.7
light_intensity = 0.25

iterations_needed = params_POMBFM.idays * seconds_per_day / params_POMBFM.dti  
num_boxes = vertical_layers - 1

def NPZrates(NPZ):
    num_boxes = vertical_layers - 1
    npz_rates = np.zeros((num_boxes,3))
    for j in range(0,num_boxes):
        # phyto_coeff = max_growth_rate * ( NPZ[2,j,i] / (nitrogen_half_saturation + NPZ[2,j,i]) ) * light_intensity
        # zoo_coeff1 = max_grazing_rate * ( 1 - np.exp(-grazing*NPZ[0,j,i]) )
        # zoo_coeff2 = zooplankton_assimiilated_nitrogen * max_grazing_rate * ( 1 - np.exp(-grazing*NPZ[0,j,i]) )
        # zoo_coeff3 = ( 1 - zooplankton_assimiilated_nitrogen ) * max_grazing_rate * ( 1 - np.exp(-grazing*NPZ[0,j,i]) )
        
        # NPZ[0,j,i+1] = NPZ[0,j,i] + params_POMBFM.dti * ( phyto_coeff * NPZ[0,j,i] - zoo_coeff1 * NPZ[1,j,i] - phytoplankton_death_rate * NPZ[0,j,i] )
        # NPZ[1,j,i+1] = NPZ[1,j,i] + params_POMBFM.dti * ( zoo_coeff2 * NPZ[1,j,i] - zooplankton_death_rate * NPZ[1,j,i] )
        # NPZ[2,j,i+1] = NPZ[2,j,i] + params_POMBFM.dti * ( -phyto_coeff * NPZ[0,j,i] + zoo_coeff3 * NPZ[1,j,i] + phytoplankton_death_rate * NPZ[0,j,i] + zooplankton_death_rate * NPZ[1,j,i] )
        
        # Updates for 4/6
        phyto_coeff = max_growth_rate * ( NPZ[j,2] / (nitrogen_half_saturation + NPZ[j,2]) ) * light_intensity
        zoo_coeff1 = max_grazing_rate * ( 1 - np.exp(-grazing*NPZ[j,0]) )
        zoo_coeff2 = zooplankton_assimiilated_nitrogen * max_grazing_rate * ( 1 - np.exp(-grazing*NPZ[j,0]) )
        zoo_coeff3 = ( 1 - zooplankton_assimiilated_nitrogen ) * max_grazing_rate * ( 1 - np.exp(-grazing*NPZ[j,0]) )
        
        npz_rates[j,0] = params_POMBFM.dti * ( phyto_coeff * NPZ[j,0] - zoo_coeff1 * NPZ[j,1] - phytoplankton_death_rate * NPZ[j,0] )
        npz_rates[j,1] = params_POMBFM.dti * ( zoo_coeff2 * NPZ[j,1] - zooplankton_death_rate * NPZ[j,1] )
        npz_rates[j,2] = params_POMBFM.dti * ( -phyto_coeff * NPZ[j,0] + zoo_coeff3 * NPZ[j,1] + phytoplankton_death_rate * NPZ[j,0] + zooplankton_death_rate * NPZ[j,1] )
    
    return npz_rates


# def NPZrates(NPZ, NPZb):
#     num_boxes = vertical_layers - 1
#     for j in range(0,num_boxes):

#         # Updates for 4/6
#         phyto_coeff = max_growth_rate * ( NPZb[j,2] / (nitrogen_half_saturation + NPZb[j,2]) ) * light_intensity
#         zoo_coeff1 = max_grazing_rate * ( 1 - np.exp(-grazing*NPZb[j,0]) )
#         zoo_coeff2 = zooplankton_assimiilated_nitrogen * max_grazing_rate * ( 1 - np.exp(-grazing*NPZb[j,0]) )
#         zoo_coeff3 = ( 1 - zooplankton_assimiilated_nitrogen ) * max_grazing_rate * ( 1 - np.exp(-grazing*NPZb[j,0]) )
        
#         NPZ[j,0] = NPZb[j,0] + params_POMBFM.dti * ( phyto_coeff * NPZb[j,0] - zoo_coeff1 * NPZb[j,0] - phytoplankton_death_rate * NPZb[j,0] )
#         NPZ[j,1] = NPZb[j,1] + params_POMBFM.dti * ( zoo_coeff2 * NPZb[j,1] - zooplankton_death_rate * NPZb[j,1] )
#         NPZ[j,2] = NPZb[j,2] + params_POMBFM.dti * ( -phyto_coeff * NPZb[j,0] + zoo_coeff3 * NPZb[j,1] + phytoplankton_death_rate * NPZb[j,0] + zooplankton_death_rate * NPZb[j,1] )

#     return NPZ


def calculate_vertical_diffusivity_NPZ(vertical_grid, diffusion, NPZ, NPZb, npz_rates, bfm_phys_vars):
   
    npz_state_var = BfmStateVariableData()
    sinking_velocity = np.zeros(vertical_layers)

    # The input general cir. vertical vel. is suppose to be in m/s
    W_ON = 1.0
    # The input eddy vertical vel. is provided in m/d
    Weddy_ON = 0.1/86400.0  # to m/s

    for M in range(0,3):
        # ZEROING
        npz_state_var.surface_flux = 0.
        npz_state_var.bottom_flux = 0.
        npz_state_var.current[:] = 0.
        npz_state_var.backward[:] = 0.
        npz_state_var.forward[:] = 0.
        sinking_velocity[:] = 0.

        # LOAD NPZ STATE VAR.
        for i in range(0,vertical_layers-1):
            npz_state_var.current[i] = NPZ[i,M]
            npz_state_var.backward[i] = NPZb[i,M]

        npz_state_var.current[vertical_layers-1] = npz_state_var.current[vertical_layers-2]
        npz_state_var.backward[vertical_layers-1] = npz_state_var.backward[vertical_layers-2]

        for i in range(0,vertical_layers):
            sinking_velocity[i] = W_ON*bfm_phys_vars.wgen[i] + Weddy_ON*bfm_phys_vars.weddy[i]

        # if M == 0:
        #     phyto_sedimentation_rates = phyto_sedimentation()
        #     for i in range(0,vertical_layers-1):
        #         sinking_velocity[i] = sinking_velocity[i] - phyto_sedimentation_rates[i,0]/seconds_per_day
        #     # FINAL SINK VALUE
        #     sinking_velocity[vertical_layers - 1] = sinking_velocity[vertical_layers - 2]

        # SINKING: UPSTREAM VERTICAL ADVECTION
        npz_state_var = calculate_vertical_advection(npz_state_var,sinking_velocity,vertical_grid)

        # SOURCE SPLITTING LEAPFROG INTEGRATION
        for i in range(0,vertical_layers-1):
            # npz_state_var.forward[i] = npz_state_var.backward[i] + twice_the_timestep*((npz_state_var.forward[i]/params_POMBFM.h) + npz_rates[i,M])
            npz_state_var.forward[i] = npz_state_var.current[i] + npz_state_var.forward[i] + npz_rates[i,M]

        # COMPUTE VERTICAL DIFFUSION
        npz_state_var = calculate_vertical_temperature_and_salinity_profiles(vertical_grid, diffusion, npz_state_var, 0, params_POMBFM.nbcbfm, params_POMBFM.umolbfm)

        for N in range(0,vertical_layers-1):
            NPZb[N,M] = npz_state_var.current[N] + 0.5*params_POMBFM.smoth*(npz_state_var.forward[N] + npz_state_var.backward[N] - 2.*npz_state_var.current[N])
            NPZ[N,M] = npz_state_var.forward[N]

    return NPZ, NPZb


def pom_npz_1d(vertical_grid, diffusion, bfm_phys_vars, NPZ, NPZb, NPZave):
    
    num_boxes = vertical_layers - 1
    npz_rates = np.zeros((num_boxes,3))

    npz_rates = NPZrates(NPZ)
    calculate_vertical_diffusivity_NPZ(vertical_grid, diffusion, NPZ, NPZb, npz_rates, bfm_phys_vars)

    NPZave.count = NPZave.count + 1
    if (NPZave.count > 0) and (NPZave.count < seconds_per_day/params_POMBFM.dti):
        NPZave = calculate_daily_average_field(NPZave,NPZ,'Accumulate')
    elif NPZave.count == seconds_per_day/params_POMBFM.dti:
        NPZave = calculate_daily_average_field(NPZave,NPZ,'Mean')
        NPZave.count = 0

    return NPZ, NPZb, NPZave



# def phyto_sedimentation():

    # FROM namelists_bfm
    # p_rPIm        [m/d]   phytoplanktom background sinking rate
    # p_burvel_PI   [m/d]   Botttom burial velocity for detritus
    p_rPIm = [0.0, 0.0, 0.0, 0.0]
    p_burvel_PI = 0.0

    # FROM MODULEMEM.F90 (338)
    iiPhytoPlankton = 4
    iiP1 = 1; iiP2 = 2; iiP3 = 3; iiP4 = 4

    # FROM PelGLobal.F90 (149-154)
    phyto_sedimentation_rates = np.zeros((vertical_layers-1,iiPhytoPlankton))
    for i in range(0,iiPhytoPlankton):
        phyto_sedimentation_rates[:,i] = p_rPIm[i]
        # try:
        #     BFM_POM
        # except NameError:
        #     BFM_POM = False
        if not BFM_POM:
            phyto_sedimentation_rates[vertical_layers-2,i] = p_burvel_PI

    return phyto_sedimentation_rates





def NPZplots(NPZave):
    # plt.subplots()
    # phyto_data = NPZ[0,:,:]
    # zoo_data = NPZ[1,:,:]
    # nut_data = NPZ[2,:,:]

    phyto_data = NPZave.daily_ave[:,0,:]
    zoo_data = NPZave.daily_ave[:,1,:]
    nut_data = NPZave.daily_ave[:,2,:]
    
    fig1 = plt.figure()
    phyto = plt.imshow(phyto_data)
    fig1.colorbar(phyto)
    plt.xlabel('Time (Days)')
    plt.ylabel('Depth (m)')
    plt.title('Phytoplankton (μmol N l^-1)')

    fig2 = plt.figure()
    zoo = plt.imshow(zoo_data)
    fig2.colorbar(zoo)
    plt.xlabel('Time (Days)')
    plt.ylabel('Depth (m)')
    plt.title('Zooplankton (μmol N l^-1)')

    fig3 = plt.figure()
    nut = plt.imshow(nut_data)
    fig3.colorbar(nut)
    plt.xlabel('Time (Days)')
    plt.ylabel('Depth (m)')
    plt.title('Nutrient (μmol N l^-1)')

    plt.show()



class AverageDataNPZ:
    def __init__(self,count=0,day=0,single_day_ave=np.zeros((num_boxes,3)),daily_ave=np.zeros((num_boxes,3,params_POMBFM.idays))):
        self.count = count
        self.day = day
        self.single_day_ave = single_day_ave
        self.daily_ave = daily_ave



def checkNPZrates(NPZcheck, j):
# Updates for 4/6
    phyto_coeff = max_growth_rate * ( NPZcheck[j,2] / (nitrogen_half_saturation + NPZcheck[j,2]) ) * light_intensity
    zoo_coeff1 = max_grazing_rate * ( 1 - np.exp(-grazing*NPZcheck[j,0]) )
    zoo_coeff2 = zooplankton_assimiilated_nitrogen * max_grazing_rate * ( 1 - np.exp(-grazing*NPZcheck[j,0]) )
    zoo_coeff3 = ( 1 - zooplankton_assimiilated_nitrogen ) * max_grazing_rate * ( 1 - np.exp(-grazing*NPZcheck[j,0]) )
        
    NPZcheck[j+1,0] = NPZcheck[j,0] + params_POMBFM.dti * ( phyto_coeff * NPZcheck[j,0] - zoo_coeff1 * NPZcheck[j,1] - phytoplankton_death_rate * NPZcheck[j,0] )
    NPZcheck[j+1,1] = NPZcheck[j,1] + params_POMBFM.dti * ( zoo_coeff2 * NPZcheck[j,1] - zooplankton_death_rate * NPZcheck[j,1] )
    NPZcheck[j+1,2] = NPZcheck[j,2] + params_POMBFM.dti * ( -phyto_coeff * NPZcheck[j,0] + zoo_coeff3 * NPZcheck[j,1] + phytoplankton_death_rate * NPZcheck[j,0] + zooplankton_death_rate * NPZcheck[j,1] )

    return NPZcheck