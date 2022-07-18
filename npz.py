from cProfile import label
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
from pom.constants import vertical_layers, seconds_per_day
from inputs import params_POMBFM

max_growth_rate = 1. / seconds_per_day              # s^-1
nitrogen_half_saturation = 1.                       # μmol N l^-1
max_grazing_rate = 1. / seconds_per_day             # s^-1
zooplankton_death_rate = 0.2 / seconds_per_day      # s^-1
grazing = 0.2                                       # μmol N l^-1
phytoplankton_death_rate = 0.1 / seconds_per_day    # s^-1
zooplankton_assimiilated_nitrogen = 0.7
light_intensity = 0.25

iterations_needed = params_POMBFM.idays * seconds_per_day / params_POMBFM.dti  
t_span = np.linspace(0,int(iterations_needed),int(iterations_needed)+1)

def dYdt(Y,t):
    p, z, n = Y
    dpdt = ( max_growth_rate * (n/ (nitrogen_half_saturation + n)) * light_intensity * p ) \
                - ( max_grazing_rate * (1 - np.exp(-grazing*p)) * z ) - ( phytoplankton_death_rate * p )
    dzdt = ( zooplankton_assimiilated_nitrogen * max_grazing_rate * (1 - np.exp(-grazing*p)) * z ) - ( zooplankton_death_rate * z )
    dndt = -( max_growth_rate * (n/ (nitrogen_half_saturation + n)) * light_intensity * p ) \
                - ( (1-zooplankton_assimiilated_nitrogen) * max_grazing_rate * (1 - np.exp(-grazing*p)) * z ) \
                    + ( phytoplankton_death_rate * p ) + ( zooplankton_death_rate * z )

    return [dpdt, dzdt,  dndt]

# t_span = np.linspace(0,int(params_POMBFM.idays)*seconds_per_day)
p0 = 2.5
z0 = 0.5
n0 = 4.
Y0 = [p0, z0, n0]

sol = odeint(dYdt, Y0, t_span)
# sol2 = solve_ivp(dYdt, t_span, Y0)

fig1, ax = plt.subplots()
# fig1 = plt.figure()
# ax = fig1.add_subplot(1,3,1)
phyto, = ax.plot(t_span,sol[:,0], label='Phytoplankton')
zoo, = ax.plot(t_span,sol[:,1], label='Zooplankton')
nut, = ax.plot(t_span,sol[:,2], label='Nutrients')
ax.legend(loc='upper right')
plt.title('odeint')
plt.xlabel('Time (seconds)')
plt.ylabel('Concentration (μmol Ν l^-1)')
# plt.show()

# fig1, ax = plt.subplots()
# ax = fig1.add_subplot(1,3,2)
# phyto, = ax.plot(t_span,sol2[0,:], label='Phytoplankton')
# zoo, = ax.plot(t_span,sol2[1,:], label='Zooplankton')
# nut, = ax.plot(t_span,sol2[2,:], label='Nutrients')
# ax.legend(loc='upper right')
# plt.title('solve_ivp')
# plt.xlabel('Time (seconds)')
# plt.ylabel('Concentration (μmol Ν l^-1)')
# plt.show()

# ---------- 0D ----------
NPZ = np.zeros((3,int(iterations_needed)+1))
NPZ[0,0] = 2.5    # μmol N l^-1
NPZ[1,0] = 0.5    # μmol N l^-1
NPZ[2,0] = 4.     # μmol N l^-1

for i in range(0,int(iterations_needed)):
    phyto_coeff = max_growth_rate * ( NPZ[2,i] / (nitrogen_half_saturation + NPZ[2,i]) ) * light_intensity
    zoo_coeff1 = max_grazing_rate * ( 1 - np.exp(-grazing*NPZ[0,i]) )
    zoo_coeff2 = zooplankton_assimiilated_nitrogen * max_grazing_rate * ( 1 - np.exp(-grazing*NPZ[0,i]) )
    zoo_coeff3 = ( 1 - zooplankton_assimiilated_nitrogen ) * max_grazing_rate * ( 1 - np.exp(-grazing*NPZ[0,i]) )
    
    NPZ[0,i+1] = NPZ[0,i] + params_POMBFM.dti * ( phyto_coeff * NPZ[0,i] - zoo_coeff1 * NPZ[1,i] - phytoplankton_death_rate * NPZ[0,i] )
    NPZ[1,i+1] = NPZ[1,i] + params_POMBFM.dti * ( zoo_coeff2 * NPZ[1,i] - zooplankton_death_rate * NPZ[1,i] )
    NPZ[2,i+1] = NPZ[2,i] + params_POMBFM.dti * ( -phyto_coeff * NPZ[0,i] + zoo_coeff3 * NPZ[1,i] + phytoplankton_death_rate * NPZ[0,i] + zooplankton_death_rate * NPZ[1,i] )
    
fig2, ax = plt.subplots()
# ax = fig1.add_subplot(1,3,3)
phyto, = ax.plot(t_span,NPZ[0,:], label='Phytoplankton')
zoo, = ax.plot(t_span,NPZ[1,:], label='Zooplankton')
nut, = ax.plot(t_span,NPZ[2,:], label='Nutrients')
ax.legend(loc='upper right')
plt.title('forward euler')
plt.xlabel('Time (days)')
plt.ylabel('Concentration (μmol Ν l^-1)')
plt.show()



# ---------- 1D ----------
# NPZ = np.zeros((3,int(vertical_layers),int(iterations_needed)))
# NPZ[0,:,0] = 2.5    # μmol N l^-1
# NPZ[1,:,0] = 0.5    # μmol N l^-1
# NPZ[2,:,0] = 4.     # μmol N l^-1

# for i in range(0,int(iterations_needed)-1):
#     for j in range(0,int(vertical_layers)):
#         phyto_coeff = max_growth_rate * ( NPZ[2,j,i] / (nitrogen_half_saturation + NPZ[2,j,i]) ) * light_intensity
#         zoo_coeff1 = max_grazing_rate * ( 1 - np.exp(-grazing*NPZ[0,j,i]) )
#         zoo_coeff2 = zooplankton_assimiilated_nitrogen * max_grazing_rate * ( 1 - np.exp(-grazing*NPZ[0,j,i]) )
#         zoo_coeff3 = ( 1 - zooplankton_assimiilated_nitrogen ) * max_grazing_rate * ( 1 - np.exp(-grazing*NPZ[0,j,i]) )
        
#         NPZ[0,j,i+1] = NPZ[0,j,i] + params_POMBFM.dti * ( phyto_coeff * NPZ[0,j,i] - zoo_coeff1 * NPZ[1,j,i] - phytoplankton_death_rate * NPZ[0,j,i] )
#         NPZ[1,j,i+1] = NPZ[1,j,i] + params_POMBFM.dti * ( zoo_coeff2 * NPZ[1,j,i] - zooplankton_death_rate * NPZ[1,j,i] )
#         NPZ[2,j,i+1] = NPZ[2,j,i] + params_POMBFM.dti * ( -phyto_coeff * NPZ[0,j,i] + zoo_coeff3 * NPZ[1,j,i] + phytoplankton_death_rate * NPZ[0,j,i] + zooplankton_death_rate * NPZ[1,j,i] )

print('done')
# xx, yy = np.meshgrid(np.linspace(0,int(params_POMBFM.idays)), np.linspace(0,int(vertical_layers)))
# fig, ax = plt.subplots()
# phyto, = ax.imshow(t,NPZ[0,:,:], label='Phytoplankton')
# zoo, = ax.imshow(t,NPZ[1,:,:], label='Zooplankton')
# nut, = ax.imshow(t,NPZ[2,:,:], label='Nutrients')
# ax.legend()
# plt.xlabel('Time (days)')
# plt.ylabel('Concentration (μmol Ν l^-1)')
# plt.colorbar()
# plt.show()