import numpy as np
from pom.initialize_variables import read_pom_input
from pom.modules import vertical_layers, water_specific_heat_times_density

wind_speed_zonal, wind_speed_meridional, surface_salinity, solar_radiation, inorganic_suspended_matter, \
    salinity_climatology, temperature_climatology, w_velocity_climatology, w_eddy_velocity_1, \
    w_eddy_velocity_2, salinity, temperature, shortwave_radiation, surface_heat_flux, kinetic_energy_loss, \
    NO3_s, NH4_s, PO4_s, SIO4_s, O2_b, NO3_b, PO4_b, PON_b                          = read_pom_input()

# Initialize arrays
sclim = np.zeros(vertical_layers)          # salinity_climatology[:][month_counter]
tclim = np.zeros(vertical_layers)          # temperature_climatology[:][month_counter]
wclim = np.zeros(vertical_layers)          # w_velocity_climatology[:][month_counter]
weddyi = np.zeros(vertical_layers)         # w_eddy_velocity_1[:][month_counter]
weddyj = np.zeros(vertical_layers)         # w_eddy_velocity_2[:][month_counter]
ism   = np.zeros(vertical_layers-1)        # inorganic_suspended_matter[0:vertical_layers-1][month_counter]

def write_forcing_data(month_counter):

    for i in range(0,vertical_layers):
        sclim[i] = salinity_climatology[i][month_counter]
        tclim[i] = temperature_climatology[i][month_counter]
        wclim[i] = w_velocity_climatology[i][month_counter]
        weddyi[i] = w_eddy_velocity_1[i][month_counter]
        weddyj[i] = w_eddy_velocity_2[i][month_counter]
        if i < vertical_layers:
            ism[i] = inorganic_suspended_matter[i][month_counter]

    wsu = wind_speed_zonal[month_counter]
    wsv = wind_speed_meridional[month_counter]
    swrad  = shortwave_radiation[month_counter]
    wtsurf = surface_heat_flux[month_counter]
    qcorr  = kinetic_energy_loss[month_counter]

    NO3_si  = NO3_s[month_counter]
    NH4_si  = NH4_s[month_counter]
    PO4_si  = PO4_s[month_counter]
    SIO4_si = SIO4_s[month_counter]

    O2_bi  = O2_b[month_counter]
    NO3_bi = NO3_b[month_counter]
    PO4_bi = PO4_b[month_counter]
    PON_bi = PON_b[month_counter]

    # Convert to POM units
    wsu = -wsu * 0.001                                      # N/m2-->m2/s2
    wsv = -wsv * 0.001                                      # N/m2-->m2/s2
    swrad  = -swrad / water_specific_heat_times_density     # W/m2-->deg.C*m/s
    wtsurf = -wtsurf / water_specific_heat_times_density    # W/m2-->deg.C*m/s

    return sclim, tclim, wclim, weddyi, weddyj, ism, wsu, wsv, swrad, wtsurf, qcorr, \
           NO3_si, NH4_si, PO4_si, SIO4_si, O2_bi, NO3_bi, PO4_bi, PON_bi
