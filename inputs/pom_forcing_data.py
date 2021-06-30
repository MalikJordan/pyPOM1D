import numpy as np
from pom.initialize_variables import read_pom_input
from pom.modules import vertical_layers, water_specific_heat_times_density
from pom.data_classes import MonthlyForcingData

wind_speed_zonal, wind_speed_meridional, surface_salinity, solar_radiation, inorganic_suspended_matter, \
    salinity_climatology, temperature_climatology, w_velocity_climatology, w_eddy_velocity_1, \
    w_eddy_velocity_2, salinity, temperature, shortwave_radiation, surface_heat_flux, kinetic_energy_loss, \
    NO3_surf, NH4_surf, PO4_surf, SIO4_surf, O2_bott, NO3_bott, PO4_bott, PON_bott                  = read_pom_input()

def write_forcing_data(month_counter):

    sclim = salinity_climatology[:,month_counter]
    tclim = temperature_climatology[:,month_counter]
    wclim = w_velocity_climatology[:,month_counter]
    weddy1 = w_eddy_velocity_1[:,month_counter]
    weddy2 = w_eddy_velocity_2[:,month_counter]
    ism = inorganic_suspended_matter[0:vertical_layers-1,month_counter]

    wsu = wind_speed_zonal[month_counter]
    wsv = wind_speed_meridional[month_counter]
    swrad  = shortwave_radiation[month_counter]
    wtsurf = surface_heat_flux[month_counter]
    qcorr  = kinetic_energy_loss[month_counter]

    NO3_s  = NO3_surf[month_counter]
    NH4_s  = NH4_surf[month_counter]
    PO4_s  = PO4_surf[month_counter]
    SIO4_s = SIO4_surf[month_counter]

    O2_b  = O2_bott[month_counter]
    NO3_b = NO3_bott[month_counter]
    PO4_b = PO4_bott[month_counter]
    PON_b = PON_bott[month_counter]

    # Convert to POM units
    wsu = -wsu * 0.001                                      # N/m2-->m2/s2
    wsv = -wsv * 0.001                                      # N/m2-->m2/s2
    swrad  = -swrad / water_specific_heat_times_density     # W/m2-->deg.C*m/s
    wtsurf = -wtsurf / water_specific_heat_times_density    # W/m2-->deg.C*m/s

    month_data = MonthlyForcingData(sclim, tclim, wclim, weddy1, weddy2, ism, wsu, wsv, swrad, wtsurf, qcorr,
                                    NO3_s, NH4_s, PO4_s, SIO4_s, O2_b, NO3_b, PO4_b, PON_b)

    return month_data

