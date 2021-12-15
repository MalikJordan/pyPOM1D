import numpy as np
from pom.constants import vertical_layers, water_specific_heat_times_density
from inputs import params_POMBFM
from bfm.bfm50.BFM50_rate_eqns import bfm50_rate_eqns
from pom_bfm_coupling.calculations import calculate_vertical_diffusivity
from pom_bfm_coupling.initialize_variables import set_initial_conditions

def pom_to_bfm(bfm_phys_vars, vertical_grid, temperature, salinity, inorganic_suspended_matter, shortwave_radiation, vertical_density_profile, wind_stress):

    """
    Description: Passes the physical variables to the BFM

    :return: seawater density, temperature and salinity, suspended sediment load,
             photosynthetically available radiation, gridpoint depth, and wind speed
    """

    # try:
    #     import NOPOINTERS
    #     NOPOINTERS = True
    # except FileNotFoundError:
    #     NOPOINTERS = False
    # if NOPOINTERS:
    #     from modules import ETW, ESW, EIR, ESS, ERHO, EWIND, Depth

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   1D ARRAYS FOR BFM
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    bfm_phys_vars.temperature = temperature.backward
    bfm_phys_vars.salinity = salinity.backward
    bfm_phys_vars.density = (vertical_density_profile * 1.E3) + 1.E3
    bfm_phys_vars.suspended_matter = inorganic_suspended_matter
    bfm_phys_vars.depth = vertical_grid.vertical_spacing * params_POMBFM.h

    bfm_phys_vars.irradiance = -1. * shortwave_radiation * water_specific_heat_times_density

    wind = np.sqrt(wind_stress.zonal**2 + wind_stress.meridional**2) * 1.E3
    bfm_phys_vars.wind = np.sqrt(wind/(1.25 * 0.0014))

    return bfm_phys_vars


def pom_bfm_1d(vertical_grid, time, diffusion, nutrients, bfm_phys_vars):

    bfm_variables = set_initial_conditions()
    bfm_rates = np.zeros((vertical_layers,50))

    dOdt_wind = np.zeros(vertical_layers)
    do3cdt_air_sea_flux = np.zeros(vertical_layers)
    
    for count in range(0,vertical_layers):
        conc = bfm_variables[count,:]
        bfm_rates[count,:], dOdt_wind[count], do3cdt_air_sea_flux[count] = bfm50_rate_eqns(count, bfm_variables, time, conc, seasonal_cycle=False)

    calculate_vertical_diffusivity(vertical_grid, diffusion, nutrients, bfm_variables, dOdt_wind, do3cdt_air_sea_flux)

