import numpy as np
from pom.constants import vertical_layers, water_specific_heat_times_density
from inputs import params_POMBFM
from bfm.bfm50.BFM50_rate_eqns import bfm50_rate_eqns
from pom_bfm_coupling.calculations import calculate_vertical_diffusivity

def pom_to_bfm(vertical_grid, temperature, salinity, inorganic_suspended_matter, shortwave_radiation, vertical_density_profile, wind_stress):

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

    etw = temperature.backward
    esw = salinity.backward
    erho = (vertical_density_profile * 1.E3) + 1.E3
    ess = inorganic_suspended_matter
    depth = vertical_grid.vertical_spacing * params_POMBFM.h

    eir = np.zeros(vertical_layers - 1)
    eir[0] = -1. * shortwave_radiation * water_specific_heat_times_density

    wind = np.sqrt(wind_stress.zonal**2 + wind_stress.meridional**2) * 1.E3
    ewind = np.sqrt(wind/(1.25 * 0.0014))

    return etw, esw, eir, ess, erho, depth, ewind


def pom_bfm_1d(vertical_grid, i, time, diffusion, nutrients, temperature, salinity, inorganic_suspended_matter,
               shortwave_radiation, vertical_density_profile, wind_stress):

    etw, esw, eir, ess, erho, depth, ewind = pom_to_bfm(vertical_grid, temperature, salinity, inorganic_suspended_matter,
                                                        shortwave_radiation, vertical_density_profile, wind_stress)

    bfm50_rate_eqns(time, conc, seasonal_cycle=False)

    calculate_vertical_diffusivity(vertical_grid, diffusion, nutrients, o2o, n1p, n3n)

# import numpy as np
#
#
# def environmental_forcing():
#
#     # Pass physical variables into bfm
#     # pom_to_bfm()
#
#     # Calculation of vertical extinction coefficient
#
#     # Calculation of the irradiation (forcing function)
