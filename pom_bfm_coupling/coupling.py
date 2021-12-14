import numpy as np
from pom.constants import vertical_layers, water_specific_heat_times_density
from inputs import params_POMBFM
from bfm.bfm50.BFM50_rate_eqns import bfm50_rate_eqns
from pom_bfm_coupling.calculations import calculate_vertical_diffusivity
from pom_bfm_coupling.initialize_variables import set_initial_conditions

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

    [disOxygen_IO_O, phospate_IO_P, nitrate_IO_N, ammonium_IO_N, o4n, silicate_IO_Si, reductEquiv_IO_R, pelBacteria_LO_C, pelBacteria_LO_N, pelBacteria_LO_P,
     diatoms_LO_C, diatoms_LO_N, diatoms_LO_P, diatoms_LO_Chl, diatoms_LO_Si, nanoflagellates_LO_C, nanoflagellates_LO_N, nanoflagellates_LO_P, nanoflagellates_LO_Chl, picophyto_LO_C, picophyto_LO_N, picophyto_LO_P, picophyto_LO_Chl, largephyto_LO_C, largephyto_LO_N, largephyto_LO_P, largephyto_LO_Chl,
     carnivMesozoo_LO_C, carnivMesozoo_LO_N, carnivMesozoo_LO_P, omnivMesozoo_LO_C, omnivMesozoo_LO_N, omnivMesozoo_LO_P, microzoo_LO_C, microzoo_LO_N, microzoo_LO_P, z6c, z6n, z6p,
     labileDOM_NO_C, labileDOM_NO_N, labileDOM_NO_P, semilabileDOC_NO_C, semirefractDOC_NO_C, particOrganDetritus_NO_C, particOrganDetritus_NO_N, particOrganDetritus_NO_P, particOrganDetritus_NO_Si, disInorgCarbon_IO_C, o3h] = set_initial_conditions()

    bfm50_rate_eqns(time, conc, seasonal_cycle=False)

    calculate_vertical_diffusivity(vertical_grid, diffusion, nutrients, disOxygen_IO_O, phospate_IO_P, nitrate_IO_N)

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
