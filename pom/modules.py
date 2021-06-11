import numpy as np
import cppdefs
from pom.initialize_variables import read_pom_input
from bfm.bfm_constants import seconds_per_day
from inputs import params_POMBFM

# From ModulePom

# SPECIFIC HEAT TIMES RHO0
water_specific_heat_times_density = 4.187E6

# 1 DAY IN SECONDS (RECIPROCAL)
DAYI = 1. / seconds_per_day

# VERTICAL LAYERS
vertical_layers = 151

def forcing_manager(time_loop_counter):

    # LENGTH OF INPUT ARRAYS
    array_length = 13

    # LOOP COUNTER
    # K = float()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL ARRAYS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # RLENGTH = float()

    # INITIALISATION AND FIRST FORCING READING
    if time_loop_counter == 0:

        # DOUBLE READING OF DATA (NEEDED TO CARRY OUT THE TIME LINEAR INTERPOLATION)
        wind_speed_zonal1, wind_speed_meridional1, surface_salinity1, solar_radiation1, inorganic_suspended_matter1, \
            salinity_climatology1, temperature_climatology1, w_velocity_climatology1, w_eddy_velocity_1, \
            w_eddy_velocity_2, salinity_backward1, temperature_backward1, \
            shortwave_radiation1, surface_heat_flux1, kinetic_energy_loss1, \
            NO3_s1, NH4_s1, PO4_s1, SIO4_s1, O2_b1, NO3_b1, PO4_b1, PON_b1                          = read_pom_input()

        # NOT ALL VARIABLES FROM SECOND CALLING WILL BE USED
        wind_speed_zonal2, wind_speed_meridional2, surface_salinity2, solar_radiation2, inorganic_suspended_matter2, \
            salinity_climatology2, temperature_climatology2, w_velocity_climatology2, w_eddy_velocity_3, \
            w_eddy_velocity_4, salinity_backward2, temperature_backward2, \
            shortwave_radiation2, surface_heat_flux2, kinetic_energy_loss2, \
            NO3_s2, NH4_s2, PO4_s2, SIO4_s2, O2_b2, NO3_b2, PO4_b2, PON_b2                          = read_pom_input()

        # DAY COUNTER
        day_counter = 1

        # MONTH COUNTER
        month_counter = 1

        # TIME STEPS TO COVER ONE DAY
        timesteps_per_day = int(seconds_per_day) / int(params_POMBFM.dti)

        # TIME STEPS TO COVER ONE MONTH
        timesteps_per_month = 30 * timesteps_per_day

        # DAY INTERPOLATOR

        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # **                                                                   **
        # **  THE DAILY CLIMATOLOGICAL FORCING DATA ARE ASSUMED TO BE          **
        # **  CENTERED AT h 00.00 OF EACH CLIMATOLOGICAL DAY.THEREFORE         **
        # **  THE MONTH INTERPOLATOR(day_interpolator) IS INITIALISED AT THE VALUE       **
        # **  CORRESPONDING TO MIDNIGHT MINUS 1 TIMESTEP.                      **
        # **                                                                   **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **

        day_interpolator = -1

        # MONTH INTERPOLATOR

        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # **                                                                   **
        # **  THE MONTHLY CLIMATOLOGICAL FORCING DATA ARE ASSUMED TO BE        **
        # **  CENTERED AT DAY 15 OF EACH CLIMATOLOGICAL MONTH. THEREFORE       **
        # **  THE MONTH INTERPOLATOR (month_interpolator) IS INITIALISED AT THE VALUE       **
        # **  (timesteps_per_month/2)-1 CORRESPONDING TO DAY 15 MINUS 1 TIMESTEP.           **
        # **                                                                   **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **

        month_interpolator = (timesteps_per_month / 2.) - 1

        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # **                                                                   **
        # **  INITIAL READING OF THE MONTHLY FORCING                           **
        # **                                                                   **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **

        # WIND STRESS CONVERTED TO POM UNITS (N/m2-->m2/s2)
        wind_speed_zonal1[:] = -wind_speed_zonal1[:]  * 1.E-03
        wind_speed_zonal2[:]  = -wind_speed_zonal2[:]  * 1.E-03
        wind_speed_meridional1[:]  = -wind_speed_meridional1[:]  * 1.E-03
        wind_speed_meridional2[:]  = -wind_speed_meridional2[:]  * 1.E-03

        # HEAT FLUX CONVERTED TO POM UNITS(W/m2-->deg.C*m/s)
        shortwave_radiation1[:] = -shortwave_radiation1[:] / water_specific_heat_times_density
        shortwave_radiation2[:] = -shortwave_radiation2[:] / water_specific_heat_times_density
        surface_heat_flux1[:] = -surface_heat_flux1[:] / water_specific_heat_times_density
        surface_heat_flux2[:] = -surface_heat_flux2[:] / water_specific_heat_times_density

        # VERTICAL VELOCITY CONVERTED TO POM UNITS (m/s-->m/s)
        # w_velocity_climatology1[:] = w_velocity_climatology1[:]
        # w_velocity_climatology2[:] = w_velocity_climatology2[:]

        # UPDATE THE DAY COUNTER
        day_counter = day_counter + 1

        # UPDATE THE MONTH COUNTER
        month_counter = month_counter + 1

    # UPDATE INTERPOLATION COUNTERS
    day_interpolator = day_interpolator + 1
    ratio_day = day_interpolator / timesteps_per_day

    month_interpolator = month_interpolator + 1
    ratio_month = month_interpolator / timesteps_per_month

    # INTERPOLATE WIND STRESS
    wind_stress_zonal = wind_speed_zonal1 + ratio_month * (wind_speed_zonal2 - wind_speed_zonal1)
    wind_stress_meridional = wind_speed_meridional1 + ratio_month * (wind_speed_meridional2 - wind_speed_meridional1)

    # INTERPOLATE HEAT FLUX
    if params_POMBFM.idiagn == 0:
        surface_heat_flux_loss = surface_heat_flux1 + ratio_month * (surface_heat_flux2 - surface_heat_flux1)
        surface_solar_radiation = shortwave_radiation1 + ratio_month * (shortwave_radiation2 - shortwave_radiation1)
    elif params_POMBFM.idiagn == 1:
        # DAILY
        # surface_solar_radiation = shortwave_radiation1 + ratio_day * (shortwave_radiation2 - shortwave_radiation1)
        # MONTHLY
        surface_solar_radiation = shortwave_radiation1 + ratio_month * (shortwave_radiation2 - shortwave_radiation1)

    # INTERPOLATE T&S PROFILES
    interpolated_temperature = np.zeros((array_length, vertical_layers), dtype=float)
    interpolated_salinity = np.zeros((array_length, vertical_layers), dtype=float)
    interpolated_temperature[:] = temperature_climatology1[:] + \
                                           ratio_month * (temperature_climatology2[:] - temperature_climatology1[:])
    interpolated_salinity[:] = salinity_climatology1[:] + \
                                        ratio_month * (salinity_climatology2[:] - salinity_climatology1[:])
    interpolated_w_velocity[:]  = w_velocity_climatology1[:] + ratio_month * (w_velocity_climatology2[:] - w_velocity_climatology1[:])

    w_eddy_velocity = np.zeros((array_length, vertical_layers), dtype=float)
    if ratio_month <= 0.5:
        w_eddy_velocity[:] = w_eddy_velocity_1[:]
    else:
        w_eddy_velocity[:] = w_eddy_velocity_2[:]

    if params_POMBFM.idiagn == 0:
        surface_temperature = interpolated_temperature[0]
        surface_salinity = interpolated_salinity[0]
    elif params_POMBFM.idiagn == 1:
        temperature_forward[:] = interpolated_temperature[:]
        salinity_forward[:] = interpolated_salinity[:]

    # INTERPOLATE SUSPENDED INORGANIC MATTER
    inorganic_suspended_matter = np.zeros((array_length, vertical_layers), dtype=float)
    inorganic_suspended_matter[:] = inorganic_suspended_matter1[:] + ratio_month * (inorganic_suspended_matter2[:] - inorganic_suspended_matter1[:])

    # INTERPOLATE SURFACE NUTRIENTS
    NO3SURF = NO3_s1 + ratio_month * (NO3_s2 - NO3_s1)
    NH4SURF = NH4_s1 + ratio_month * (NH4_s2 - NH4_s1)
    PO4SURF = PO4_s1 + ratio_month * (PO4_s2 - PO4_s1)
    SIO4SURF = SIO4_s1 + ratio_month * (SIO4_s2 - SIO4_s1)

    # INTERPOLATE BOTTOM NUTRIENTS
    O2BOTT = O2_b1 + ratio_month * (O2_b2 - O2_b1)
    NO3BOTT = NO3_b1 + ratio_month * (NO3_b2 - NO3_b1)
    PO4BOTT = PO4_b1 + ratio_month * (PO4_b2 - PO4_b1)
    PONBOTTgrad = PON_b1 + ratio_month * (PON_b2 - PON_b1)

    if month_interpolator == timesteps_per_month:

        # A MONTH HAS GONE...IT IS NECESSARY TO...
        # ....UPDATE MONTH COUNTER....
        month_counter = month_counter + 1
        print('month_counter = ',month_counter)

        # ....RESET INTERPOLATOR....
        month_interpolator = 0

        # ....SHIFT THE MONTHLY DATA....
        wind_speed_zonal1 = wind_speed_zonal2
        wind_speed_meridional1 = wind_speed_meridional2
        shortwave_radiation1 = shortwave_radiation2
        surface_heat_flux1 = surface_heat_flux2
#         SLUX1 = SLUX2
        NO3_s1 = NO3_s2
        NH4_s1 = NH4_s2
        PO4_s1 = PO4_s2
        SIO4_s1 = SIO4_s2
        NO3_b1 = NO3_b2
        O2_b1 = O2_b2
        PO4_b1 = PO4_b2
        PON_b1 = PON_b2
        inorganic_suspended_matter1[:]     = inorganic_suspended_matter2[:]
        temperature_climatology1[:]   = temperature_climatology2[:]
        salinity_climatology1[:]   = salinity_climatology2[:]
        w_velocity_climatology1[:]   = w_velocity_climatology2[:]
        w_eddy_velocity_1[:]   = w_eddy_velocity_3[:]
        w_eddy_velocity_2[:]   = w_eddy_velocity_4[:]

        if month_counter > 13:
            # IF 12 MONTHS HAVE GONE, RESTART THE READING SEQUENCE
            month_counter = 2

            wind_speed_zonal1, wind_speed_meridional1, surface_salinity1, solar_radiation1, inorganic_suspended_matter1, \
                salinity_climatology1, temperature_climatology1, w_velocity_climatology1, w_eddy_velocity_1, \
                w_eddy_velocity_2, salinity_backward1, temperature_backward1, \
                shortwave_radiation1, surface_heat_flux1, kinetic_energy_loss1, \
                NO3_s1, NH4_s1, PO4_s1, SIO4_s1, O2_b1, NO3_b1, PO4_b1, PON_b1                       = read_pom_input()
            wind_speed_zonal1 = -wind_speed_zonal1 * 1.E-03
            wind_speed_meridional1 = -wind_speed_meridional1 * 1.E-03
            shortwave_radiation1 = -shortwave_radiation1 / water_specific_heat_times_density
            surface_heat_flux1 = -surface_heat_flux1 / water_specific_heat_times_density

            wind_speed_zonal2, wind_speed_meridional2, surface_salinity2, solar_radiation2, inorganic_suspended_matter2, \
                salinity_climatology2, temperature_climatology2, w_velocity_climatology2, w_eddy_velocity_3, \
                w_eddy_velocity_4, salinity_backward2, temperature_backward2, \
                shortwave_radiation2, surface_heat_flux2, kinetic_energy_loss2, \
                NO3_s2, NH4_s2, PO4_s2, SIO4_s2, O2_b2, NO3_b2, PO4_b2, PON_b2                       = read_pom_input()

            wind_speed_zonal2 = -wind_speed_zonal2 * 1.E-03
            wind_speed_meridional2 = -wind_speed_meridional2 * 1.E-03
            shortwave_radiation2 = -shortwave_radiation2 / water_specific_heat_times_density
            surface_heat_flux2 = -surface_heat_flux2 / water_specific_heat_times_density

    return

