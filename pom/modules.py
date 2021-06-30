import numpy as np
import cppdefs
from pom.initialize_variables import read_pom_input
from bfm.bfm_constants import seconds_per_day
from inputs import params_POMBFM
from inputs.pom_forcing_data import write_forcing_data
from pom.data_classes import MonthlyForcingData, Stresses

# From ModulePom

# SPECIFIC HEAT TIMES RHO0
water_specific_heat_times_density = 4.187E6

# 1 DAY IN SECONDS (RECIPROCAL)
DAYI = 1. / seconds_per_day

# VERTICAL LAYERS
vertical_layers = 151


def forcing_manager(time_loop_counter,counters,month1_data,month2_data):

    # INITIALISATION AND FIRST FORCING READING
    if time_loop_counter == 0:

        # DAY COUNTER
        counters.day_counter = 1

        # MONTH COUNTER
        counters.month_counter = 0

        # TIME STEPS TO COVER ONE DAY
        counters.timesteps_per_day = seconds_per_day / params_POMBFM.dti

        # TIME STEPS TO COVER ONE MONTH
        counters.timesteps_per_month = 30 * counters.timesteps_per_day

        # DAY INTERPOLATOR
        counters.day_interpolator = -1
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

        # MONTH INTERPOLATOR
        counters.month_interpolator = (counters.timesteps_per_month / 2.) - 1
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

        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # **                                                                   **
        # **  INITIAL READING OF THE MONTHLY FORCING                           **
        # **                                                                   **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
        # ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **

        # READING FOR FIRST MONTH
        month1_data = write_forcing_data(counters.month_counter)

        # UPDATE THE DAY COUNTER
        counters.day_counter = counters.day_counter + 1

        # UPDATE THE MONTH COUNTER
        counters.month_counter = counters.month_counter + 1

        # READING FOR SECOND MONTH
        month2_data = write_forcing_data(counters.month_counter)

    # UPDATE INTERPOLATION COUNTERS
    counters.day_interpolator = counters.day_interpolator + 1
    counters.ratio_day = counters.day_interpolator / counters.timesteps_per_day

    counters.month_interpolator = counters.month_interpolator + 1
    counters.ratio_month = counters.month_interpolator / counters.timesteps_per_month

    # INTERPOLATE WIND STRESS
    wusurf = month1_data.wsu + counters.ratio_month * (month2_data.wsu - month1_data.wsu)
    wvsurf = month1_data.wsv + counters.ratio_month * (month2_data.wsv - month1_data.wsv)
    wind_stress = Stresses(wusurf,wvsurf)

    # INTERPOLATE HEAT FLUX
    if params_POMBFM.idiagn == 0:
        wtsurf = month1_data.wtsurf + counters.ratio_month * (month2_data.wtsurf - month1_data.wtsurf)
        swrad = month1_data.swrad + counters.ratio_month * (month2_data.swrad - month1_data.swrad)
    else:
        wtsurf = 0  # not needed for diagnostic mode (idiagn = 1), see 4.6.5 in manual
        swrad = month1_data.swrad + counters.ratio_month * (month2_data.swrad - month1_data.swrad)



    # INTERPOLATE T&S PROFILES
    tstar = np.zeros(vertical_layers)
    sstar = np.zeros(vertical_layers)
    wgen  = np.zeros(vertical_layers)
    weddy = np.zeros(vertical_layers)
    tf    = np.zeros(vertical_layers)
    sf    = np.zeros(vertical_layers)

    tstar[:] = month1_data.tclim[:] + counters.ratio_month * (month2_data.tclim[:] - month1_data.tclim[:])
    sstar[:] = month1_data.sclim[:] + counters.ratio_month * (month2_data.sclim[:] - month1_data.sclim[:])
    wgen[:]  = month1_data.wclim[:] + counters.ratio_month * (month2_data.wclim[:] - month1_data.wclim[:])

    if counters.ratio_month <= 0.5:
        weddy[:] = month1_data.weddy1[:]
    else:
        weddy[:] = month1_data.weddy2[:]

    if params_POMBFM.idiagn == 0:
        tsurf = tstar[0]
        ssurf = sstar[0]
    elif params_POMBFM.idiagn == 1:
        tf[:] = tstar[:]
        sf[:] = sstar[:]

    # INTERPOLATE SUSPENDED INORGANIC MATTER
    ism = np.zeros(vertical_layers-1)
    ism[:] = month1_data.ism[:] + counters.ratio_month * (month2_data.ism[:] - month1_data.ism[:])

    # INTERPOLATE SURFACE NUTRIENTS
    NO3surf = month1_data.NO3_s + counters.ratio_month * (month2_data.NO3_s - month1_data.NO3_s)
    NH4surf = month1_data.NH4_s + counters.ratio_month * (month2_data.NH4_s - month1_data.NH4_s)
    PO4surf = month1_data.PO4_s + counters.ratio_month * (month2_data.PO4_s - month1_data.PO4_s)
    SIO4surf = month1_data.SIO4_s + counters.ratio_month * (month2_data.SIO4_s - month1_data.SIO4_s)

    # INTERPOLATE BOTTOM NUTRIENTS
    O2bott = month1_data.O2_b + counters.ratio_month * (month2_data.O2_b - month1_data.O2_b)
    NO3bott = month1_data.NO3_b + counters.ratio_month * (month2_data.NO3_b - month1_data.NO3_b)
    PO4bott = month1_data.PO4_b + counters.ratio_month * (month2_data.PO4_b - month1_data.PO4_b)
    PONbott_grad = month1_data.PON_b + counters.ratio_month * (month2_data.PON_b - month1_data.PON_b)

    if counters.month_interpolator == counters.timesteps_per_month:

        # A MONTH HAS GONE...IT IS NECESSARY TO...
        # ....UPDATE MONTH COUNTER....
        counters.month_counter = counters.month_counter + 1
        print('month_counter = ',counters.month_counter)

        # ....RESET INTERPOLATOR....
        counters.month_interpolator = 0

        # ....SHIFT THE MONTHLY DATA....
        month1_data = month2_data

        # IF 12 MONTHS HAVE GONE, RESTART THE READING SEQUENCE
        if counters.month_counter > 12:
            counters.month_counter = 0
            month1_data = write_forcing_data(counters.month_counter)

            counters.month_counter = counters.month_counter + 1

        # READ FOLLOWING MONTH
        month2_data = write_forcing_data(counters.month_counter)

    return tf, tstar, sf, sstar, swrad, wtsurf, wind_stress, month1_data, month2_data, counters

