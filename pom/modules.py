import numpy as np
import cppdefs
from pom.initialize_variables import read_pom_input
from bfm.bfm_constants import seconds_per_day
from inputs import params_POMBFM
from inputs.pom_forcing_data import write_forcing_data

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

    # INITIALISATION AND FIRST FORCING READING
    if time_loop_counter == 0:

        # DAY COUNTER
        day_counter = 1

        # MONTH COUNTER
        month_counter = 0

        # TIME STEPS TO COVER ONE DAY
        timesteps_per_day = seconds_per_day / params_POMBFM.dti

        # TIME STEPS TO COVER ONE MONTH
        timesteps_per_month = 30 * timesteps_per_day

        # DAY INTERPOLATOR
        day_interpolator = -1
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
        month_interpolator = (timesteps_per_month / 2.) - 1
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
        sclim1, tclim1, wclim1, weddy1, weddy2, ism1, wsu1, wsv1, swrad1, wtsurf1, qcorr1, \
            NO3_s1, NH4_s1, PO4_s1, SIO4_s1, O2_b1, NO3_b1, PO4_b1, PON_b1 = write_forcing_data(month_counter)

        # UPDATE THE DAY COUNTER
        day_counter = day_counter + 1

        # UPDATE THE MONTH COUNTER
        month_counter = month_counter + 1

        # READING FOR SECOND MONTH
        sclim2, tclim2, wclim2, weddy3, weddy4, ism2, wsu2, wsv2, swrad2, wtsurf2, qcorr2, \
            NO3_s2, NH4_s2, PO4_s2, SIO4_s2, O2_b2, NO3_b2, PO4_b2, PON_b2 = write_forcing_data(month_counter)

    # UPDATE INTERPOLATION COUNTERS
    day_interpolator = day_interpolator + 1
    ratio_day = day_interpolator / timesteps_per_day

    month_interpolator = month_interpolator + 1
    ratio_month = month_interpolator / timesteps_per_month

    # INTERPOLATE WIND STRESS
    wusurf = wsu1 + ratio_month * (wsu2 - wsu1)
    wvsurf = wsv1 + ratio_month * (wsv2 - wsv1)

    # INTERPOLATE HEAT FLUX
    if params_POMBFM.idiagn == 0:
        wtsurf = wtsurf1 + ratio_month * (wtsurf2 - wtsurf1)
        swrad = swrad1 + ratio_month * (swrad2 - swrad1)
    elif params_POMBFM.idiagn == 1:
        # DAILY
        # surface_solar_radiation = shortwave_radiation1 + ratio_day * (shortwave_radiation2 - shortwave_radiation1)
        # MONTHLY
        swrad = swrad1 + ratio_month * (swrad2 - swrad1)

    # INTERPOLATE T&S PROFILES
    tstar = np.zeros(vertical_layers)
    sstar = np.zeros(vertical_layers)
    wgen  = np.zeros(vertical_layers)
    weddy = np.zeros(vertical_layers)
    tf    = np.zeros(vertical_layers)
    sf    = np.zeros(vertical_layers)

    tstar[:] = tclim1[:] + ratio_month * (tclim2[:] - tclim1[:])
    sstar[:] = sclim1[:] + ratio_month * (sclim2[:] - sclim1[:])
    wgen[:]  = wclim1[:] + ratio_month * (wclim2[:] - wclim1[:])

    if ratio_month <= 0.5:
        weddy[:] = weddy1[:]
    else:
        weddy[:] = weddy2[:]

    if params_POMBFM.idiagn == 0:
        tsurf = tstar[0]
        ssurf = sstar[0]
    elif params_POMBFM.idiagn == 1:
        tf[:] = tstar[:]
        sf[:] = sstar[:]

    # INTERPOLATE SUSPENDED INORGANIC MATTER
    ism = np.zeros(vertical_layers-1)
    ism[:] = ism1[:] + ratio_month * (ism2[:] - ism1[:])

    # INTERPOLATE SURFACE NUTRIENTS
    NO3surf = NO3_s1 + ratio_month * (NO3_s2 - NO3_s1)
    NH4surf = NH4_s1 + ratio_month * (NH4_s2 - NH4_s1)
    PO4surf = PO4_s1 + ratio_month * (PO4_s2 - PO4_s1)
    SIO4surf = SIO4_s1 + ratio_month * (SIO4_s2 - SIO4_s1)

    # INTERPOLATE BOTTOM NUTRIENTS
    O2bott = O2_b1 + ratio_month * (O2_b2 - O2_b1)
    NO3bott = NO3_b1 + ratio_month * (NO3_b2 - NO3_b1)
    PO4bott = PO4_b1 + ratio_month * (PO4_b2 - PO4_b1)
    PONbott_grad = PON_b1 + ratio_month * (PON_b2 - PON_b1)

    if month_interpolator == timesteps_per_month:

        # A MONTH HAS GONE...IT IS NECESSARY TO...
        # ....UPDATE MONTH COUNTER....
        month_counter = month_counter + 1
        print('month_counter = ',month_counter)

        # ....RESET INTERPOLATOR....
        month_interpolator = 0

        # ....SHIFT THE MONTHLY DATA....
        wsu1 = wsu2
        wsv1 = wsv2
        swrad1 = swrad2
        wtsurf1 = wtsurf2
#         SLUX1 = SLUX2
        NO3_s1 = NO3_s2
        NH4_s1 = NH4_s2
        PO4_s1 = PO4_s2
        SIO4_s1 = SIO4_s2
        NO3_b1 = NO3_b2
        O2_b1 = O2_b2
        PO4_b1 = PO4_b2
        PON_b1 = PON_b2
        ism1[:]     = ism2[:]
        tclim1[:]   = tclim2[:]
        sclim1[:]   = sclim2[:]
        wclim1[:]   = wclim2[:]
        weddy1[:]   = weddy3[:]
        weddy2[:]   = weddy4[:]

        # IF 12 MONTHS HAVE GONE, RESTART THE READING SEQUENCE
        if month_counter > 12:
            month_counter = 0
            sclim1, tclim1, wclim1, weddy1, weddy2, ism1, wsu1, wsv1, swrad1, wtsurf1, qcorr1, \
                NO3_s1, NH4_s1, PO4_s1, SIO4_s1, O2_b1, NO3_b1, PO4_b1, PON_b1 = write_forcing_data(month_counter)

            month_counter = month_counter + 1

        # READ FOLLOWING MONTH
        sclim2, tclim2, wclim2, weddy3, weddy4, ism2, wsu2, wsv2, swrad2, wtsurf2, qcorr2, \
            NO3_s2, NH4_s2, PO4_s2, SIO4_s2, O2_b2, NO3_b2, PO4_b2, PON_b2 = write_forcing_data(month_counter + 1)

    return

