from inputs import params_POMBFM

current_path = '/Users/malikjordan/Desktop/pyPOM1D'
seconds_per_day = 86400.
earth_angular_velocity = 7.29E-5                                                        # OMEGA
vertical_layers = 151
DAYI = 1. / seconds_per_day
water_specific_heat_times_density = 4.187E6
twice_the_timestep = 2. * params_POMBFM.dti
