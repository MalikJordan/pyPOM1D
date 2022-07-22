import numpy as np
from bfm.variable_info import set_variable_info_bfm

variable_abbrev, variable_names, variable_units, index, pelagic_index, seaice_index, benthic_index = set_variable_info_bfm()

# def calculate_average_field(d3state,case):
#     ave_count = 0
#     do_3ave = True
#     if case == 'initialize':
#         i = 1
#     elif case == 'accumulate':
    
#         ave_count = ave_count + 1

#         if (pelagic_index.start >= 0) and (do_3ave):
#             k = 0
#             j = -1

#             for i in range(pelagic_index.state_start,pelagic_index.state_end):
#                 j = j+1






#     elif case == 'mean':
#         x = 1
    

# def calculate_daily_average_field(ave_count, case, day, d3source, d3ave):
def calculate_daily_average_field(d3ave,d3state,case):
    if case == 'Initialize':
        d3ave.single_day_ave[:,:] = 0
        # d3ave.count = d3ave.count + 1

    elif case == 'Accumulate':
        d3ave.single_day_ave[:,:] = d3ave.single_day_ave[:,:] + d3state[:,:]
        # d3ave.count = d3ave.count + 1

    elif case == 'Mean':
        d3ave.single_day_ave[:,:] = d3ave.single_day_ave[:,:] + d3state[:,:]
        d3ave.single_day_ave[:,:] = d3ave.single_day_ave[:,:]/d3ave.count
        d3ave.daily_ave[:,:,d3ave.day] = d3ave.single_day_ave[:,:]
        d3ave.monthly_ave[:,:,d3ave.month] = d3ave.monthly_ave[:,:,d3ave.month] + d3ave.daily_ave[:,:,d3ave.day]
        if (d3ave.day != 0) & ((d3ave.day+1) % 30 == 0):
            # should actually be when (d3ave.day+1) % 30 == 0 because of python '0 first' indexing
            d3ave.monthly_ave[:,:,d3ave.month] = d3ave.monthly_ave[:,:,d3ave.month]/30
            d3ave.month += 1

        d3ave.day = d3ave.day + 1
        d3ave.single_day_ave[:,:] = 0
    
    elif case == 'Reset':
        # Output writing and print functions will go here
        d3ave.count = 0
        d3ave.day = 0
        d3ave.month = 0
        d3ave.single_day_ave[:,:] = 0
        d3ave.daily_ave[:,:,:] = 0
        d3ave.monthly_ave[:,:,:] = 0
        


    return d3ave 