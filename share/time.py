import numpy as np
from cppdefs import *
from pom.modules_old import RLEN, bfm_lwp, LOGUNIT, seconds_per_day


class TimeInfo:

    def __init__(self, date0, time0, timeEnd, step0, timestep, stepnow, stepEnd):
        self.date0 = date0
        self.time0 = time0
        self.timeEnd = timeEnd
        self.step0 = step0
        self.timestep = timestep
        self.stepnow = stepnow
        self.stepEnd = stepEnd


# PUBLIC DATA MEMBERS:
timestr = ' ' * 19
start = '2000-01-01 00:00:00'
stop = ' ' * 19

# PRIVATE DATA MEMBERS
jul0 = -1
secs0 = -1


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: INITIALIZE THE TIME SYSTEM
#
# DESCRIPTION:  The subroutine {\tt init\_time()} initialises the time module by reading
#               a namelist and take actions according to the specifications.
#               On exit from this subroutine the two variables MinN and MaxN have well
#               defined values and can be used in the time loop.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def initialize_time_system(MinN,MaxN,timestep,timefmt,simdays,start,stop):

    """
    Description: Initializes the time module by reading a namelist.
                 Takes actions accrding to the specifications.

    :param MinN: minimum number of time steps
    :param MaxN: maximum number of time steps
    :param timestep: size of time step
    :param timefmt: time format
    :param simdays: number of days for simulation
    :param start: start date
    :param stop: end date
    :return:
    """

    jul1 = -1
    secs1 = -1

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   READ TIME SPECIFIC THINGS FROM THE NAMELIST
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    LEVEL1(); print('initialize_time_system')

    # CALCULATE MaxN -> MinN IS 1 IF NOT CHANGED BY HOTSTART
    MinN = 1
    LEVEL2(); print('Time step:      ', timestep, ' seconds')
    LEVEL2(); print('Time format:    ', timefmt)

    if timefmt == 1:
        HasRealTime = False
        LEVEL2(); print('# of timesteps: ', MaxN)
        # start = '2000-01-01 00:00:00'
        LEVEL2(); print('Fake start:     ', start)

    elif timefmt == 2:
        HasRealTime = True
        LEVEL2(); print('Start:          ', start)
        LEVEL2(); print('Stop:           ', stop)
        jul1, secs1 = read_time_string(start)
        jul2, secs2 = read_time_string(stop)

        nsecs = time_diff(jul2, secs2, jul1, secs1)
        MaxN = round(nsecs / timestep)

        ndays = jul2 - jul1
        if nsecs < 86400 and jul1 != jul2:
            ndays = ndays - 1

        nsecs = nsecs - 86400 * ndays
        string = '  ==> ' + str(ndays) + ' day(s) and ' + str(nsecs) + ' seconds ==> ' + str(MaxN) + ' time steps'
        STDERR(string)

    elif timefmt == 3:
        HasRealTime = True
        LEVEL2(); print('Start:          ', start)
        LEVEL2(); print('# of timesteps: ', MaxN)

        jul1, secs1 = read_time_string(start)

        nsecs = np.rint(MaxN * timestep) + secs1
        ndays = nsecs / 86400
        jul2 = jul1 + ndays
        secs2 = nsecs % 86400

        write_time_string(jul2, secs2)
        LEVEL2(); print('Stop:           ', stop)

    elif timefmt == 4:
        HasRealTime = False
        nsecs = simdays * 86400
        MaxN = np.rint(nsecs / timestep)
        LEVEL2(); print('# of timesteps: ', MaxN)
        # start = '2000-01-01 00:00:00'
        LEVEL2(); print('Fake start:     ', start)

    else:
        STDERR('Fatal error: A non valid input format has been chosen')
        return

    jul0 = jul1
    secs0 = secs1

    julianday = jul0
    secondsofday = secs0

    simtime = timestep * (MaxN - MinN + 1)

    # SET BFM TIME
    jday = float(jul0)
    yy, mm, dd, hh, nn = calendar_date(jday)

    bfmtime = TimeInfo(str(yy) + '-' + str(mm) + '-' + str(dd) + ' ' + str(hh) + ':' + str(nn),
                       jday,
                       jday + (float(MaxN) * timestep) / seconds_per_day,
                       MinN - 1,
                       timestep,
                       MinN - 1,
                       MaxN)

    LEVEL2(); print('bfmtime : ', bfmtime)

    return bfmtime


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: CONVERT TRUE JULIAN DAY TO CALENDAR DATE
#
# DESCRIPTION:  Converts a Julian day to a calendar date --- year, month and day.
#               Based on a similar routine in \emph{Numerical Recipes}.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def calendar_date(julian):

    """
    Description: Converts a Julian day to a calender date --- year, month, and day

    :param julian: Julian day
    :return: calender date
    """

    # LOCAL VARIABLES
    IGREG = 2299161

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    jday = np.floor(julian)
    if jday >= IGREG:
        x = ((jday - 1867216) - 0.25) / 36524.25
        ja = jday + 1 + int(x) - int(0.25 * x)
    else:
        ja = jday

    jb = ja + 1524
    jc = int(6680 + ((jb - 2439870) - 122.1) / 365.25)
    jd = int(365 * jc + (0.25 * jc))
    je = int((jb - jd) / 30.6001)

    dd = jb - jd - int(30.6001 * je)
    mm = je - 1

    if mm > 12:
        mm = mm - 12

    yyyy = jc - 4715

    if mm > 2:
        yyyy = yyyy - 1

    if yyyy <= 0:
        yyyy = yyyy - 1

    res = julian - float(jday)
    hh = np.floor(res * 24)
    nn = np.floor(((res * 24) - float(hh)) * 60)

    return yyyy, mm, dd, hh, nn


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: CONVERT CALENDAR DATE TO JULIAN DAY
#
# DESCRIPTION:  Converts a calendar date to a Julian day.
#               Based on a similar routine in \emph{Numerical Recipes}.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def julian_day(yyyy,mm,dd,hh,nn):

    """
    Description: Converts a calender date to a Julian day.

    :param yyyy: year
    :param mm: month
    :param dd: day
    :param hh: hour
    :param nn: minute
    :return: Julian day
    """

    # LOCAL VARIABLES
    IGREG = 15 + 31 * (10 + 12 * 1582)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    jy = yyyy
    if jy < 0:
        jy = jy + 1

    if mm > 2:
        jm = mm + 1
    else:
        jy = jy - 1
        jm = mm + 13

    jday = int(np.floor(365.25 * jy) + np.floor(30.6001 * jm) + dd + 1720995)
    if dd + 31 * (mm + 12 * yyyy) >= IGREG:
        ja = int(0.01 * jy)
        jday = jday + 2 - ja + int(0.25 * ja)

    jh = hh
    jn = nn

    if jn >= 60:
        jh = jh + 1
        jn = jn - 60

    if jh >= 24:
        jday = jday + 1
        jh = jh - 24

    julian = float(jday) + float(jh) / 24 + float(jn) / (24 * 60)

    return julian


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: KEEP TRACK OF TIME (JULIAN DAYS AND SECONDS)
#
# DESCRIPTION:  Based on a starting time this routine calculates the actual time
#               in a model integration using the number of time steps, {\tt n},
#               and the size of the time step, {\tt timestep}. More public variables
#               can be updated here if necessary.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def update_time(n,timestep,secs0):

    """
    Description: Based on a starting time this routine calculates the actual time in a model integration
                 using the number of time steps and the size of the time step.

    :param n: number of time steps
    :param timestep: size of time step
    :param secs0:
    :return: actual time in a model integration
    """

    nsecs = np.rint(n * timestep) + secs0
    fsecs = n * timestep + secs0
    julianday = jul0 + nsecs / 86400
    secondsofday = nsecs % 86400

    return fsecs, julianday, secondsofday


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: CONVERT A TIME STRING TO JULIAN DAY AND SECONDS
#
# DESCRIPTION:  Converts a time string to the true Julian day and seconds of that day.
#               The format of the time string must be: {\tt yyyy-mm-dd hh:hh:ss }.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def read_time_string(timestr):

    """
    Description: Converts a time string to the true Julian day and the seconds of that day.

    :param timestr: time string
    :return: Julian day and seconds of that day
    """

    yy = timestr[0:4]
    mm = timestr[5:7]
    dd = timestr[8:10]
    hh = timestr[11:13]
    mins = timestr[14:16]
    ss = timestr[17:19]

    julian, jday = julian_day(yy,mm,dd,0,0)
    jul = int(jday)
    secs = 3600*hh + 60*mins + ss

    return jul, secs


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: CONVERT JULIAN DAY AND SECONDS TO A TIME STRING
#
# DESCRIPTION:  Formats Julian day and seconds of that day to a nice looking
#               character string.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


def write_time_string(jul,secs):

    """
    Description: Formats Julian day and seconds of that day into a character string.

    :param jul: Julian day
    :param secs: seconds of Julian day
    :return: time string
    """

    jday = float(jul)
    yy, mm, dd, hh, nn = calendar_date(jday)

    hh = secs / 3600
    mins = (secs - hh * 3600) / 60
    ss = secs - 3600 * hh - 60 * mins

    timestr = str(yy) + '-' + str(mm) + '-' + str(dd) + ' ' + str(hh) + ':' + str(mins) + ':' + str(ss)

    return timestr


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: RETURN THE TIME DIFFERENCE IN SECONDS
#
# DESCRIPTION:  This functions returns the time difference between two
#               dates in seconds. The dates are given as Julian day and seconds
#               of that day.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def time_diff(jul2,secs2,jul1,secs1):

    """
    Description: Returns the time difference between two dates in seconds.

    :param jul2: Julian day 2
    :param secs2: seconds of julian day 2
    :param jul1: Julian day 1
    :param secs1: seconds of julian day 1
    :return: difference between Julian day 2 and Julian day 1
    """

    time_diff = 86400*(jul2-jul1) + (secs2-secs1)

    return time_diff


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: CONVERT A CALENDAR DATE TO TRUE JULIAN DAY
#
# DESCRIPTION:  Converts a Julian day to the day number in the current year
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def dayofyear(julian,yy):

    """
    Description: Converts a Julian day to the day number in the curent year.

    :param julian: Julian day
    :param yy: current year
    :return: day number of the current year
    """

    jday = float(julian)
    calendar_date(jday)
    julian0, jday = julian_day(yy,1,1,0,0)

    ddyear = julian - int(julian0) + 1

    return ddyear


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE:
#
# DESCRIPTION:
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def eomdays(Year,Month):

    """
    Description: Provides system with the number of days in a given month.
                 Accounts for leap year using Year input.

    :param Year: current year
    :param Month: current month
    :return: numbr of days in current month
    """

    eomdays = 0
    if Month < 1 or Month > 12:
        print('eomdays: Invalid Month!!')
    elif Month in [1,3,4,7,8,10,12]:  # Jan, Mar, May, July, Aug, Oct, Dec
        eomdays = 31
    elif Month == 2:  # Feb
        eomdays = 28
        if Year % 4 == 0:  # Leap Year
            eomdays = 29
    else:  # Apr, June, Sep, Nov
        eomdays = 30

    return eomdays


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE:
#
# DESCRIPTION:
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


def yeardays(Year):

    """
    Description: Provides system with the number of days in the current year.

    :param Year: current year
    :return: number of days in current year
    """

    yeardays = 0
    for i in range(0,12):
        yeardays = yeardays + float(eomdays(Year,i))

    if yeardays == 0 or yeardays > 366:
        print('yeardays out of bounds!!')

    return yeardays



# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   Copyright by the GOTM-team under the GNU Public License - www.gnu.org
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
