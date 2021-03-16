import numpy as np
from cppdefs import *
from pom.modules import RLEN, bfm_lwp, LOGUNIT, seconds_per_day


# TimeInfo = dict(date0=' '*25,       # CALENDER DATE OF START (EMPTY STRING OF LENGTH 25)
#                 time0=float(),      # JULIAN DAY START OF RUN
#                 timeEnd=float(),    # JULIAN DAY END OF RUN
#                 step0=int(),        # INITIAL STEP #
#                 timestep=int(),     # DELTA t
#                 stepnow=int(),      # ACTUAL STEP #
#                 stepEnd=int())      # ACTUAL STEP #


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
timestep = float()
fsecs = float()
simtime = float()
julianday = int()
secondsofday = int()
timefmt = int()
simdays = int()
# MinN = int()
# MaxN = int()
HasRealTime = True

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

def initialize_time_system(MinN, MaxN):

    bfmtime = TimeInfo
    # INPUT/OUTPUT PARAMETERS
    # MinN = int()
    # MaxN = int()

    # LOCAL VARIABLES
    jul1 = -1
    secs1 = -1
    jul2 = int()
    secs2 = int()
    ndays = int()
    nsecs = int()
    dd = int()
    mm = int()
    yy = int()
    hh = int()
    nn = int()
    jday = float

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
        LEVEL2(); print('Start:          ', start)
        LEVEL2(); print('Stop:           ', stop)
        read_time_string(start, jul1, secs1)
        read_time_string(stop, jul2, secs2)

        nsecs = time_diff(jul2, secs2, jul1, secs1)
        MaxN = round(nsecs / timestep)

        ndays = jul2 - jul1
        if nsecs < 86400 and jul1 != jul2:
            ndays = ndays - 1

        nsecs = nsecs - 86400 * ndays
        string = '  ==> ', ndays, ' day(s) and ', nsecs, ' seconds ==> ', MaxN, ' time steps'
        STDERR(string)

    elif timefmt == 3:
        LEVEL2(); print('Start:          ', start)
        LEVEL2(); print('# of timesteps: ', MaxN)

        read_time_string(start, jul1, secs1)

        nsecs = np.rint(MaxN * timestep) + secs1
        ndays = nsecs / 86400
        jul2 = jul1 + ndays
        secs2 = nsecs % 86400

        write_time_string(jul2, secs2, stop)
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
    yy, mm, dd, hh, nn = calendar_date(jday, yy, mm, dd, hh, nn)

    bfmtime.date0 = yy, '-', mm, '-', dd, ' ', hh, ':', nn
    bfmtime.time0 = jday
    bfmtime.timeEnd = jday + (float(MaxN) * timestep) / seconds_per_day
    bfmtime.step0 = MinN - 1
    bfmtime.timestep = timestep
    bfmtime.stepnow = MinN - 1
    bfmtime.stepEnd = MaxN

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

def calendar_date(julian, yyyy, mm, dd, hh, nn):
    # INPUT PARAMETERS
    # julian = float()

    # OUTPUT PARAMETERS
    yyyy = int()
    mm = int()
    dd = int()
    hh = int()
    nn = int()

    # LOCAL VARIABLES
    IGREG = 2299161
    ja = int()
    jb = int()
    jc = int()
    jd = int()
    je = int()
    jday = int()
    x = float()
    res = float()

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
# # ROUTINE: CONVERT CALENDAR DATE TRUE JULIAN DAY
#
# DESCRIPTION:  Converts a calendar date to a Julian day.
#               Based on a similar routine in \emph{Numerical Recipes}.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def julian_day(yyyy,mm,dd,hh,nn,julian):

    # INPUT PARAMETERS
    # yyyy = int()
    # mm = int()
    # dd = int()
    # hh = int()
    # nn = int()

    # OUTPUT PARAMETERS
    # julian = float()

    # LOCAL VARIABLES
    IGREG = 15 + 31 * (10 + 12 * 1582)
    ja = int()
    jy = int()
    jm = int()
    jh = int()
    jn = int()
    jday = int()

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

    return julian, jday


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

def update_time(n):

    # INPUT PARAMETERES
    # n = int()

    # LOCAL VARIABLES
    nsecs = int()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    nsecs = np.rint(n * timestep) + secs0
    fsecs = n * timestep + secs0
    julianday = jul0 + nsecs / 86400
    secondsofday = nsecs % 86400

    return


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: CONVERT A TIME STRING TO JULINA DAY AND SECONDS
#
# DESCRIPTION:  Converts a time string to the true Julian day and seconds of that day.
#               The format of the time string must be: {\tt yyyy-mm-dd hh:hh:ss }.
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def read_time_string(timestr,jul,secs):

    # INPUT PARAMETERS
    # timestr = ' '*19

    # OUTPUT PARAMETERS
    jul = int()
    secs = int()

    # LOCAL VARIABLES
    yy = int()
    mm = int()
    dd = int()
    hh = int()
    mins = int()
    ss = int()
    jday = float()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    yy = timestr[0:4]
    mm = timestr[5:7]
    dd = timestr[8:10]
    hh = timestr[11:13]
    mins = timestr[14:16]
    ss = timestr[17:19]

    julian, jday = julian_day(yy,mm,dd,0,0,jday)
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


def write_time_string(jul,secs,timestr):

    # INPUT PARAMETERS
    # jul = int()
    # secs = int()

    # OUTPUT PARAMETERS
    timestr = ' '*19

    # LOCAL VARIABLES
    ss = int()
    mins = int()
    hh = int()
    dd = int()
    mm = int()
    yy = int()
    jh = int()
    jn = int()
    jday = float()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    hh = secs / 3600
    mins = (secs - hh * 3600) / 60
    ss = secs - 3600 * hh - 60 * mins
    jday = float(jul)

    calendar_date(jday,yy,mm,dd,jh,jn)

    timestr = yy,'-',mm,'-',dd,' ',hh,':',mins,':',ss

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

def time_diff(jul1,secs1,jul2,secs2):

    # INPUT PARAMETERS
    # jul1 = int()
    # secs1 = int()
    # jul2 = int()
    # secs2 = int()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    time_diff = 86400*(jul1-jul2) + (secs1-secs2)

    return time_diff


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: CONVERT A CALENDAR DATE TO TRUE JULIAN DAY
#
# DESCRIPTION:  Converts a Julian day to the day number in the current year
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def dayofyear(julian,ddyear):

    # INPUT PARAMETERS
    # julian = int()

    # OUTPUT PARAMETERS
    ddyear = int()

    # LOCAL VARIABLES
    yy = int()
    mm = int()
    dd = int()
    hh = int()
    nn = int()
    julian0 = float()
    jday = float()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    jday = float(julian)
    calendar_date(jday,yy,mm,dd,hh,nn)
    julian0, jday = julian_day(yy,1,1,0,0,julian0)

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

    # Month = int()
    # Year = int()

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

    im = int()
    # Year = int()
    yeardays = 0
    for im in range(0,12):
        yeardays = yeardays + float(eomdays(Year,im))

    if yeardays == 0 or yeardays > 366:
        print('yeardays out of bounds!!')

    return yeardays



# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   Copyright by the GOTM-team under the GNU Public License - www.gnu.org
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
