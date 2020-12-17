# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  BFM - Biogeochemical Flux Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: ModuleGlobalMem
#
# DESCRIPTION
#
#   Definiton of the runtime error messages.
#   This module contains global settings:
#   -general constants for controlling prescision,
#   -parameters defining fle streams and error message numbers
#   -the subroutine for printing the message
#   and aborting the simulation
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from decimal import *

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBAL CONSTANTS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# DOUBLE PRECISION
getcontext().prec = 12
RLEN = Decimal()
NMLUNIT = 310

# THE UNIT OF THE LOG FILE IS NOT A PARAMETER TO ALLOW PARALLEL WRITING
LOGUNIT = 0
bef_lwp = True
ZERO = Decimal(0.0)
ONE = Decimal(1.0)
PI = Decimal(3.14159265359)
BASETEMP = Decimal(20.0)

#   NEXT PARAMETERS ARE DEFINED TO CONTROL TYPES OF STATE VARIABLES
OFF = -100
SINKSOURCE = -1
NOTRANSPORT = 0
NOOBCSTATES = 1
HORTRANSPORT = 10
ALLTRANSPORT = 20
DONE = Decimal(1.0)

#   ERROR CODES:
ALLOC = 10
NML_OPEN = 11
NML_READ = 12
DIM_MISMATCH = 13


def error_msg_prn(code,infile,what):

    f = open("LOGUNIT","w")
    f.write("*********** RUN TIME ERROR BEGIN ***********")

    if code == ALLOC:
        f.write("Unable to allocate " + what + " in " + infile)
    elif code == NML_OPEN:
        f.write("Unable to open " + what + " in " + infile)
    elif code == NML_READ:
        f.write("Namelist mismatch in " + what + " opened by " + infile)
    elif code == DIM_MISMATCH:
        f.write("Dimension mismatch while reading " + what + " in " + infile)

    f.write("***********  RUN TIME ERROR END  ***********")
    print("BFM error (see logfile)")

    return

# EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  BFM - Biogeochemical Flux Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
