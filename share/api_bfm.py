from cppdefs import *
import numpy as np
from pom.modules_old import RLEN, zero, bfm_lwp, LOGUNIT, no_d3_box_states

dim = int()
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: API_BFM
#
# DESCRIPTION:  API for the BFM.
#               Storage of variables and diagnostics
#               To be used in all the coupled applications except
#               GOTM, where it actually originated from.
#               The GOTM module netcdfout is needed
#               Appropriate functions are already available in the GOTM structure
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# PUBLIC DATA MEMBERS
bio_calc = bool()
bioshade = bool()
feedback = bool()
bfm_rstctl = bool()

bio_setup = 1
bfm_init  = 0
surface_flux_method = -1
bottom_flux_method = -1
n_surface_fluxes = -1
calc_init_bennut_states = int()

out_dir = ''
out_fname = ''
out_title = ''

out_units = int()
out_delta = int()
out_secs = int()
save_delta = int()
time_delta = int()

rst_fname = ''
rst_fname_3d = ''

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   PARAMETERS FOR MASSIVE PARALLEL COMPUTATION
#   THE FOLLOWING ARE THE DEFAULT VALUES IF THE MACRO BFM_PARALLEL IS NOT DEFINED
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
parallel = False
parallel_log = False
parallel_rank = 0
# str = ' '*4

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   DIMENSION LENGTHS FOR OUTPUT
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
lon_len = int()
lat_len = int()
depth_len = int()
ocepoint_len = int()
surfpoint_len = int()
botpoint_len = int()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   BFM VARIABLE INFORMATION FOR OUTPUT
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
var_ids = np.empty(dim,dtype=int)
var_ave = np.empty(dim,dtype=bool)
ave_count = float()
ave_ct1 = False
D3ave = np.empty((dim,dim),dtype=float)
D2ave = np.empty((dim,dim),dtype=float)

try:
    import INCLUDE_SEAICE
    INCLUDE_SEAICE = True
except FileNotFoundError:
    INCLUDE_SEAICE = False
if INCLUDE_SEAICE:
    D2ave_ice = np.empty((dim,dim),dtype=float)

try:
    import INCLUDE_BEN
    INCLUDE_BEN = True
except FileNotFoundError:
    INCLUDE_BEN = False
if INCLUDE_BEN:
    D2ave_ben = np.empty((dim,dim),dtype=float)

var_names = np.empty(dim,dtype=str)
var_units = np.empty(dim,dtype=str)
var_long = np.empty(dim,dtype=str)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   INCICES OF THE VARIOUS OUTPUT VARIABLES
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
stStart = 0
stEnd = 0

stPelStateS = 0
stPelDiagS = 0
stPelFluxS = 0

stPelDiag2dS = 0

stPelSurS = 0
stPelBotS = 0
stPelRivS = 0

stPelStateE = 0
stPelDiagE = 0
stPelFluxE = 0

stPelDiag2dE = 0

stPelSurE = 0
stPelBotE = 0
stPelRivE = 0

stPelStart = 0
stPelEnd = 0

try:
    import INCLUDE_SEAICE
    INCLUDE_SEAICE = True
except FileNotFoundError:
    INCLUDE_SEAICE = False
if INCLUDE_SEAICE:
    stIceStateS = 0
    stIceDiag2dS = 0
    stIceFlux2dS = 0

    stIceStateE = 0
    stIceDiag2dE = 0
    stIceFlux2dE = 0

    stIceStart = 0
    stIceEnd = 0

try:
    import INCLUDE_BEN
    INCLUDE_BEN = True
except FileNotFoundError:
    INCLUDE_BEN = False
if INCLUDE_BEN:
    stBenStateS = 0
    stBenDiag2dS = 0
    stBenFlux2dS = 0

    stBenStateE = 0
    stBenDiag2dE = 0
    stBenFlux2dE = 0

    stBenStart = 0
    stBenEnd = 0

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   ADDITIONAL OUTPUT VARIABLES
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
c1dim = np.empty(dim,dtype=float)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   BFM variable information for data input
#   integer init: select the initialization
#                 0 = homogeneous
#                 1 = analytical
#                 2 = from file
#   options for init==1
#   real anv1: value in the surface layer
#   real anz1: depth of the surface layer
#   real anv2: value in the bottom layer
#   real anz2: depth of the bottom layer
#   options for init==2
#   char filename: name of the input file
#   char  varname: name of the var in input file
#   Options currently used when coupled with NEMO
#   logical obc: variable has open boundary data
#   logical sbc: variable has surface boundary data
#   logical cbc: variable has coastal boundary data
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
class InputInfo:

    def __init__(self,init,filename,varname,anz1,anv1,anz2,anv2,obc,sbc,cbc):
        self.init = init
        self.filename = filename
        self.varname = varname
        self.anz1 = anz1
        self.anv1 = anv1
        self.anz2 = anz2
        self.anv2 = anv2
        self.obc = obc
        self.sbc = sbc
        self.cbc = cbc


InitVar = InputInfo

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   ADDITIONAL 1D ARRAYS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# INDICES OF BOTTOM AND SURFACE POINTS
BOTindices = np.empty(dim,dtype=int)
SRFindices = np.empty(dim,dtype=int)

# REAL MASK OF RIVER POINTS AT SURFACE
RIVmask = np.empty(dim,dtype=float)

# TOTAL AMOUNT FOR EACH VARIABLE
D3STATE_tot = np.empty(dim,dtype=float)

try:
    import INCLUDE_SEAICE
    INCLUDE_SEAICE = True
except FileNotFoundError:
    INCLUDE_SEAICE = False
if INCLUDE_SEAICE:
    D2STATE_ICE_tot = np.empty(dim,dtype=float)

try:
    import INCLUDE_BEN
    INCLUDE_BEN = True
except FileNotFoundError:
    INCLUDE_BEN = False
if INCLUDE_BEN:
    D2STATE_BEN_tot = np.empty(dim,dtype=float)


try:
    import BFM_NEMO
    BFM_NEMO = True
except FileNotFoundError:
    BFM_NEMO = False
if BFM_NEMO:
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   ADDITIONAL 3D ARAYS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ZEROS = np.empty((dim,dim,dim),dtype=float)
    # 3D BOOLEAN LAND-SEA MASK
    SEAmask = np.empty((dim,dim,dim),dtype=bool)
    # 3D BOOLEAN SEA-BOTTOM MASK
    BOTmask = np.empty((dim,dim,dim),dtype=bool)
    # 3D BOOLEAN MASK OF THE SURFACE POINTS
    SRFmask = np.empty((dim,dim,dim),dtype=bool)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   3D INDICES OF THE WET POINTS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    iwet = np.empty(dim,dtype=int)
    jwet = np.empty(dim,dtype=int)
    kwet = np.empty(dim,dtype=int)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   ADDITIONAL INTEGRATION ARRAYS FOR LEAPFROG SCHEME
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    D3STATEB = np.empty((dim,dim),dtype=float)
    try:
        import INCLUDE_SEAICE
        INCLUDE_SEAICE = True
    except FileNotFoundError:
        INCLUDE_SEAICE = False
    if INCLUDE_SEAICE:
        D2STATEB_ICE = np.empty((dim,dim),dtype=float)

    try:
        import INCLUDE_BEN
        INCLUDE_BEN = True
    except FileNotFoundError:
        INCLUDE_BEN = False
    if INCLUDE_BEN:
        D2STATEB_BEN = np.empty((dim, dim), dtype=float)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=











