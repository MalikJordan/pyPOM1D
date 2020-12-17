# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# # ROUTINE: Service
#
# DESCRIPTION
#
#   List of Fortran parameters.
#   ONE-DIMENSIONAL BFM-POM SYSTEM
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from BFM17_POM1D_VrsFnl.src.BFM.General.ModuleGlobalMem import RLEN
from BFM17_POM1D_VrsFnl.src.pom.phys.POMModule import KB
from decimal import *
import numpy as np

getcontext().prec = 12

# SURFACE NUTRIENTS
PO4SURF, NO3SURF, NH4SURF, SIO4SURF = Decimal()
PO4BOTT, NO3BOTT, O2BOTT, PONBOTTgrad = Decimal()

# SUSPENDED INORGANIC MATTER PROFILE
ISM = np.empty(KB-1, dtype=float)
WGEN, WEDDY = np.empty(KB, dtype=float)

# FREQUENCY OF AVERAGING FOR OUTPUTi (IN DAYS)
savef, nitend = Decimal()
deltat = Decimal()

# THESE ARE THE PATHWAYS FOR THE IC AND FORCING FILES (READ TROUGH NML)
wind_input, \
surfaceS_input, \
radiance_input, \
ism_input, \
Sal_input, \
Temp_input, \
W_input, \
Weddy_input1, \
Weddy_input2, \
Sprofile_input, \
Tprofile_input, \
heat_input, \
surfNut_input, \
bottNut_input, \
read_restart = np.chararray((1,200))

# THESE ARE THE PATHWAYS FOR THE BFM17 IC (READ TROUGH NML)
phyto_input, \
zoop_input, \
poc_input, \
doc_input, \
phos_input, \
nit_input, \
am_input, \
oxy_input = np.chararray((1,200))

# EOC
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MODEL  POM - Princeton Ocean Model
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
