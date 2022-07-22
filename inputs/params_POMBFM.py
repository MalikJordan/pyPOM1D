# params_pombfm
h = 150.0           # Depth [m]
# dti = 100.0         # timestep [s]
dti = 360.
# dti = 1
alat = 45.0         # latitude [degrees]
idiagn = 1          # switch between prognostic (idiagn = 0) and diagnostic (idiagn = 1) mode
# idays = 3600        # length of run [days]
idays = 366
# idays = 100
smoth = 0.1         # parameter for hasselin filter
ihotst = 0          # switch for cold start (ihotst = 0) and hot start, ie reading restart (ihotst = 1)
kl1 = 2             # surface logarithmic layers distribution
kl2 = 150           # bottom logarithmic layers distribution
savef = 1           # output averaging and saving frequency [days]
nrt_o2o = 0.06      # relaxation velocity for oxygen [m/d]
nrt_n1p = 0.06      # relaxation velocity for phosphate [m/d]
nrt_n3n = 0.06      # relaxation velocity for nitrate [m/d]
nrt_n4n = 0.05      # relaxation velocity for ammonium [m/d]
nbct = 2            # flag for temperature boundary conditions
nbcs = 1            # flag for salinity boundary conditions
nbcbfm = 1          # flag for bfm boundary conditions
umol = 1e-06        # background diffusion
umolt = 1e-07       # background diffusion for temperature
umols = 1.3e-07     # background diffusion for salalinity
umolbfm = 0.0001    # background diffusion for bfm
ntp = 2             # flag for jerlov water type ( 1 = I, 2 = IA, 3 = IB, 4 = II, 5 = III)
trt = 0             # relaxation time for lateral temperature advection
srt = 1             # relaxation time for lateral salinity advection
upperh = 5.0        # depth where lateral advection starts
ssrt = 5.68         # relaxation time for surface salinity flux

