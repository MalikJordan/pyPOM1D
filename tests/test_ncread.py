import netCDF4 as nc

f = '/Users/malikjordan/Desktop/BFM17_POM1D/BFM17_POM1D_VrsFnl/bfm_run/bfm17_pom1d/bfm17_pom1d.nc'

data = nc.Dataset(f)
oxy = data['O2o']
x=1