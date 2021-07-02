import numpy as np
from pom.pom_constants import vertical_layers

class DiffusionCoefficients:
    def __init__(self, tracers=np.zeros(vertical_layers), momentum=np.zeros(vertical_layers), kinetic_energy=np.zeros(vertical_layers)):
        self.tracers = tracers
        self.momentum = momentum
        self.kinetic_energy = kinetic_energy


class ForcingManagerCounters:
    def __init__(self, day_counter=0, day_interpolator=0, ratio_day=0, month_counter=0,
                 month_interpolator=0, ratio_month=0, timesteps_per_day=0, timesteps_per_month=0):
        self.day_counter = day_counter
        self.day_interpolator = day_interpolator
        self.ratio_day = ratio_day
        self.month_counter = month_counter
        self.month_interpolator = month_interpolator
        self.ratio_month = ratio_month
        self.timesteps_per_day = timesteps_per_day
        self.timesteps_per_month = timesteps_per_month


class LeapFrogTimeLevels:
    def __init__(self, current=np.zeros(vertical_layers), forward=np.zeros(vertical_layers), backward=np.zeros(vertical_layers)):
        self.current = current
        self.forward = forward
        self.backward = backward


class MonthlyForcingData:
    def __init__(self, sclim=np.zeros(vertical_layers), tclim=np.zeros(vertical_layers), wclim=np.zeros(vertical_layers),
                 weddy1=np.zeros(vertical_layers), weddy2=np.zeros(vertical_layers), ism=np.zeros(vertical_layers-1),
                 wsu=0, wsv=0, swrad=0, wtsurf=0, qcorr=0, NO3_s=0, NH4_s=0,
                 PO4_s=0, SIO4_s=0, O2_b=0, NO3_b=0, PO4_b=0, PON_b=0):
        self.sclim = sclim
        self.tclim = tclim
        self.wclim = wclim
        self.weddy1 = weddy1
        self.weddy2 = weddy2
        self.ism = ism
        self.wsu = wsu
        self.wsv = wsv
        self.swrad = swrad
        self.wtsurf = wtsurf
        self.qcorr = qcorr
        self.NO3_s = NO3_s
        self.NH4_s = NH4_s
        self.PO4_s = PO4_s
        self.SIO4_s = SIO4_s
        self.O2_b = O2_b
        self.NO3_b = NO3_b
        self.PO4_b = PO4_b
        self.PON_b = PON_b


class Stresses:
    def __init__(self, zonal=0, meridional=0):
        self.zonal = zonal
        self.meridional = meridional


class TemperatureSalinityData:
    def __init__(self, current=np.zeros(vertical_layers), forward=np.zeros(vertical_layers), backward=np.zeros(vertical_layers),
                 interpolated=np.zeros(vertical_layers), surface_value=0, surface_flux=0, bottom_flux=0,
                 lateral_advection=np.zeros(vertical_layers)):
        self.current = current
        self.forward = forward
        self.backward = backward
        self.interpolated = interpolated
        self.surface_value = surface_value
        self.surface_flux = surface_flux
        self.bottom_flux = bottom_flux
        self.lateral_advection = lateral_advection


class VerticalGridData:
    def __init__(self, length_scale=np.zeros(vertical_layers), vertical_coordinates=np.zeros(vertical_layers),
                 vertical_coordinates_staggered=np.zeros(vertical_layers), vertical_spacing=np.zeros(vertical_layers),
                 vertical_spacing_staggered=np.zeros(vertical_layers), vertical_spacing_reciprocal=np.zeros(vertical_layers)):
        self.length_scale = length_scale
        self.vertical_coordinates = vertical_coordinates
        self.vertical_coordinates_staggered = vertical_coordinates_staggered
        self.vertical_spacing = vertical_spacing
        self.vertical_spacing_staggered = vertical_spacing_staggered
        self.vertical_spacing_reciprocal = vertical_spacing_reciprocal
