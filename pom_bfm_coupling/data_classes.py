from calendar import WEDNESDAY
import numpy as np
from pom.constants import vertical_layers
from inputs.params_POMBFM import idays
from bfm.constants import num_d3_box_states

num_boxes = vertical_layers - 1

class BfmStateVariableData:
    def __init__(self, current=np.zeros(vertical_layers), forward=np.zeros(vertical_layers), backward=np.zeros(vertical_layers),
                 surface_value=0, surface_flux=0, bottom_flux=0):
        self.current = current
        self.forward = forward
        self.backward = backward
        self.surface_value = surface_value
        self.surface_flux = surface_flux
        self.bottom_flux = bottom_flux


class BfmPhysicalVariableData:
    # def __init__(self, temperature=np.zeros(vertical_layers), salinity=np.zeros(vertical_layers), density=np.zeros(vertical_layers), 
    #              suspended_matter=np.zeros(vertical_layers-1), depth=np.zeros(vertical_layers), irradiance=np.zeros(vertical_layers-1), 
    #              vertical_extinction=np.zeros(vertical_layers-1), wind=0, wgen=np.zeros(vertical_layers), weddy=np.zeros(vertical_layers)):
    def __init__(self, temperature=np.zeros(vertical_layers), salinity=np.zeros(vertical_layers), density=np.zeros(vertical_layers), 
                 suspended_matter=np.zeros(num_boxes), depth=np.zeros(vertical_layers), irradiance=np.zeros((num_boxes,4)), 
                 vertical_extinction=np.zeros((num_boxes,4)), wind=0, wgen=np.zeros(vertical_layers), weddy=np.zeros(vertical_layers)):
        self.temperature = temperature
        self.salinity = salinity
        self.density = density
        self.suspended_matter = suspended_matter
        self.depth = depth
        self.irradiance = irradiance
        self.vertical_extinction = vertical_extinction
        self.wind = wind
        self.wgen = wgen
        self.weddy = weddy


class NutrientData:
    def __init__(self,NO3surf,NH4surf,PO4surf,SIO4surf,O2bott,NO3bott,PO4bott,PONbott_grad):
        self.NO3surf = NO3surf
        self.NH4surf = NH4surf
        self.PO4surf = PO4surf
        self.SIO4surf = SIO4surf
        self.O2bott = O2bott
        self.NO3bott = NO3bott
        self.PO4bott = PO4bott
        self.PONbott_grad = PONbott_grad


class AverageData:
    def __init__(self,count=0,day=0,single_day_ave=np.zeros((num_boxes,num_d3_box_states)),daily_ave=np.zeros((num_boxes,num_d3_box_states,idays))):
        self.count = count
        self.day = day
        self.single_day_ave = single_day_ave
        self.daily_ave = daily_ave
