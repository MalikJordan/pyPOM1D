import numpy as np

def insw_vector(parameter):
    switch = 0.0
    if parameter > 0.0:
        switch = 1.0

    return switch


def eTq_vector(temp, basetemp, q10):
    eTq = np.exp(np.log(q10)*(temp-basetemp)/basetemp)
    
    return eTq

def get_concentration_ratio(numerator, denominator, p_small):
    concentration_ratio = 0.0
    if numerator > 0:
        concentration_ratio = numerator/(denominator + p_small)
    
    return concentration_ratio