import numpy as np
from pom.constants import num_boxes
from bfm.bfm50.Functions.other_functions import eTq_vector
from inputs import params_POMBFM

def pel_chem_eqns(pel_chem_parameters, environmental_parameters, constant_parameters, del_z, temper, conc, flPTn6r):
    """ calculates the non-living equations for DOM, POM, and nutrients """
    
    # State variables
    o2o = conc[:,0]              # Dissolved oxygen (mg O_2 m^-3)
    n3n = conc[:,2]              # Nitrate (mmol N m^-3)
    n4n = conc[:,3]              # Ammonium (mmol N m^-3)
    n6r = conc[:,6]              # Reduction equivalents (mmol S m^-3)
    r1c = conc[:,39]             # Labile dissolved organic carbon (mg C m^-3)
    r1n = conc[:,40]             # Labile dissolved organic nitrogen (mmol N m^-3)
    r1p = conc[:,41]             # Labile dissolved organic phosphate (mmol P m^-3)
    r6c = conc[:,44]             # Particulate organic carbon (mg C m^-3)
    r6n = conc[:,45]             # Particulate organic nitrogen (mmol N m^-3)
    r6p = conc[:,46]             # Particulate organic phosphate (mmol P m^-3)
    r6s = conc[:,47]             # Particulate organic silicate (mmol Si m^-3)
    
    # Regulating factors
    eo = np.zeros(num_boxes)
    for i in range(0,num_boxes):
        eo[i] = max(constant_parameters["p_small"], o2o[i])/(max(constant_parameters["p_small"], o2o[i])+ pel_chem_parameters["h_o"])
    er = n6r/(n6r + pel_chem_parameters["h_r"])
    
    # Temperature regulating factors
    fTn = eTq_vector(temper, environmental_parameters["basetemp"], environmental_parameters["q10n"])
    fTr6 = eTq_vector(temper, environmental_parameters["basetemp"], environmental_parameters["q10n5"])
    
    # Nitrification in the water  [mmol N m^-3 s^-1]   
    dn4ndt_nit_n3n = np.zeros(num_boxes)
    for i in range(0,num_boxes):
        dn4ndt_nit_n3n[i] = max(0.0, pel_chem_parameters["lambda_n4nit"]*n4n[i]*fTn[i]*eo[i])

    # Denitrification flux [mmol N m^-3 s^-1] from PelChem.F90 line 134
    rPAo = flPTn6r/constant_parameters["omega_r"]
    dn3ndt_denit = np.zeros(num_boxes)
    for i in range(0,num_boxes):
        dn3ndt_denit[i] = max(0.0, pel_chem_parameters["lambda_N3denit"]*fTn[i]*er[i]*rPAo[i]/pel_chem_parameters["m_o"]*n3n[i])
    
    # Reoxidation of reduction equivalents [mmol S m^-3 s^-1]
    dn6rdt_reox = pel_chem_parameters["lambda_n6reox"]*eo*n6r
    
    # Dissolution of biogenic silicate [mmol Si m^-3 s^-1]
    dr6sdt_rmn_n5s = pel_chem_parameters["lambda_srmn"]*fTr6*r6s
    # dr6sdt_rmn_n5s = np.zeros(num_boxes)

    # Constant organic matter remineralization from PelChem.F90
    # bacprof = (1/np.log(params_POMBFM.h))*(np.log(params_POMBFM.h) - del_z)
    p_sR6O3 = 0.1
    p_sR1O3 = 0.05
    dr6cdt_remin_o3c = p_sR6O3*r6c
    dr1cdt_remin_o3c = p_sR1O3*r1c
    dr6cdt_remin_o2o = dr6cdt_remin_o3c
    dr1cdt_remin_o2o = dr1cdt_remin_o3c

    sR6N1 = 0.1
    sR1N1 = 0.05
    dr6pdt_remin_n1p = sR6N1*r6p
    # dr1pdt_remin_n1p = sR1N1*bacprof*r1p
    dr1pdt_remin_n1p = sR1N1*r1p

    p_sR6N4 = 0.1
    p_sR1N4 = 0.05
    dr6ndt_remin_n4n = p_sR6N4*r6n
    dr1ndt_remin_n4n = p_sR1N4*r1n
    
    return (dn4ndt_nit_n3n, dn3ndt_denit, dn6rdt_reox, dr6sdt_rmn_n5s, dr6cdt_remin_o3c, dr1cdt_remin_o3c, dr6cdt_remin_o2o, dr1cdt_remin_o2o, 
            dr6pdt_remin_n1p, dr1pdt_remin_n1p, dr6ndt_remin_n4n, dr1ndt_remin_n4n)
