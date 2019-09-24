#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 09:15:46 2018

@author: emh
"""

#Ceres_Filtration_v1.
#
#Tool for selection of nutrients recovery technologies.
#
#Edgar Martín Hernández.
#
#Cincinnati 2019.


# ===============================================================================================================================================================================================================================
# ###############################################################################################################################################################################################################################
# ===============================================================================================================================================================================================================================
# IMPORT LIBRARIES MODULE
import numpy as np
import pandas as pd
import scipy.optimize as opt
from scipy.integrate import odeint, solve_bvp
import matplotlib.pyplot as plt

pd.options.display.max_columns = 50

#F_tot=1 #kg/s
#fc_P_PO4 = 0.00031919000000000001
#fc_N_NH4 = 0.0022656
#fc_Ca_ion = 0.00018018
#fc_MgCl2 = 0.0019606709090322582
#fc_struvite = 0.0019285752044940962
#NH4_molar = 0.04124722 
#Mg_molar = 0.002562
#PO4_molar = 0.00043357894736842102

def MULTIFORM_cost_module(fc_P_PO4_feed):
    
    # IMPORT MODULES
    from global_parameters_module import UnitConv, MW, c_p_liq_sol, dH_vap_0, Tc, Tb, dH_f, dH_c, c_p_v_1, c_p_v_2, c_p_v_3, c_p_v_4, coef_vapor_pressure_1, coef_vapor_pressure_2, coef_vapor_pressure_3, CEI, price, nu_p, k_p, n_watson, epsilon, T_amb, P_ref, density, latent_heat_evap, nat_gas_heat_value
    #from feedstock_input_module import elements_wet, elements_dry, nutrients, feedstock_parameters, elements_dry_comp

    ec_param_matrix = pd.read_csv('economic_datasheets/operation_parameters.csv', sep=",", header=0)
    #ec_param_matrix = pd.read_csv('cereslibrary/techmodels/economic_datasheets/operation_parameters.csv', sep=",", header=0)
    ec_param        = dict(zip(np.array(ec_param_matrix["Item"].dropna()), np.array(ec_param_matrix["Value"].dropna())))

    # ====================================================
    # ####################################################
    # ====================================================
    #INPUT CONCENTRATIONS
    fc_P_PO4_feed_day = fc_P_PO4_feed*3600*24 # to kgPO4/day (feed)


    # ====================================================
    # ####################################################
    # ====================================================
    #SIZE DETERMINATION
    MULTIFORM_size = 38.5 # kgP-PO4feed/dayunit

        
    #NUMBER OF UNITS
    n_MULTIFORM = np.ceil(fc_P_PO4_feed_day/MULTIFORM_size)

    # ====================================================
    # ####################################################
    # ====================================================
    #Unit cost(OSTARA) includes all the system
    MULTIFORM_cost = 1045000 #USD/unit

    MULTIFORM_equipment_cost = MULTIFORM_cost*n_MULTIFORM #Total installation

    # ====================================================
    # ####################################################
    # ====================================================
    #OPERATION COST (except chemicals)
    MULTIFORM_OperatingCost = 15.419 # USD/kgP-PO4feed

    MULTIFORM_operating_cost_partial = MULTIFORM_OperatingCost*fc_P_PO4_feed_day*365

    return {'MULTIFORM_size':MULTIFORM_size, 'n_MULTIFORM':n_MULTIFORM, 'MULTIFORM_equipment_cost':MULTIFORM_equipment_cost, 'MULTIFORM_operating_cost_partial':MULTIFORM_operating_cost_partial}



# plot results
#plt.figure()
#plt.plot(H_inverse, L)
#plt.xlabel('Height (m)')
#plt.ylabel('L(H) (m)')
#plt.show()

#plt.figure()
#plt.plot(H, C)
#plt.xlabel('Height (m)')
#plt.ylabel('C(H) (mol/L)')
#plt.show()

#plt.figure()
#plt.plot(H, result.y[0])
#plt.xlabel('Height (m)')
#plt.ylabel('L(H) (m)')
#plt.show()
