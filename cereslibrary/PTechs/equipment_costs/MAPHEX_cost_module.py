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

def MAPHEX_cost_module(fc_P_PO4_feed):
    import pandas as pd
    import numpy as np
    
    # IMPORT MODULES
    from global_parameters_module import UnitConv, MW, c_p_liq_sol, dH_vap_0, Tc, Tb, dH_f, dH_c, c_p_v_1, c_p_v_2, c_p_v_3, c_p_v_4, coef_vapor_pressure_1, coef_vapor_pressure_2, coef_vapor_pressure_3, CEI, price, nu_p, k_p, n_watson, epsilon, T_amb, P_ref, density, latent_heat_evap, nat_gas_heat_value
    from feedstock_input_module import elements_wet, elements_dry, nutrients, feedstock_parameters, elements_dry_comp

    ec_param_matrix = pd.read_csv('economic_datasheets/operation_parameters.csv', sep=",", header=0)
    #ec_param_matrix = pd.read_csv('cereslibrary/techmodels/economic_datasheets/operation_parameters.csv', sep=",", header=0)
    ec_param        = dict(zip(np.array(ec_param_matrix["Item"].dropna()), np.array(ec_param_matrix["Value"].dropna())))

    # ====================================================
    # ####################################################
    # ====================================================
    #NUMBER OF UNITS
    daily_PO4_feed = fc_P_PO4_feed*3600*24
    n_MAPHEX = np.ceil(daily_PO4_feed/18.54)

    # ====================================================
    # ####################################################
    # ====================================================
    #Unit cost(MAPHEX) includes all the system
    MAPHEX_cost = 291000 #USD
    MAPHEX_equipment_cost = MAPHEX_cost*n_MAPHEX #Total installation

    # ====================================================
    # ####################################################
    # ====================================================
    #OPERATION COST 
    MAPHEX_operating_cost = 110.8*daily_PO4_feed*365#USD/year

    return {'n_MAPHEX':n_MAPHEX, 'MAPHEX_equipment_cost':MAPHEX_equipment_cost, 'MAPHEX_operating_cost':MAPHEX_operating_cost}
