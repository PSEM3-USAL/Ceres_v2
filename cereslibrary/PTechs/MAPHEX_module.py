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
# ==============================================================================================================================================================================================================================
def MAPHEX_module(F_ini,fc_ini,x_ini):
    import pandas as pd
    import numpy as np


    # IMPORT MODULES
    from global_parameters_module import UnitConv, MW, c_p_liq_sol, dH_vap_0, Tc, Tb, dH_f, dH_c, c_p_v_1, c_p_v_2, c_p_v_3, c_p_v_4, coef_vapor_pressure_1, coef_vapor_pressure_2, coef_vapor_pressure_3, CEI, price, nu_p, k_p, n_watson, epsilon, T_amb, P_ref, density
    from feedstock_input_module import elements_wet, elements_dry, nutrients, feedstock_parameters, elements_dry_comp
    from equipment_costs.MAPHEX_cost_module import MAPHEX_cost_module
    from economic_parameters_module import ec_param
    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================
    # SETS, PARAMETERS AND NODES DATA ACQUISITION
    nodes_matrix            = pd.read_csv('nodes/nodes_MAPHEX.csv', sep=",", header=None)
    nodes       = np.array(nodes_matrix[0].dropna())   
    #process_elements_matrix = pd.read_csv('process_elements/process_elements_FBR.csv', sep=",", header=0)
  
    
    
    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================
    # MANURE DATA ACQUISITION AND COMPOSITION SETTLEMENT
    #chemicals_comp   = np.array(process_elements_matrix["Component"][3:8].dropna())
    #chemicals_conc   = np.full((len(chemicals_comp)), 0)
    #chemicals        = dict(zip(chemicals_comp, chemicals_conc))
    
    #products_comp       = np.array(process_elements_matrix["Component"][0:3].dropna())
    #product_conc      = np.full((len(products_comp)), 0)
    #product            = dict(zip(products_comp, product_conc))
    
    
    
    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================
    # TOTAL ELEMENTS
    #total_elements = {**elements_wet,**chemicals,**product} # Merge the two dictionaries
    total_elements = elements_wet

    # VARIABLES DEFINITION (IN NESTED DICTIONARIES) (INITIALIZATION)
    nodes_list              = nodes.tolist()
    initialization_comp     = total_elements
    initialization_nan      = np.full((len(initialization_comp)), 0)
    
    fc = {key: dict(zip(initialization_comp,initialization_nan)) for key in nodes_list}
    
    x = {key: dict(zip(initialization_comp,initialization_nan)) for key in nodes_list}
    
    F = {key: np.nan for key in nodes_list}
    
    
    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================
    
    
    # MASS BALANCE
    for i in total_elements.keys():
#        x["Src1Mixer"][i]  = total_elements[i]/100
#        fc["Src1Mixer"][i] = x["Src1Mixer"][i]*F_ini
        x["Src_MAPHEX"][i]  = x_ini[i]
        fc["Src_MAPHEX"][i] = fc_ini[i]
    
    fc["MAPHEX_SolEff"]["N"]         = fc["Src_MAPHEX"]["N"]*0.93
    
    if (fc["Src_MAPHEX"]["N"]+fc["Src_MAPHEX"]["N-NH3"])*0.9-fc["Src_MAPHEX"]["N"] < fc["Src_MAPHEX"]["N-NH3"]:
        fc["MAPHEX_SolEff"]["N-NH3"]     = (fc["Src_MAPHEX"]["N"]+fc["Src_MAPHEX"]["N-NH3"])*0.9-fc["Src_MAPHEX"]["N"]
    else:
        fc["MAPHEX_SolEff"]["N-NH3"] = fc["Src_MAPHEX"]["N-NH3"]*0.9
    
    fc["MAPHEX_LiqEff"]["N"]         = fc["Src_MAPHEX"]["N"]-fc["MAPHEX_SolEff"]["N"]
    fc["MAPHEX_LiqEff"]["N-NH3"]     = fc["Src_MAPHEX"]["N-NH3"]-fc["MAPHEX_SolEff"]["N-NH3"]
    
    fc["MAPHEX_SolEff"]["P"]         = fc["Src_MAPHEX"]["P"]*0.93
    
    if (fc["Src_MAPHEX"]["P"]+fc["Src_MAPHEX"]["P-PO4"])*0.9-fc["Src_MAPHEX"]["P"] < fc["Src_MAPHEX"]["P-PO4"]:
        fc["MAPHEX_SolEff"]["P-PO4"]     = (fc["Src_MAPHEX"]["P"]+fc["Src_MAPHEX"]["P-PO4"])*0.9-fc["Src_MAPHEX"]["P"]
    else:
        fc["MAPHEX_SolEff"]["P-PO4"] = fc["Src_MAPHEX"]["P-PO4"]*0.9
    
    fc["MAPHEX_LiqEff"]["P"]         = fc["Src_MAPHEX"]["P"]-fc["MAPHEX_SolEff"]["P"]
    fc["MAPHEX_LiqEff"]["P-PO4"]     = fc["Src_MAPHEX"]["P-PO4"]-fc["MAPHEX_SolEff"]["P-PO4"]
    
    #It is considered that all the components are separated in the same amount than solids
    fc["MAPHEX_SolEff"]["Rest"]      = fc["Src_MAPHEX"]["Rest"]*0.93
    fc["MAPHEX_SolEff"]["C"]         = fc["Src_MAPHEX"]["C"]*0.93
    fc["MAPHEX_SolEff"]["Ca"]        = fc["Src_MAPHEX"]["Ca"]*0.93
    fc["MAPHEX_SolEff"]["K"]         = fc["Src_MAPHEX"]["K"]*0.93
    fc["MAPHEX_SolEff"]["Ca_ion"]    = fc["Src_MAPHEX"]["Ca_ion"]*0.93
    fc["MAPHEX_SolEff"]["K_ion"]     = fc["Src_MAPHEX"]["K_ion"]*0.93
    
    fc["MAPHEX_LiqEff"]["Rest"]      = fc["Src_MAPHEX"]["Rest"]*(1-0.93)
    fc["MAPHEX_LiqEff"]["C"]         = fc["Src_MAPHEX"]["C"]*(1-0.93)
    fc["MAPHEX_LiqEff"]["Ca"]        = fc["Src_MAPHEX"]["Ca"]*(1-0.93)
    fc["MAPHEX_LiqEff"]["K"]         = fc["Src_MAPHEX"]["K"]*(1-0.93)
    fc["MAPHEX_LiqEff"]["Ca_ion"]    = fc["Src_MAPHEX"]["Ca_ion"]*(1-0.93)
    fc["MAPHEX_LiqEff"]["K_ion"]     = fc["Src_MAPHEX"]["K_ion"]*(1-0.93)
    
    #Solids containt 75% moisture (https://doi.org/10.13031/aea.12632)
    fc["MAPHEX_SolEff"]["Wa"] = sum(fc["MAPHEX_LiqEff"][i] for i in elements_dry)*(0.75/0.25)
    fc["MAPHEX_LiqEff"]["Wa"] = fc["Src_MAPHEX"]["Wa"] - fc["MAPHEX_SolEff"]["Wa"]
    
    for i in nodes_list:
        F[i] = sum(fc[i][ii] for ii in total_elements)
        
    for i in nodes_list:
        if i!="Src_MAPHEX":
            for ii in total_elements.keys():
                x[i][ii] = fc[i][ii]/F[i]
    
    #Checks
    checks_store = ['OK', 'FAIL']
    checks_F = []
    checks_x = []
    
    #    for i in total_elements.keys():
    #        if abs(fc["Src1Mixer"][i] + fc["Src2Mixer"][i] + fc["Src3FBR"][i] - (fc["HydrocycloneSink1"][i] + fc["HydrocycloneSink2"][i] + fc["DryerSink3"][i] + fc["DryerSink4"][i])) <= 0.005:
    #            print("fc check", i, "OK")
    #        else:
    #            print("fc check", i, "FAIL")
    
    if abs(F["Src_MAPHEX"] - (F["MAPHEX_LiqEff"] + F["MAPHEX_SolEff"])) <= 0.005: #F["HydrocycloneSink1"] +
    #        print("F check OK")
        checks_F = checks_store[0]
    else:
    #        print("F check FAIL")
        checks_F = checks_store[1]
    
    for i in nodes_list:
        if abs(1-sum(x[i][ii] for ii in total_elements)) <= 0.005:
    #            print("x check", i, "OK")
            checks_x = checks_store[0]
        else:
    #            print("x check", i, "FAIL")
            checks_x = checks_store[1]
    print("\n\n")
        
    
    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================
    
    # ECONOMICS/DESIGN
   
    #    
    MAPHEX_results              = MAPHEX_cost_module(fc["Src_MAPHEX"]["P-PO4"])
    n_MAPHEX                    = MAPHEX_results['n_MAPHEX']
    MAPHEX_equipment_cost       = MAPHEX_results['MAPHEX_equipment_cost']
    MAPHEX_operating_cost       = MAPHEX_results['MAPHEX_operating_cost']
    
    equipment_cost              = MAPHEX_equipment_cost
    #    physical_plant_cost_2016    = 3.15*equipment_cost
    #    fixed_capital_cost_2016     = 1.4*physical_plant_cost_2016
    fixed_capital_cost_2016     = equipment_cost
        
    #operation_cost_2016 = (chemicals_cost_2016+0.3*fixed_capital_cost_2016+FBR_operating_cost_partial) #+1.5*labour_cost_2016
    operation_cost_2016_amortized = (equipment_cost/ec_param['plant_lifetime'])+MAPHEX_operating_cost #+1.5*labour_cost_2016
    operation_cost_2016_non_amortized = MAPHEX_operating_cost
    
    recovered_P = fc["MAPHEX_SolEff"]["P"]
    recovered_PO4 = fc["MAPHEX_SolEff"]["P-PO4"]
    released_P = fc["MAPHEX_LiqEff"]["P"]
    released_PO4 = fc["MAPHEX_LiqEff"]["P-PO4"]
    released_NH3 = (fc["MAPHEX_LiqEff"]["N-NH3"]) #kgN as NH4/s
    released_N = (fc["MAPHEX_LiqEff"]["N"]) #kg rest of N/s
    
    fraction_recoved_P = recovered_P/(fc["Src_MAPHEX"]["P"])
    fraction_recoved_PO4 = recovered_PO4/fc["Src_MAPHEX"]["P-PO4"]
    fraction_recoved_TP = (recovered_P+recovered_PO4)/(fc["Src_MAPHEX"]["P-PO4"]+fc["Src_MAPHEX"]["P"])
    fraction_released_P = released_P/(fc["Src_MAPHEX"]["P"])
    fraction_released_PO4 = released_PO4/(fc["Src_MAPHEX"]["P-PO4"])
    fraction_released_TP = (released_P+released_PO4)/(fc["Src_MAPHEX"]["P-PO4"]+fc["Src_MAPHEX"]["P"])
    PO4_conc_released = released_PO4/MW['P']*MW['PO4']*UnitConv['K_to_mili']/(F["MAPHEX_LiqEff"]/feedstock_parameters['Wa_density'])
    Eutrophication_potential = (released_N*0.42+released_NH3/MW['N']*MW['NH4']*0.33+released_PO4/MW['P']*MW['PO4']*1+released_P*3.06)*365*3600*24/1000 #ton phosphate equivalent / year , potency factors from 'The sustainability metrics. Sustainable Development Progess Metrics recommended for use in the Process Idustries' IChenE
    
    return {'tech':'MAPHEX',
            'fixed_capital_cost_2016':fixed_capital_cost_2016,
            'n_MAPHEX':n_MAPHEX,
            'MAPHEX_size':'Unique',
            'investment_cost':fixed_capital_cost_2016,
            'fc':fc,
            'x':x,
            'F':F,
            'check_F':checks_F,
            'checks_x':checks_x,
            'equipment_cost':MAPHEX_equipment_cost,
            'operation_cost_2016_amortized':operation_cost_2016_amortized,
            'operation_cost_2016_non_amortized':operation_cost_2016_non_amortized,
            'recovered_P':recovered_P,
            'recovered_PO4':recovered_PO4,
            'released_P':released_P,
            'released_PO4':released_PO4,
            'fraction_recoved_P':fraction_recoved_P,
            'fraction_recoved_PO4':fraction_recoved_PO4,
            'fraction_released_P':fraction_released_P,
            'fraction_released_PO4':fraction_released_PO4,
            'fraction_recoved_TP':fraction_recoved_TP,
            'fraction_released_TP':fraction_released_TP,
            'PO4_conc_released':PO4_conc_released,          
            'Eutrophication_potential':Eutrophication_potential,  
            }
    
    #print("fc: \n"),  pprint.pprint(fc), print("\n\n")
    #print("x: \n"),  pprint.pprint(x), print("\n\n")
    #print("F: \n"),  pprint.pprint(F), print("\n\n")
    #print("Cost_prec_tank_2016: ", Cost_prec_tank_2016)
    #print("Centrifuge_cost_2016: ", Centrifuge_cost_2016)
    #print("Chemicals_cost: ", Chemicals_cost)
    #print("Labour_cost: ", Labour_cost)
    #print("Operational_cost: ", Operational_cost)
    #print("Nutrients_benefits: ", Nutrients_benefits)
    #print("Benefits: ", Benefits)
    #print("\n\n**************************************************************** \n\n")



# ===============================================================================================================================================================================================================================
# ###############################################################################################################################################################################################################################
# ===============================================================================================================================================================================================================================
    
    
# GRAPH

#G=pgv.AGraph(directed=True,)
#G.edge_attr.update(len='2.0',color='black')
##G.node_attr['style']='filled'
##nodelist=['Src1','Filter','Sink1', 'Sink2']
##G.add_nodes_from(nodelist, color='white')
#G.add_node('Src1', style='filled', color='darkgoldenrod2',shape='diamond')
#G.add_node('Src2', style='filled', color='darkgoldenrod2',shape='diamond')
#G.add_node('PrecTank', style='filled', color='darkgoldenrod2',shape='diamond')
#G.add_node('Centrif', style='filled', color='darkgoldenrod2',shape='diamond')
#G.add_node('Sink1', style='filled', color='darkgoldenrod2',shape='diamond')
#G.add_node('Sink2', style='filled', color='darkgoldenrod2',shape='diamond')
#G.add_edge('Src1','PrecTank')
#G.add_edge('Src2','PrecTank')
#G.add_edge('PrecTank','Centrif')
#G.add_edge('Centrif','Sink1')
#G.add_edge('Centrif','Sink2')
#
#G.graph_attr['label']='Centrifugation'
##G.node_attr['shape']='circle'
#
##G.string()
#G.layout(prog='dot') 
#G.draw('Centrifugation.png')
#print("Centrifugation.png")



