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
#def FBR_module(F_ini):

# IMPORT LIBRARIES MODULE
#import numpy as np
#import pandas as pd

#pd.options.display.max_columns = 50
##F_ini=0.12683916793505834
def PROC_module(F_ini,fc_ini,x_ini):
    import pandas as pd
    import numpy as np
    
    # IMPORT MODULES
    from global_parameters_module import UnitConv, MW, c_p_liq_sol, dH_vap_0, Tc, Tb, dH_f, dH_c, c_p_v_1, c_p_v_2, c_p_v_3, c_p_v_4, coef_vapor_pressure_1, coef_vapor_pressure_2, coef_vapor_pressure_3, CEI, price, nu_p, k_p, n_watson, epsilon, T_amb, P_ref, density
    from feedstock_input_module import elements_wet, elements_dry, nutrients, feedstock_parameters, elements_dry_comp
    from equipment_costs.CSTR_cost_module import CSTR_investment_cost
    from equipment_costs.BeltFilt_design_cost import BeltFilt_design_cost
    from equipment_costs.ConveyorDryer_design_cost import ConveyorDryer_design_cost
    from equipment_costs.vessel_design_cost import vessel_design_cost
    from economic_parameters_module import ec_param
    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================
    # SETS, PARAMETERS AND NODES DATA ACQUISITION
    nodes_matrix            = pd.read_csv('nodes/nodes_PROC.csv', sep=",", header=None)
    #nodes_matrix            = pd.read_csv('cereslibrary/techmodels/nodes/nodes_FBR.csv', sep=",", header=None)
    nodes       = np.array(nodes_matrix[0].dropna())   
    process_elements_matrix = pd.read_csv('process_elements/process_elements_PROC.csv', sep=",", header=0)
    #process_elements_matrix = pd.read_csv('cereslibrary/techmodels/process_elements/process_elements_FBR.csv', sep=",", header=0)
    
    #Specefic technology parameters
    CSH_P_ratio          = (9+12)/2 # #kg/kg_Precovered doi:10.1016/j.scitotenv.2016.07.019
    yield_PO4            = 0.6 # P-RoC technology - field of application
    humidity_CSTR        = 1-0.2693 # %mass Struvite humidity in the moment of leave the reactor (%mass) 10.2166/wst.2014.236
    humidity_BeltFilt        = 0.5 # %mass Struvite humidity out belt filter (%mass) Walas 2012
    
    
    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================
    
    
    # MANURE DATA ACQUISITION AND COMPOSITION SETTLEMENT
    chemicals_comp   = np.array(process_elements_matrix["Component"][3:8].dropna())
    chemicals_conc   = np.full((len(chemicals_comp)), 0)
    chemicals        = dict(zip(chemicals_comp, chemicals_conc))
    
    products_comp       = np.array(process_elements_matrix["Component"][0:3].dropna())
    product_conc      = np.full((len(products_comp)), 0)
    product            = dict(zip(products_comp, product_conc))
    
    
    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================
    # TOTAL ELEMENTS
    total_elements = {**elements_wet,**chemicals,**product} # Merge the two dictionaries
    
    
    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================
    
    
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
#        x["Src1_PROC"][i]  = total_elements[i]/100
#        fc["Src1_PROC"][i] = x["Src1_MULTIFORM"][i]*F_ini
        x["Src1_PROC"][i]  = x_ini[i]
        fc["Src1_PROC"][i] = fc_ini[i]
        

    fc["Src2_PROC"]["CSH"]         = (fc["Src1_PROC"]["P-PO4"]*yield_PO4)*CSH_P_ratio #UnitOperatingCost(USD/kgP-PO4feed day)
    
    
    
    #Total suspended solids penalization in precipitates formation
    #DM_Src1_PROC_percentage      = 100*(fc["Src1_PROC"]["P"]+fc["Src1_PROC"]["N"]+fc["Src1_PROC"]["C"]+fc["Src1_PROC"]["Rest"]+fc["Src1_PROC"]["Ca"]+fc["Src1_PROC"]["K"])/sum(fc["Src1_PROC"][ii] for ii in total_elements)
    #TSS_Src1_PROC_percentage     = 0.078*DM_Src1_PROC_percentage + 0.051
    penalization_index = 1
    #if TSS_Src1_MULTIFORM_percentage <= 0.215:
        #penalization_index = 1
    #else:
        #penalization_index = -7.335E-01*TSS_Src1_MULTIFORM_percentage + 1.158
        
    fc["PROC_Sedimentator"]["P-PO4"] = fc["Src1_PROC"]["P-PO4"]*(1-yield_PO4)
    fc["PROC_Sedimentator"]["CSH-P"] = fc["Src2_PROC"]["CSH"] + fc["Src1_PROC"]["P-PO4"]*(yield_PO4)
    
    fc["PROC_Sedimentator"]["Rest"]      = fc["Src1_PROC"]["Rest"]
    fc["PROC_Sedimentator"]["C"]         = fc["Src1_PROC"]["C"]
    fc["PROC_Sedimentator"]["P"]         = fc["Src1_PROC"]["P"]
    fc["PROC_Sedimentator"]["N"]         = fc["Src1_PROC"]["N"]
    fc["PROC_Sedimentator"]["Ca"]        = fc["Src1_PROC"]["Ca"]
    fc["PROC_Sedimentator"]["K"]         = fc["Src1_PROC"]["K"]
    fc["PROC_Sedimentator"]["N-NH3"]     = fc["Src1_PROC"]["N-NH3"]
    fc["PROC_Sedimentator"]["Ca_ion"]    = fc["Src1_PROC"]["Ca_ion"]
    fc["PROC_Sedimentator"]["K_ion"]     = fc["Src1_PROC"]["K_ion"]
    fc["PROC_Sedimentator"]["Wa"]        = fc["Src1_PROC"]["Wa"]
    
    fc["Sedimentator_BeltFilter"]["CSH-P"]     = fc["PROC_Sedimentator"]["CSH-P"]
    fc["Sedimentator_BeltFilter"]["C"]         = fc["PROC_Sedimentator"]["C"]
    fc["Sedimentator_BeltFilter"]["P"]         = fc["PROC_Sedimentator"]["P"]
    fc["Sedimentator_BeltFilter"]["N"]         = fc["PROC_Sedimentator"]["N"]
    fc["Sedimentator_BeltFilter"]["Ca"]        = fc["PROC_Sedimentator"]["Ca"]
    fc["Sedimentator_BeltFilter"]["K"]         = fc["PROC_Sedimentator"]["K"]
    
    fc["Sedimentator_BeltFilter"]["Wa"]        = (fc["Sedimentator_BeltFilter"]["CSH-P"]+fc["Sedimentator_BeltFilter"]["C"]+fc["Sedimentator_BeltFilter"]["P"]+fc["Sedimentator_BeltFilter"]["N"]+fc["Sedimentator_BeltFilter"]["Ca"]+fc["Sedimentator_BeltFilter"]["K"])*humidity_CSTR/(1-humidity_CSTR)   #Computing the amount of water contained in solids
    
    Relative_water = fc["Sedimentator_BeltFilter"]["Wa"]/fc["PROC_Sedimentator"]["Wa"] #Ratio water solids / water liquid phase
    
    fc["Sedimentator_BeltFilter"]["P-PO4"]     = fc["Src1_PROC"]["P-PO4"]*Relative_water
    fc["Sedimentator_BeltFilter"]["N-NH3"]     = fc["Src1_PROC"]["N-NH3"]*Relative_water
    fc["PROC_Sedimentator"]["Ca_ion"]    = fc["Src1_PROC"]["Ca_ion"]*Relative_water
    fc["Sedimentator_BeltFilter"]["K_ion"]     = fc["Src1_PROC"]["K_ion"]*Relative_water
    
    
    fc["Sedimentator_LiqEff"]["CSH-P"]     = fc["PROC_Sedimentator"]["CSH-P"] - fc["Sedimentator_BeltFilter"]["CSH-P"]
    fc["Sedimentator_LiqEff"]["Rest"]      = fc["PROC_Sedimentator"]["Rest"] - fc["Sedimentator_BeltFilter"]["Rest"]
    fc["Sedimentator_LiqEff"]["C"]         = fc["PROC_Sedimentator"]["C"] - fc["Sedimentator_BeltFilter"]["C"]
    fc["Sedimentator_LiqEff"]["P"]         = fc["PROC_Sedimentator"]["P"] - fc["Sedimentator_BeltFilter"]["P"]
    fc["Sedimentator_LiqEff"]["N"]         = fc["PROC_Sedimentator"]["N"] - fc["Sedimentator_BeltFilter"]["N"]
    fc["Sedimentator_LiqEff"]["Ca"]        = fc["PROC_Sedimentator"]["Ca"] - fc["Sedimentator_BeltFilter"]["Ca"]
    fc["Sedimentator_LiqEff"]["K"]         = fc["PROC_Sedimentator"]["K"] - fc["Sedimentator_BeltFilter"]["K"]
    fc["Sedimentator_LiqEff"]["P-PO4"]     = fc["PROC_Sedimentator"]["P-PO4"] - fc["Sedimentator_BeltFilter"]["P-PO4"]
    fc["Sedimentator_LiqEff"]["N-NH3"]     = fc["PROC_Sedimentator"]["N-NH3"] - fc["Sedimentator_BeltFilter"]["N-NH3"]
    fc["Sedimentator_LiqEff"]["Ca_ion"]    = fc["PROC_Sedimentator"]["Ca_ion"] - fc["Sedimentator_BeltFilter"]["Ca_ion"]
    fc["Sedimentator_LiqEff"]["K_ion"]     = fc["PROC_Sedimentator"]["K_ion"] - fc["Sedimentator_BeltFilter"]["K_ion"]
    fc["Sedimentator_LiqEff"]["Wa"]        = fc["PROC_Sedimentator"]["Wa"] - fc["Sedimentator_BeltFilter"]["Wa"]  
    
    
    
    fc["BeltFilter_Dryer"]["CSH-P"]     = fc["Sedimentator_BeltFilter"]["CSH-P"]
    fc["BeltFilter_Dryer"]["Rest"]      = fc["Sedimentator_BeltFilter"]["Rest"]
    fc["BeltFilter_Dryer"]["C"]         = fc["Sedimentator_BeltFilter"]["C"]
    fc["BeltFilter_Dryer"]["P"]         = fc["Sedimentator_BeltFilter"]["P"]
    fc["BeltFilter_Dryer"]["N"]         = fc["Sedimentator_BeltFilter"]["N"]
    fc["BeltFilter_Dryer"]["Ca"]        = fc["Sedimentator_BeltFilter"]["Ca"]
    fc["BeltFilter_Dryer"]["K"]         = fc["Sedimentator_BeltFilter"]["K"]
    
    fc["BeltFilter_Dryer"]["Wa"]        = (fc["BeltFilter_Dryer"]["CSH-P"]+fc["BeltFilter_Dryer"]["C"]+fc["BeltFilter_Dryer"]["P"]+fc["BeltFilter_Dryer"]["N"]+fc["BeltFilter_Dryer"]["Ca"]+fc["BeltFilter_Dryer"]["K"])*humidity_BeltFilt/(1-humidity_BeltFilt)   #Computing the amount of water contained in solids
    
    Relative_water_BeltFilter = fc["BeltFilter_Dryer"]["Wa"]/fc["Sedimentator_BeltFilter"]["Wa"] #Ratio water solids / water liquid phase
    
    fc["BeltFilter_Dryer"]["P-PO4"]     = fc["Sedimentator_BeltFilter"]["P-PO4"]*Relative_water_BeltFilter
    fc["BeltFilter_Dryer"]["N-NH3"]     = fc["Sedimentator_BeltFilter"]["N-NH3"]*Relative_water_BeltFilter
    fc["BeltFilter_Dryer"]["Ca_ion"]    = fc["Sedimentator_BeltFilter"]["Ca_ion"]*Relative_water_BeltFilter
    fc["BeltFilter_Dryer"]["K_ion"]     = fc["Sedimentator_BeltFilter"]["K_ion"]*Relative_water_BeltFilter
    
    fc["BeltFilter_LiqEff"]["CSH-P"]     = fc["Sedimentator_BeltFilter"]["CSH-P"] - fc["BeltFilter_Dryer"]["CSH-P"]
    fc["BeltFilter_LiqEff"]["Rest"]      = fc["Sedimentator_BeltFilter"]["Rest"] - fc["BeltFilter_Dryer"]["Rest"]
    fc["BeltFilter_LiqEff"]["C"]         = fc["Sedimentator_BeltFilter"]["C"] - fc["BeltFilter_Dryer"]["C"]
    fc["BeltFilter_LiqEff"]["P"]         = fc["Sedimentator_BeltFilter"]["P"] - fc["BeltFilter_Dryer"]["P"]
    fc["BeltFilter_LiqEff"]["N"]         = fc["Sedimentator_BeltFilter"]["N"] - fc["BeltFilter_Dryer"]["N"]
    fc["BeltFilter_LiqEff"]["Ca"]        = fc["Sedimentator_BeltFilter"]["Ca"] - fc["BeltFilter_Dryer"]["Ca"]
    fc["BeltFilter_LiqEff"]["K"]         = fc["Sedimentator_BeltFilter"]["K"] - fc["BeltFilter_Dryer"]["K"]
    fc["BeltFilter_LiqEff"]["P-PO4"]     = fc["Sedimentator_BeltFilter"]["P-PO4"] - fc["BeltFilter_Dryer"]["P-PO4"]
    fc["BeltFilter_LiqEff"]["N-NH3"]     = fc["Sedimentator_BeltFilter"]["N-NH3"] - fc["BeltFilter_Dryer"]["N-NH3"]
    fc["BeltFilter_LiqEff"]["Ca_ion"]    = fc["Sedimentator_BeltFilter"]["Ca_ion"] - fc["BeltFilter_Dryer"]["Ca_ion"]
    fc["BeltFilter_LiqEff"]["K_ion"]     = fc["Sedimentator_BeltFilter"]["K_ion"] - fc["BeltFilter_Dryer"]["K_ion"]
    fc["BeltFilter_LiqEff"]["Wa"]        = fc["Sedimentator_BeltFilter"]["Wa"] - fc["BeltFilter_Dryer"]["Wa"]  
    
    
    fc["Dryer_VapEff"]["Wa"]       = fc["BeltFilter_LiqEff"]["Wa"]
    
    fc["Dryer_PEff"]["Wa"]        = fc["BeltFilter_Dryer"]["Wa"]-fc["Dryer_VapEff"]["Wa"] 
    fc["Dryer_PEff"]["CSH-P"]     = fc["BeltFilter_Dryer"]["CSH-P"]
    fc["Dryer_PEff"]["Rest"]      = fc['BeltFilter_Dryer']["Rest"]
    fc["Dryer_PEff"]["C"]         = fc["BeltFilter_Dryer"]["C"]
    fc["Dryer_PEff"]["P"]         = fc["BeltFilter_Dryer"]["P"]
    fc["Dryer_PEff"]["N"]         = fc["BeltFilter_Dryer"]["N"]
    fc["Dryer_PEff"]["Ca"]        = fc["BeltFilter_Dryer"]["Ca"]
    fc["Dryer_PEff"]["K"]         = fc["BeltFilter_Dryer"]["K"]
    fc["Dryer_PEff"]["P-PO4"]     = fc["BeltFilter_Dryer"]["P-PO4"]
    fc["Dryer_PEff"]["N-NH3"]     = fc["BeltFilter_Dryer"]["N-NH3"]
    fc["Dryer_PEff"]["Ca_ion"]    = fc["BeltFilter_Dryer"]["Ca_ion"]
    fc["Dryer_PEff"]["K_ion"]     = fc["BeltFilter_Dryer"]["K_ion"]

    
    
    
    #for i in total_elements.keys():
        ##if i!="Struvite":
        #fc["MULTIFORM_LiqEff"][i] =fc["MULTIFORM_LiqEff"][i]
    
    
    for i in nodes_list:
        F[i] = sum(fc[i][ii] for ii in total_elements)
        
    for i in nodes_list:
        if i!="Src1_PROC":
            for ii in total_elements.keys():
                x[i][ii] = fc[i][ii]/F[i]
    
    #Checks
    checks_store = ['OK', 'FAIL']
    checks_F = []
    checks_x = []
    
    #    for i in total_elements.keys():
    #        if abs(fc["Src1_MULTIFORM"][i] + fc["Src2Mixer"][i] + fc["Src3FBR"][i] - (fc["HydrocycloneSink1"][i] + fc["MULTIFORM_LiqEff"][i] + fc["DryerMULTIFORM_VaporEff"][i] + fc["DryerMULTIFORM_Peff"][i])) <= 0.005:
    #            print("fc check", i, "OK")
    #        else:
    #            print("fc check", i, "FAIL")
    
    if abs(F["Src1_PROC"] + F["Src2_PROC"] - (F["Sedimentator_LiqEff"] + F["BeltFilter_LiqEff"] + F["Dryer_VapEff"] + F["Dryer_PEff"])) <= 0.005: #F["HydrocycloneSink1"] +
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
    # ECONOMICS/DESIGN
    #-------------------------Residence time----------------------------------
    total_time              = 3600*2 #doi: 10.2166/wst.2011.061
    
    #-------------------------CSTR cost----------------------------------
    mixing_operation='Slurries' #Operation types: 'Blending','Homogeneous reaction','Reaction with heat transfer','Liquid-liquid mixtures','Liquid-gas mixtures','Slurries'
    CSTR_results = CSTR_investment_cost (F_ini, total_time, mixing_operation)
    
    CSTR_V = CSTR_results['reactor_V']
    CSTR_D = CSTR_results['reactor_D']
    CSTR_L = CSTR_results['reactor_L']
    agitator_power = CSTR_results['agitator_power']
    n_CSTR = CSTR_results['n_reactors']
    CSTR_cost_2016 = CSTR_results['CSTR_cost_2016']
    CSTR_operation_cost_2016 = CSTR_results['CSTR_operation_cost_2016']
    
    
    #-------------------------Clarifier----------------------------------
    HRT = 3600 #s
    Clarifier_results = vessel_design_cost(F["PROC_Sedimentator"], HRT)
    
    Clarifier_V_total = Clarifier_results['Volume']
    #Results
    Clarifier_V = Clarifier_results['V_design']
    Clarifier_D = Clarifier_results['D_design']
    Clarifier_L = Clarifier_results['L_design']
    n_Clarifier = Clarifier_results['n_vessels']
    Clarifier_cost_2016 =  Clarifier_results['vessels_cost_2016']

    
    #--------------------------------Belt Filt-------------------------------------
    BeltFilt_results = BeltFilt_design_cost (fc["BeltFilter_LiqEff"]["Wa"])
    BeltFilt_area = BeltFilt_results['Area']
    n_BeltFilt = BeltFilt_results['n_filt']
    BeltFilt_cost_2016 = BeltFilt_results['BeltFilt_cost_2016']
    
    
    #--------------------------------Dryer-------------------------------------
    ConveyorDryer_results = ConveyorDryer_design_cost (F['BeltFilter_Dryer'], fc["BeltFilter_Dryer"]["Wa"])
    ConveyorDryer_area = ConveyorDryer_results['Area']
    n_ConveyorDryer = ConveyorDryer_results['n_dryer']
    ConveyorDryer_cost_2016 = ConveyorDryer_results['ConveyorDryer_cost_2016']
    ConveyorDryer_operating_cost_2016 = ConveyorDryer_results['ConveyorDryer_operating_cost_2016']
    
    
    
    # ====================================================
    # ####################################################
    # ====================================================
    #SIZE DETERMINATION
    PROC_CSTR_size = CSTR_V # m3

        
    #NUMBER OF UNITS
    n_PROC = n_CSTR

    # ====================================================
    # ####################################################
    # ====================================================
    #Unit cost
    equipment_cost              = (CSTR_cost_2016+ Clarifier_cost_2016 + BeltFilt_cost_2016 + ConveyorDryer_cost_2016)
    PROC_equipment_cost    = 3.15*equipment_cost
    #fixed_capital_cost_2016     = 1.4*physical_plant_cost_2016
    fixed_capital_cost_2016     = PROC_equipment_cost

    # ====================================================
    # ####################################################
    # ====================================================
    #OPERATION COST (except chemicals)
    if fc["Src1_PROC"]["P-PO4"] < 135:
        PROC_operating_cost_partial = 115.5*fc["Src1_PROC"]["P-PO4"]*3600*24*365
        
    elif fc["Src1_PROC"]["P-PO4"] >662:
        PROC_operating_cost_partial = 67.9*fc["Src1_PROC"]["P-PO4"]*3600*24*365
    
    else:
        PROC_operating_cost_partial = (-0.09*fc["Src1_PROC"]["P-PO4"]+127.19)*3600*24*365
    
    operation_cost_2016_amortized = (PROC_operating_cost_partial + (fixed_capital_cost_2016/ec_param['plant_lifetime']))
    operation_cost_2016_non_amortized = PROC_operating_cost_partial #+ 1.5*labor_cost_2016['labor_cost']

    recovered_P = fc["Dryer_PEff"]["P"]
    recovered_PO4 = fc["Src1_PROC"]["P-PO4"]*(yield_PO4)+fc["Dryer_PEff"]["P-PO4"]
    released_P = fc["Sedimentator_LiqEff"]["P"] + fc["BeltFilter_LiqEff"]["P"]
    released_PO4 = fc["Sedimentator_LiqEff"]["P-PO4"] + fc["BeltFilter_LiqEff"]["P-PO4"]
    released_NH3 = (fc["Sedimentator_LiqEff"]["N-NH3"]+fc["BeltFilter_LiqEff"]["N-NH3"]) #kgN as NH4/s
    released_N = (fc["Sedimentator_LiqEff"]["N"]+fc["BeltFilter_LiqEff"]["N"]) #kg rest of N/s
    
    fraction_recoved_P = recovered_P/(fc["Src1_PROC"]["P"])
    fraction_recoved_PO4 = recovered_PO4/fc["Src1_PROC"]["P-PO4"]
    fraction_recoved_TP = (recovered_P+recovered_PO4)/(fc["Src1_PROC"]["P-PO4"]+fc["Src1_PROC"]["P"])
    fraction_released_P = released_P/(fc["Src1_PROC"]["P"])
    fraction_released_PO4 = released_PO4/(fc["Src1_PROC"]["P-PO4"])
    fraction_released_TP = (released_P+released_PO4)/(fc["Src1_PROC"]["P-PO4"]+fc["Src1_PROC"]["P"])
    PO4_conc_released = (fc["Sedimentator_LiqEff"]["P-PO4"] + fc["BeltFilter_LiqEff"]["P-PO4"])/MW['P']*MW['PO4']*UnitConv['K_to_mili']/(F["Sedimentator_LiqEff"]/feedstock_parameters['Wa_density'])
    Eutrophication_potential = (released_N*0.42+released_NH3/MW['N']*MW['NH4']*0.33+released_PO4/MW['P']*MW['PO4']*1+released_P*3.06)*365*3600*24/1000 #ton phosphate equivalent / year , potency factors from 'The sustainability metrics. Sustainable Development Progess Metrics recommended for use in the Process Idustries' IChenE
    
    return {'fc':fc,
            'x':x,
            'F':F,
            'check_F':checks_F,
            'checks_x':checks_x,
            'tech':'P-RoC',
            'fixed_capital_cost_2016':fixed_capital_cost_2016,
            'PROC_size':PROC_CSTR_size, 'n_PROC':n_PROC,
            'investment_cost':fixed_capital_cost_2016,
            'equipment_cost':PROC_equipment_cost,
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



