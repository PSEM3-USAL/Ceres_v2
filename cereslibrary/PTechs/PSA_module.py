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
def PSA_module(F_ini,fc_ini,x_ini, RINs_methane,elements_wet, elements_dry, nutrients, feedstock_parameters):
    import numpy as np
    import pandas as pd


    # SETS, PARAMETERS AND NODES DATA ACQUISITION
    nodes_matrix_PSA            = pd.read_csv('nodes/nodes_PSA.csv', sep=",", header=None)
    nodes_PSA                   = np.array(nodes_matrix_PSA[0].dropna())
    process_elements_matrix_PSA = pd.read_csv('process_elements/process_elements_digestor.csv', sep=",", header=0)

    from global_parameters_module import UnitConv, MW, c_p_liq_sol, dH_vap_0, Tc, Tb, dH_f, dH_c, c_p_v_1, c_p_v_2, c_p_v_3, c_p_v_4, coef_vapor_pressure_1, coef_vapor_pressure_2, coef_vapor_pressure_3, CEI, price, nu_p, k_p, n_watson, epsilon, T_amb, P_ref, density, latent_heat_evap, nat_gas_heat_value
    #from feedstock_input_module import elements_wet, elements_dry, nutrients, feedstock_parameters
    from economic_parameters_module import ec_param



    # MANURE DATA ACQUISITION AND COMPOSITION SETTLEMENT   
    gases_comp_PSA   = np.array(process_elements_matrix_PSA["Component"].dropna())
    gases_conc_PSA   = np.full((len(gases_comp_PSA)), 0)
    gases_PSA        = dict(zip(gases_comp_PSA, gases_conc_PSA))

    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================  
    # TOTAL ELEMENTS
    total_elements_PSA = {**elements_wet,**gases_PSA} # Merge the two dictionaries

    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================
    # PARAMETERS
    T_biogas = 25
    Eff_comp_PSA = 0.8
    P_PSA = 3.185
    T_PSA_bed = 34.572
    eff_stabilized =0.65 #efficiency of zeolite when stabilized over initial performance /0.65/ ;
    Eff_PSA_CH4 = 0.05
    Eff_PSA_CO2 = 0.99
    t_cycle = 20 # tiempo ciclo zeolita min /20/
    n_units_parallel = 2 # units in parallel /2/
    zeo_lifetime = 5 # zeolites lifetime years /5/



    # VARIABLES DEFINITION (IN NESTED DICTIONARIES) (INITIALIZATION)
    nodes_list_PSA              = nodes_PSA.tolist()
    initialization_comp_PSA     = total_elements_PSA #["Wa", "C", "NH3", "PO4", "Ca_ion", "K_ion"]
    initialization_nan_PSA      = np.full((len(initialization_comp_PSA)), 0.00)

    fc_PSA = {key: dict(zip(initialization_comp_PSA,initialization_nan_PSA)) for key in nodes_list_PSA}

    x_PSA = {key: dict(zip(initialization_comp_PSA,initialization_nan_PSA)) for key in nodes_list_PSA}

    F_PSA = {key: np.nan for key in nodes_list_PSA}

    T_PSA = {key: np.nan for key in nodes_list_PSA}


    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================     
    # MASS BALANCE
    for i in total_elements_PSA.keys():
        x_PSA["Purification_Comp1"][i]  = x_ini[i]
        fc_PSA["Purification_Comp1"][i] = fc_ini[i]
        #x_PSA["Purification_Comp1"][i]  = x["ADEval_Biogas"][i]
        #fc_PSA["Purification_Comp1"][i] = fc_kgs["ADEval_Biogas"][i]
        
    T_PSA["Comp1_HX1"] = ((T_biogas+273)+(1/Eff_comp_PSA)*(T_biogas+273)*((P_PSA+0.001)**(0.4/1.4)-1))-273

    y_gas_CH4 = fc_PSA["Purification_Comp1"]['CH4']/MW['CH4']/(fc_PSA["Purification_Comp1"]['CH4']/MW['CH4']+fc_PSA["Purification_Comp1"]['CO2']/MW['CO2'])
    y_gas_CO2 = fc_PSA["Purification_Comp1"]['CO2']/MW['CO2']/(fc_PSA["Purification_Comp1"]['CH4']/MW['CH4']+fc_PSA["Purification_Comp1"]['CO2']/MW['CO2'])
    MW_gas = MW['CO2']*y_gas_CO2+MW['CH4']*y_gas_CH4

    #  W('Comp5')*(MW_gas_cp5+0.001)=E=sum(J,fc(J,'Spl1','Comp5'))*(8.314*1.4*(1/Eff_comp_PSA)*(T('Spl1','Comp5')+273))*((P_PSA+0.001)**(0.4/1.4)-1)/((1.4-1))
    W_Comp1 = (sum(fc_PSA["Purification_Comp1"][i] for i in total_elements_PSA)*(8.314*1.4*(1/Eff_comp_PSA)*(T_biogas+273))*((P_PSA+0.001)**(0.4/1.4)-1)/((1.4-1)))/MW_gas

    for i in total_elements_PSA.keys():
        x_PSA["Comp1_HX1"][i]  = x_PSA["Purification_Comp1"][i]
        fc_PSA["Comp1_HX1"][i] = fc_PSA["Purification_Comp1"][i]
        
    for i in total_elements_PSA.keys():
        x_PSA["HX1_MS1"][i]  = x_PSA["Comp1_HX1"][i]
        fc_PSA["HX1_MS1"][i] = fc_PSA["Comp1_HX1"][i]
        
    Q_HX1 = fc_PSA["Comp1_HX1"]["CO2"]*(1/MW['CO2'])*(c_p_v_1['CO2']*(T_PSA_bed-T_PSA["Comp1_HX1"] )+
    (1/2)*c_p_v_2['CO2']*((T_PSA_bed+273)**2-(T_PSA["Comp1_HX1"]+273)**2)+
    (1/3)*c_p_v_3['CO2']*((T_PSA_bed+273)**3-(T_PSA["Comp1_HX1"]+273)**3)+
    (1/4)*c_p_v_4['CO2']*((T_PSA_bed+273)**4-(T_PSA["Comp1_HX1"]+273)**4))+fc_PSA["Comp1_HX1"]["CH4"]*(1/MW['CH4'])*(c_p_v_1['CH4']*(T_PSA_bed-T_PSA["Comp1_HX1"] )+
    (1/2)*c_p_v_2['CH4']*((T_PSA_bed+273)**2-(T_PSA["Comp1_HX1"]+273)**2)+
    (1/3)*c_p_v_3['CH4']*((T_PSA_bed+273)**3-(T_PSA["Comp1_HX1"]+273)**3)+
    (1/4)*c_p_v_4['CH4']*((T_PSA_bed+273)**4-(T_PSA["Comp1_HX1"]+273)**4))

    P_CO2 =E= y_gas_CO2*P_PSA;

    T_PSA["HX1_MS1"] = T_PSA_bed

    qm = -3.15551E-02*T_PSA_bed + 5.02915
    Klagmuir = 1.63070E-03*T_PSA_bed*T_PSA_bed - 3.68662E-01*T_PSA_bed + 2.73737E+01
    q_ads = qm*Klagmuir*P_CO2/(1+ Klagmuir*P_CO2)

    MZeolita = n_units_parallel*t_cycle*60*Eff_PSA_CO2*fc_PSA["HX1_MS1"]['CO2']/MW['CO2']/(0.001*q_ads*eff_stabilized)

    fc_PSA["MS1_CO2cap"]['CO2'] = Eff_PSA_CO2*fc_PSA["HX1_MS1"]['CO2']
    fc_PSA["MS1_CO2cap"]['CH4'] = Eff_PSA_CO2*fc_PSA["HX1_MS1"]['CH4']

    fc_PSA["MS1_CH4pur"]['CO2'] = (1-Eff_PSA_CO2)*fc_PSA["HX1_MS1"]['CO2']
    fc_PSA["MS1_CH4pur"]['CH4'] = (1-Eff_PSA_CO2)*fc_PSA["HX1_MS1"]['CH4']
    
    for i in nodes_list_PSA:
        F_PSA[i] = sum(fc_PSA[i][ii] for ii in total_elements_PSA)
        
    for i in nodes_list_PSA:
        if i!="Purification_Comp1":
            for ii in total_elements_PSA.keys():
                x_PSA[i][ii] = fc_PSA[i][ii]/F_PSA[i]
                
    #Checks
    checks_store = ['OK', 'FAIL']
    checks_F = []
    checks_x = []
    
    #    for i in total_elements.keys():
    #        if abs(fc["Src1_MULTIFORM"][i] + fc["Src2Mixer"][i] + fc["Src3FBR"][i] - (fc["HydrocycloneSink1"][i] + fc["MULTIFORM_LiqEff"][i] + fc["DryerMULTIFORM_VaporEff"][i] + fc["DryerMULTIFORM_Peff"][i])) <= 0.005:
    #            print("fc check", i, "OK")
    #        else:
    #            print("fc check", i, "FAIL")
    
    if abs(F_PSA["Purification_Comp1"] - (F_PSA["MS1_CO2cap"] + F_PSA["MS1_CH4pur"])) <= 0.005: #F["HydrocycloneSink1"] +
    #        print("F check OK")
        checks_F = checks_store[0]
    else:
    #        print("F check FAIL")
        checks_F = checks_store[1]
    
    for i in nodes_list_PSA:
        if abs(1-sum(x_PSA[i][ii] for ii in total_elements_PSA)) <= 0.005:
    #            print("x check", i, "OK")
            checks_x = checks_store[0]
        else:
    #            print("x check", i, "FAIL")
            checks_x = checks_store[1]
    print("\n\n")

    #ECONOMICS
    Pelect = 0.06 #            por kwh      /0.06/
    #RINs_methane = 1.2502995248799291 #USD per kg (liquified)
    PZeolita = 5 #    per kg             /5/
    Cost_Comp1 = 335.27*W_Comp1+36211
    Cost_PSA = ec_param['plant_lifetime']/zeo_lifetime*MZeolita*PZeolita
    OperationCost_PSA = ec_param['plant_lifetime']/zeo_lifetime*MZeolita*PZeolita/ec_param['plant_lifetime']

    EquipmentCost = Cost_Comp1+Cost_PSA

    PSA_OperatingCost = W_Comp1*3600*24*365/3600*Pelect + Cost_PSA/ec_param['plant_lifetime']
    PSA_OperatingCost_amortized = ((EquipmentCost/ec_param['plant_lifetime'])+PSA_OperatingCost)
    Methane_benefits = fc_PSA["MS1_CH4pur"]['CH4']*RINs_methane
        
    return {'fc':fc_PSA,
            'x':x_PSA,
            'F':F_PSA,
            'check_F':checks_F,
            'checks_x':checks_x,
            'tech':'PSA for biomethane generation',
            'fixed_capital_cost_2016':EquipmentCost,
            'investment_cost':EquipmentCost,
            'equipment_cost':EquipmentCost,
            'operation_cost_2016_amortized':PSA_OperatingCost_amortized,
            'operation_cost_2016_non_amortized':PSA_OperatingCost,
            'Methane_benefits':Methane_benefits,
            'Struvite_benefits':0,
            'recovered_P':0,
            'recovered_PO4':0,
            'released_P':0,
            'released_PO4':0,
            'fraction_recoved_P':0,
            'fraction_recoved_PO4':0,
            'fraction_released_P':0,
            'fraction_released_PO4':0,
            'fraction_recoved_TP':0,
            'fraction_released_TP':0,
            'PO4_conc_released':0,   
            'Eutrophication_potential':0,  
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



