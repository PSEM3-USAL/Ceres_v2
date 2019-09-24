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
def turbine_module(F_ini,fc_ini,x_ini, MW_biogas_HX2, REC_MWh,elements_wet, elements_dry, nutrients, feedstock_parameters):
    import numpy as np
    import pandas as pd
    
    REC_kwh = REC_MWh/1000

    # MW_biogas_HX2
    # SETS, PARAMETERS AND NODES DATA ACQUISITION
    nodes_matrix_turb            = pd.read_csv('nodes/nodes_turbine.csv', sep=",", header=None)
    nodes_turb                   = np.array(nodes_matrix_turb[0].dropna())
    process_elements_matrix_turb = pd.read_csv('process_elements/process_elements_digestor.csv', sep=",", header=0)

    from global_parameters_module import UnitConv, MW, c_p_liq_sol, dH_vap_0, Tc, Tb, dH_f, dH_c, c_p_v_1, c_p_v_2, c_p_v_3, c_p_v_4, coef_vapor_pressure_1, coef_vapor_pressure_2, coef_vapor_pressure_3, CEI, price, nu_p, k_p, n_watson, epsilon, T_amb, P_ref, density, latent_heat_evap, nat_gas_heat_value
    #from feedstock_input_module import elements_wet, elements_dry, nutrients, feedstock_parameters
    from economic_parameters_module import ec_param


    # MANURE DATA ACQUISITION AND COMPOSITION SETTLEMENT   
    gases_comp_turb   = np.array(process_elements_matrix_turb["Component"].dropna())
    gases_conc_turb   = np.full((len(gases_comp_turb)), 0)
    gases_turb        = dict(zip(gases_comp_turb, gases_conc_turb))

    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================  
    # TOTAL ELEMENTS
    total_elements_turb = {**elements_wet,**gases_turb} # Merge the two dictionaries

    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================
    # AD PARTICULAR PARAMETERS
    T_PSA = 25 #          Working temperature of PSA systhems /25/
    P_PSA = 4.5
    Pressure_GT_in = 8.216
    Pressure_GT_out = 1
    excess = 1.2

    # VARIABLES DEFINITION (IN NESTED DICTIONARIES) (INITIALIZATION)
    nodes_list_turb              = nodes_turb.tolist()
    initialization_comp_turb     = total_elements_turb #["Wa", "C", "NH3", "PO4", "Ca_ion", "K_ion"]
    initialization_nan_turb      = np.full((len(initialization_comp_turb)), 0.00)

    fc_turb = {key: dict(zip(initialization_comp_turb,initialization_nan_turb)) for key in nodes_list_turb}

    x_turb = {key: dict(zip(initialization_comp_turb,initialization_nan_turb)) for key in nodes_list_turb}

    y_turb  = {key: dict(zip(initialization_comp_turb,initialization_nan_turb)) for key in nodes_list_turb}

    F_turb = {key: np.nan for key in nodes_list_turb}

    T_turb = {key: np.nan for key in nodes_list_turb}


    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================     
    # MASS BALANCE
    for i in total_elements_turb.keys():
        x_turb["Src1_Comp1"][i]  = x_ini[i]
        fc_turb["Src1_Comp1"][i] = fc_ini[i]
        #fc_turb["Src1_Comp1"][i] = fc_purif["React2_Sink4Biogas"][i]
        
    for i in total_elements_turb.keys():
        fc_turb["Comp1_Furnance"][i] = fc_turb["Src1_Comp1"][i]
        
    W_Comp1 = (3*(1/0.85)*sum(fc_turb["Src1_Comp1"][i] for i in total_elements_turb)*(8.314*1.4*(T_PSA+273))*(((Pressure_GT_in)/P_PSA)**(0.33*(0.4/1.4))-1))/((MW_biogas_HX2)*(1.4-1))

    T_turb["Comp1_Furnance"] = ((T_PSA+273)+(1/0.85)*(T_PSA+273)*(((Pressure_GT_in)/P_PSA)**(0.33*(0.4/1.4))-1))-273

    fc_turb["Src2_Comp2"]['O2'] = excess*2*(fc_turb["Comp1_Furnance"]['CH4']/MW['CH4'])*MW['O2']
    fc_turb["Src2_Comp2"]['N2'] = 3.76*fc_turb["Src2_Comp2"]['O2']/MW['O2']*MW['N2']

    T_turb["Src2_Comp2"] = 20

    for i in total_elements_turb.keys():
        fc_turb["Comp2_Furnance"][i] = fc_turb["Src2_Comp2"][i]
        
    fc_turb["Furnace_TurbGas"]['CO2'] = fc_turb["Comp1_Furnance"]['CO2'] + fc_turb["Comp1_Furnance"]['CH4']/MW['CH4']*MW['CO2']
    fc_turb["Furnace_TurbGas"]['Wa'] = 2*fc_turb["Comp1_Furnance"]['CH4']/MW['CH4']*MW['Wa']
    fc_turb["Furnace_TurbGas"]['N2'] = fc_turb["Comp2_Furnance"]['N2']
    fc_turb["Furnace_TurbGas"]['O2'] = fc_turb["Comp2_Furnance"]['O2']*(excess-1)

    W_Comp2 = 3*(1/0.85)*(sum(fc_turb["Src2_Comp2"][i] for i in total_elements_turb))*(8.314*1.4*(20+273))*(((Pressure_GT_in)/1)**(0.33*(0.4/1.4))-1)/((29)*(1.4-1))

    T_turb["Comp2_Furnance"] = ((20+273)+(1/0.85)*(20+273)*(((Pressure_GT_in)/1)**(0.33*(0.4/1.4))-1))-273

    Q_react1 = sum(fc_turb["Comp1_Furnance"][ii]*(dH_f[ii]+1/MW[ii]*(c_p_v_1[ii]*(T_turb["Comp1_Furnance"]-25)+
                                                    1/2*c_p_v_2[ii]*((T_turb["Comp1_Furnance"]+273)**2-(25+273)**2 )+
                                                    1/3*c_p_v_3[ii]*((T_turb["Comp1_Furnance"]+273)**3-(25+273)**3 )+
                                                    1/4*c_p_v_4[ii]*((T_turb["Comp1_Furnance"]+273)**4-(25+273)**4))) for ii in gases_turb)+fc_turb["Comp1_Furnance"]['Wa']*(dH_f['Wa']+1/MW['Wa']*(c_p_v_1['Wa']*(T_turb["Comp1_Furnance"]-25)+
                                                    1/2*c_p_v_2['Wa']*((T_turb["Comp1_Furnance"]+273)**2-(25+273)**2 )+
                                                    1/3*c_p_v_3['Wa']*((T_turb["Comp1_Furnance"]+273)**3-(25+273)**3 )+
                                                    1/4*c_p_v_4['Wa']*((T_turb["Comp1_Furnance"]+273)**4-(25+273)**4)))

    Q_react2 = (fc_turb["Comp1_Furnance"]['O2']*(dH_f['O2']+1/MW['O2']*(c_p_v_1['O2']*(T_turb["Comp2_Furnance"]-25)+
                                                    1/2*c_p_v_2['O2']*((T_turb["Comp2_Furnance"]+273)**2-(25+273)**2 )+
                                                    1/3*c_p_v_3['O2']*((T_turb["Comp2_Furnance"]+273)**3-(25+273)**3 )+
                                                    1/4*c_p_v_4['O2']*((T_turb["Comp2_Furnance"]+273)**4-(25+273)**4))))+fc_turb["Comp2_Furnance"]['N2']*(dH_f['N2']+1/MW['N2']*(c_p_v_1['N2']*(T_turb["Comp2_Furnance"]-25)+
                                                    1/2*c_p_v_2['N2']*((T_turb["Comp2_Furnance"]+273)**2-(25+273)**2 )+
                                                    1/3*c_p_v_3['N2']*((T_turb["Comp2_Furnance"]+273)**3-(25+273)**3 )+
                                                    1/4*c_p_v_4['N2']*((T_turb["Comp2_Furnance"]+273)**4-(25+273)**4)))


    Q_prod = Q_react1+Q_react2

    # Q_prod = sum(fc_turb["Furnace_TurbGas"][ii]*(dH_f[ii]+1/MW[ii]*(c_p_v_1[ii]*(T_turb["Furnace_TurbGas"]-25)+
    #                                                   1/2*c_p_v_2[ii]*((T_turb["Furnace_TurbGas"]+273)**2-(25+273)**2 )+
    #                                                   1/3*c_p_v_3[ii]*((T_turb["Furnace_TurbGas"]+273)**3-(25+273)**3 )+
    #                                                   1/4*c_p_v_4[ii]*((T_turb["Furnace_TurbGas"]+273)**4-(25+273)**4))) for ii in gases_turb)+fc_turb["Furnace_TurbGas"]['Wa']*(dH_f['Wa']+1/MW['Wa']*(c_p_v_1['Wa']*(T_turb["Furnace_TurbGas"]-25)+
    #                                                   1/2*c_p_v_2['Wa']*((T_turb["Furnace_TurbGas"]+273)**2-(25+273)**2 )+
    #                                                   1/3*c_p_v_3['Wa']*((T_turb["Furnace_TurbGas"]+273)**3-(25+273)**3 )+
    #                                                   1/4*c_p_v_4['Wa']*((T_turb["Furnace_TurbGas"]+273)**4-(25+273)**4)))

    def fun_T_x(T_x):
        return sum(fc_turb["Furnace_TurbGas"][ii]*(dH_f[ii]+1/MW[ii]*(c_p_v_1[ii]*(T_x-25)+
                                                    1/2*c_p_v_2[ii]*((T_x+273)**2-(25+273)**2 )+
                                                    1/3*c_p_v_3[ii]*((T_x+273)**3-(25+273)**3 )+
                                                    1/4*c_p_v_4[ii]*((T_x+273)**4-(25+273)**4))) for ii in gases_turb)+fc_turb["Furnace_TurbGas"]['Wa']*(dH_f['Wa']+1/MW['Wa']*(c_p_v_1['Wa']*(T_x-25)+
                                                    1/2*c_p_v_2['Wa']*((T_x+273)**2-(25+273)**2 )+
                                                    1/3*c_p_v_3['Wa']*((T_x+273)**3-(25+273)**3 )+
                                                    1/4*c_p_v_4['Wa']*((T_x+273)**4-(25+273)**4)))-Q_prod

    # Q_prod = sum(fc_turb["Comp1_Furnance"][ii]*(dH_f[ii]+1/MW[ii]*(c_p_v_1[ii]*(T_x-25)+
    #                                                   1/2*c_p_v_2[ii]*((T_x+273)**2-(25+273)**2 )+
    #                                                   1/3*c_p_v_3[ii]*((T_x+273)**3-(25+273)**3 )+
    #                                                   1/4*c_p_v_4[ii]*((T_x+273)**4-(25+273)**4))) for ii in gases_turb)+fc_turb["Furnace_TurbGas"]['Wa']*(dH_f['Wa']+1/MW['Wa']*(c_p_v_1['Wa']*(T_x-25)+
    #                                                   1/2*c_p_v_2['Wa']*((T_x+273)**2-(25+273)**2 )+
    #                                                   1/3*c_p_v_3['Wa']*((T_x+273)**3-(25+273)**3 )+
    #                                                   1/4*c_p_v_4['Wa']*((T_x+273)**4-(25+273)**4)))-Q_prod

    from scipy import optimize
    sol_T_x = optimize.root(fun_T_x, 0)

    if sol_T_x.x < 1100:
        T_turb["Furnace_TurbGas"] = sol_T_x.x
    else:
        T_turb["Furnace_TurbGas"] = 1100
        
    W_TurbGas = (-0.85*sum(fc_turb["Furnace_TurbGas"][i] for i in total_elements_turb)*(8.314*1.3*(T_turb["Furnace_TurbGas"]*10+273))*(((Pressure_GT_out)/(Pressure_GT_in))**(0.3/1.3)-1))/((28)*(1.3-1))

    Q_HX1 = -2*sum(fc_turb["Comp1_Furnance"][ii]*(1/MW[ii]*(c_p_v_1[ii]*(T_turb["Comp1_Furnance"]-25)+
                                                    1/2*c_p_v_2[ii]*((T_turb["Comp1_Furnance"]+273)**2-(25+273)**2 )+
                                                    1/3*c_p_v_3[ii]*((T_turb["Comp1_Furnance"]+273)**3-(25+273)**3 )+
                                                    1/4*c_p_v_4[ii]*((T_turb["Comp1_Furnance"]+273)**4-(25+273)**4))) for ii in gases_turb)+fc_turb["Comp1_Furnance"]['Wa']*(dH_f['Wa']+1/MW['Wa']*(c_p_v_1['Wa']*(T_turb["Comp1_Furnance"]-25)+
                                                    1/2*c_p_v_2['Wa']*((T_turb["Comp1_Furnance"]+273)**2-(25+273)**2 )+
                                                    1/3*c_p_v_3['Wa']*((T_turb["Comp1_Furnance"]+273)**3-(25+273)**3 )+
                                                    1/4*c_p_v_4['Wa']*((T_turb["Comp1_Furnance"]+273)**4-(25+273)**4)))

    Q_HX2 = -2*(fc_turb["Comp1_Furnance"]['O2']*(dH_f['O2']+1/MW['O2']*(c_p_v_1['O2']*(T_turb["Comp2_Furnance"]-25)+
                                                    1/2*c_p_v_2['O2']*((T_turb["Comp2_Furnance"]+273)**2-(25+273)**2 )+
                                                    1/3*c_p_v_3['O2']*((T_turb["Comp2_Furnance"]+273)**3-(25+273)**3 )+
                                                    1/4*c_p_v_4['O2']*((T_turb["Comp2_Furnance"]+273)**4-(25+273)**4))))+fc_turb["Comp2_Furnance"]['N2']*(dH_f['N2']+1/MW['N2']*(c_p_v_1['N2']*(T_turb["Comp2_Furnance"]-25)+
                                                    1/2*c_p_v_2['N2']*((T_turb["Comp2_Furnance"]+273)**2-(25+273)**2 )+
                                                    1/3*c_p_v_3['N2']*((T_turb["Comp2_Furnance"]+273)**3-(25+273)**3 )+
                                                    1/4*c_p_v_4['N2']*((T_turb["Comp2_Furnance"]+273)**4-(25+273)**4)))
                                                    
    for i in nodes_list_turb:
        F_turb[i] = sum(fc_turb[i][ii] for ii in total_elements_turb)
        
    for i in nodes_list_turb:
        if i!="Src1_Comp1":
            for ii in total_elements_turb.keys():
                x_turb[i][ii] = fc_turb[i][ii]/F_turb[i]
                
    #Checks
    checks_store = ['OK', 'FAIL']
    checks_F = []
    checks_x = []
    
    #    for i in total_elements.keys():
    #        if abs(fc["Src1_MULTIFORM"][i] + fc["Src2Mixer"][i] + fc["Src3FBR"][i] - (fc["HydrocycloneSink1"][i] + fc["MULTIFORM_LiqEff"][i] + fc["DryerMULTIFORM_VaporEff"][i] + fc["DryerMULTIFORM_Peff"][i])) <= 0.005:
    #            print("fc check", i, "OK")
    #        else:
    #            print("fc check", i, "FAIL")
    
    if abs(F_turb["Src1_Comp1"] + F_turb["Src2_Comp2"] - (F_turb["Furnace_TurbGas"])) <= 0.005: #F["HydrocycloneSink1"] +
    #        print("F check OK")
        checks_F = checks_store[0]
    else:
    #        print("F check FAIL")
        checks_F = checks_store[1]
    
    for i in nodes_list_turb:
        if abs(1-sum(x_turb[i][ii] for ii in total_elements_turb)) <= 0.005:
    #            print("x check", i, "OK")
            checks_x = checks_store[0]
        else:
    #            print("x check", i, "FAIL")
            checks_x = checks_store[1]
    print("\n\n")


    #ECONOMIC
    Pelect = 0.06 #           por kwh      /0.06/

    #Compresors
    Cost_Comp1 = 335.27*W_Comp1+36211
    Cost_Comp2 = 335.27*W_Comp2+36211

    # Turbine
    Cost_TurbGas = 2281.9*W_TurbGas**0.8001

    #HXs
    T_in_HOT_HX1 = T_turb["Comp1_Furnance"]+273
    T_out_HOT_HX1 = 298+273
    T_in_COLD_HX1 = 293
    T_out_COLD_HX1 = 303
    T_logmeam_HX1 = ((T_in_HOT_HX1-T_out_COLD_HX1)-(T_out_HOT_HX1-T_in_COLD_HX1))/np.log((T_in_HOT_HX1-T_out_COLD_HX1)/(T_out_HOT_HX1-T_in_COLD_HX1))

    h_HOT_HX1 = 0.028
    h_COLD_HX1 = 0.68
    U_HX1 = (1/h_HOT_HX1+1/h_COLD_HX1)**(-1)
    A_HX1=abs(Q_HX1)/(U_HX1*T_logmeam_HX1) #m2
    if A_HX1 < 25: #Almena and Martin, 2016
        Cost_HX1 = 293.59*A_HX1**1.4915
    elif 25 < A_HX1 < 140:
        Cost_HX1 = 1593.8*A_HX1+2584.2
    elif A_HX1 > 140:
        Cost_HX1 = 22234*A_HX1**0.4671

    T_in_HOT_HX2 = T_turb["Comp2_Furnance"]+273
    T_out_HOT_HX2 = 25+273
    T_in_COLD_HX2 = 293
    T_out_COLD_HX2 = 303
    T_logmeam_HX2 = ((T_in_HOT_HX2-T_out_COLD_HX2)-(T_out_HOT_HX2-T_in_COLD_HX2))/np.log((T_in_HOT_HX2-T_out_COLD_HX2)/(T_out_HOT_HX2-T_in_COLD_HX2))

    h_HOT_HX2 = 0.028
    h_COLD_HX2 = 0.68
    U_HX2 = (1/h_HOT_HX2+1/h_COLD_HX2)**(-1)
    A_HX2=abs(Q_HX2)/(U_HX2*T_logmeam_HX2) #m2
    if A_HX2 < 25: #Almena and Martin, 2016
        Cost_HX2 = 293.59*A_HX2**1.4915
    elif 25 < A_HX2 < 140:
        Cost_HX2 = 1593.8*A_HX2+2584.2
    elif A_HX2 > 140:
        Cost_HX2 = 22234*A_HX2**0.4671
        
    Turbine_equipment_cost = Cost_TurbGas + Cost_Comp1 + Cost_Comp2 + Cost_HX1 + Cost_HX2
    Turbine_OperatingCost = (W_Comp1+W_Comp2)*3600*24*365/3600*Pelect
    Turbine_OperatingCost_amortized = ((Turbine_equipment_cost/ec_param['plant_lifetime'])+Turbine_OperatingCost)
    Electricity_benefits = (W_TurbGas)*3600*24*365/3600*REC_kwh
    
    return {'fc':fc_turb,
            'x':x_turb,
            'F':F_turb,
            'check_F':checks_F,
            'checks_x':checks_x,
            'tech':'Gas turbine',
            'fixed_capital_cost_2016':Turbine_equipment_cost,
            'investment_cost':Turbine_equipment_cost,
            'equipment_cost':Turbine_equipment_cost,
            'operation_cost_2016_amortized':Turbine_OperatingCost_amortized,
            'operation_cost_2016_non_amortized':Turbine_OperatingCost,
            'Electricity_benefits':Electricity_benefits,
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



