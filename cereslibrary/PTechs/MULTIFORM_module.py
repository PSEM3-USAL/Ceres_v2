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
def MULTIFORM_module(F_ini,fc_ini,x_ini, elements_wet, elements_dry, nutrients, feedstock_parameters):
    import pandas as pd
    import numpy as np
    
    # IMPORT MODULES
    from global_parameters_module import UnitConv, MW, c_p_liq_sol, dH_vap_0, Tc, Tb, dH_f, dH_c, c_p_v_1, c_p_v_2, c_p_v_3, c_p_v_4, coef_vapor_pressure_1, coef_vapor_pressure_2, coef_vapor_pressure_3, CEI, price, nu_p, k_p, n_watson, epsilon, T_amb, P_ref, density
    #from feedstock_input_module import elements_wet, elements_dry, nutrients, feedstock_parameters, elements_dry_comp
    from equipment_costs.MULTIFORM_cost_module import MULTIFORM_cost_module
    from economic_parameters_module import ec_param
    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================
    # SETS, PARAMETERS AND NODES DATA ACQUISITION
    nodes_matrix            = pd.read_csv('nodes/nodes_MULTIFORM.csv', sep=",", header=None)
    #nodes_matrix            = pd.read_csv('cereslibrary/techmodels/nodes/nodes_FBR.csv', sep=",", header=None)
    nodes       = np.array(nodes_matrix[0].dropna())   
    process_elements_matrix = pd.read_csv('process_elements/process_elements_MULTIFORM.csv', sep=",", header=0)
    #process_elements_matrix = pd.read_csv('cereslibrary/techmodels/process_elements/process_elements_FBR.csv', sep=",", header=0)
    
    #Specefic technology parameters
    Mg_P_ratio          = 2 # Mg/P molar ratio. MgCl2 dose. S Bhuiyan, M. I. H., Mavinic, D. S. and Koch, F. A. Phosphorus recovery from wastewater through struvite formation in fluidized bed reactors: a sustainable approach. Water Science and Technology, 57.2, 175-181, 2008. http://dx.doi.org/ 10.2166/wst.2008.002
    fines_struvite      = 0 #Struvite fines (% in mass of struvite formed)
    struvite_seeds_dose = 0.5E-3 #Stuvite seeds amount kg/L Zhang, T., Li, P., Fang, C. and Jiang, R. Phosphate recovery from animal manure wastewater by struvite crystallization and CO 2 degasificcation reactor. Ecol Chem Eng, 21, 89-99, 2014. http://dx.doi.org/10.2478/eces-2014-0008
    humidity_FBR        = 5 #mass
    g                   = 9.81 #m/s  
    
    
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
#        x["Src1_MULTIFORM"][i]  = total_elements[i]/100
#        fc["Src1_MULTIFORM"][i] = x["Src1_MULTIFORM"][i]*F_ini
        x["Src1_MULTIFORM"][i]  = x_ini[i]
        fc["Src1_MULTIFORM"][i] = fc_ini[i]
        

    fc["Src2_MULTIFORM"]["MgCl2"]         = (fc["Src1_MULTIFORM"]["P-PO4"]/MW["P"])*Mg_P_ratio*MW["MgCl2"]
    
    #fc["Src3Mixer"]["NaOH"]=0 #TEMPORAL

    
    # It is considered that the materia volume is contributed by the digestate. Other contributions are negligible
    # It is considered that volume refered in seeds_amount is the reactor volume
    fc["Src3_MULTIFORM"]["Struvite_seeds"] = (sum(fc["Src1_MULTIFORM"][i] for i in total_elements)/feedstock_parameters['digestate_density'])*struvite_seeds_dose
    

    fc["MULTIFORM_Dryer"]["Struvite_seeds"] = fc["Src3_MULTIFORM"]["Struvite_seeds"]
    # It is consider that it is dry struvite (dry base)
    #fc["FBRSink2"]["Struvite"] = (fc["FBRSink2"]["PO4"]/MW["PO4"])*MW["Struvite"]+fc["FBRSink2"]["Struvite_seeds"]
    
    #Total suspended solids penalization in precipitates formation
    DM_Src1_MULTIFORM_percentage      = 100*(fc["Src1_MULTIFORM"]["P"]+fc["Src1_MULTIFORM"]["N"]+fc["Src1_MULTIFORM"]["C"]+fc["Src1_MULTIFORM"]["Rest"]+fc["Src1_MULTIFORM"]["Ca"]+fc["Src1_MULTIFORM"]["K"])/sum(fc["Src1_MULTIFORM"][ii] for ii in total_elements)
    TSS_Src1_MULTIFORM_percentage     = 0.078*DM_Src1_MULTIFORM_percentage + 0.051
    penalization_index = 1
    #if TSS_Src1_MULTIFORM_percentage <= 0.215:
        #penalization_index = 1
    #else:
        #penalization_index = -7.335E-01*TSS_Src1_MULTIFORM_percentage + 1.158
    
    # Calcium correlations
    Ca_PO4_ratio=(fc["Src1_MULTIFORM"]["Ca_ion"]/MW["Ca"])/(fc["Src1_MULTIFORM"]["P-PO4"]/MW["PO4"])
    fc["MULTIFORM_Dryer"]["Struvite"]   = penalization_index*(fc["Src1_MULTIFORM"]["P-PO4"]*(1-(fines_struvite/100)))*(0.798/(1+(Ca_PO4_ratio*0.576)**2.113))/MW['P']*MW['Struvite']
    fc["MULTIFORM_Dryer"]["CaCO3"]      = penalization_index*(fc["Src1_MULTIFORM"]["Ca_ion"])*(1.020/(1+((Ca_PO4_ratio)*4.097E-01)**(1.029)))/MW['Ca']*MW['CaCO3']
    fc["MULTIFORM_Dryer"]["HAP"]        = penalization_index*(fc["Src1_MULTIFORM"]["Ca_ion"])*(-4.321E-02*Ca_PO4_ratio*Ca_PO4_ratio + 3.128E-01*Ca_PO4_ratio - 3.619E-02)/MW['Ca']/5*MW['HAP']
    
    fc["MULTIFORM_Dryer"]["Rest"]   = fc["Src1_MULTIFORM"]["Rest"]
    fc["MULTIFORM_Dryer"]["C"]      = fc["Src1_MULTIFORM"]["C"]
    fc["MULTIFORM_Dryer"]["N"]      = fc["Src1_MULTIFORM"]["N"]
    fc["MULTIFORM_Dryer"]["P"]      = fc["Src1_MULTIFORM"]["P"]
    fc["MULTIFORM_Dryer"]["Ca"]     = fc["Src1_MULTIFORM"]["Ca"]
    fc["MULTIFORM_Dryer"]["K"]      = fc["Src1_MULTIFORM"]["K"]
    
    fc["MULTIFORM_Dryer"]["Wa"] = ((fc["MULTIFORM_Dryer"]["Struvite"]+fc["MULTIFORM_Dryer"]["CaCO3"] +fc["MULTIFORM_Dryer"]["HAP"]+fc["MULTIFORM_Dryer"]["Rest"]+fc["MULTIFORM_Dryer"]["C"] +fc["MULTIFORM_Dryer"]["N"] +fc["MULTIFORM_Dryer"]["P"]+fc["MULTIFORM_Dryer"]["Ca"] +fc["MULTIFORM_Dryer"]["K"])/(1-(humidity_FBR/100)))*(humidity_FBR/100)
    
    fc["MULTIFORM_LiqEff"]["Wa"]     = fc["Src1_MULTIFORM"]["Wa"]-fc["MULTIFORM_Dryer"]["Wa"]
    Relative_water                  = fc["MULTIFORM_Dryer"]["Wa"]/fc["Src1_MULTIFORM"]["Wa"]
    
    
    fc["MULTIFORM_Dryer"]["P-PO4"]    = Relative_water*(fc["Src1_MULTIFORM"]["P-PO4"]/MW['P']-fc["MULTIFORM_Dryer"]["Struvite"]/MW['Struvite']-3*fc["MULTIFORM_Dryer"]["HAP"]/MW['HAP'])*MW['P']
    if fc["MULTIFORM_Dryer"]["P-PO4"] >= 0:
        pass
    else:
        print('ERROR negative fc["MULTIFORM_Dryer"]["P-PO4"] variable')
    
    fc["MULTIFORM_Dryer"]["N-NH3"]    = Relative_water*(fc["Src1_MULTIFORM"]["N-NH3"]/MW['N']-fc["MULTIFORM_Dryer"]["Struvite"]/MW['Struvite'])*MW['N']
    if fc["MULTIFORM_Dryer"]["N-NH3"] >= 0:
        pass
    else:
        print('ERROR negative fc["MULTIFORM_Dryer"]["N-NH3"] variable')
    
    fc["MULTIFORM_Dryer"]["Mg"]    = Relative_water*(fc["Src2_MULTIFORM"]["MgCl2"]/MW['MgCl2']-fc["MULTIFORM_Dryer"]["Struvite"]/MW['Struvite'])*MW['Mg']
    if fc["MULTIFORM_Dryer"]["Mg"] >= 0:
        pass
    else:
        print('ERROR negative fc["MULTIFORM_Dryer"]["Mg"] variable')
        
    fc["MULTIFORM_Dryer"]["Cl"]     = Relative_water*(fc["Src2_MULTIFORM"]["MgCl2"]/MW["MgCl2"])*2*MW["Cl"]
    
    fc["MULTIFORM_Dryer"]["Ca_ion"] = Relative_water*(fc["Src1_MULTIFORM"]["Ca_ion"]-fc["MULTIFORM_Dryer"]["CaCO3"]/MW['CaCO3']*MW['Ca']-5*fc["MULTIFORM_Dryer"]["HAP"]/MW['HAP']*MW['Ca'])
    if fc["MULTIFORM_Dryer"]["Ca_ion"] >= 0:
        pass
    else:
        print('ERROR negative fc["MULTIFORM_Dryer"]["Ca_ion"] variable')
        
    fc["MULTIFORM_Dryer"]["K_ion"]  = Relative_water*fc["Src1_MULTIFORM"]["K_ion"]
    
    fc["DryerMULTIFORM_VaporEff"]["Wa"]          = fc["MULTIFORM_Dryer"]["Wa"] 
    
    fc["DryerMULTIFORM_Peff"]["Struvite"]  = fc["MULTIFORM_Dryer"]["Struvite"]
    fc["DryerMULTIFORM_Peff"]["HAP"]       = fc["MULTIFORM_Dryer"]["HAP"]
    fc["DryerMULTIFORM_Peff"]["CaCO3"]     = fc["MULTIFORM_Dryer"]["CaCO3"]
    fc["DryerMULTIFORM_Peff"]["Rest"]      = fc["MULTIFORM_Dryer"]["Rest"]
    fc["DryerMULTIFORM_Peff"]["C"]         = fc["MULTIFORM_Dryer"]["C"]
    fc["DryerMULTIFORM_Peff"]["P"]         = fc["MULTIFORM_Dryer"]["P"]
    fc["DryerMULTIFORM_Peff"]["N"]         = fc["MULTIFORM_Dryer"]["N"]
    fc["DryerMULTIFORM_Peff"]["Ca"]        = fc["MULTIFORM_Dryer"]["Ca"]
    fc["DryerMULTIFORM_Peff"]["K"]         = fc["MULTIFORM_Dryer"]["K"]
    fc["DryerMULTIFORM_Peff"]["P-PO4"]     = fc["MULTIFORM_Dryer"]["P-PO4"]
    fc["DryerMULTIFORM_Peff"]["N-NH3"]     = fc["MULTIFORM_Dryer"]["N-NH3"]
    fc["DryerMULTIFORM_Peff"]["Mg"]        = fc["MULTIFORM_Dryer"]["Mg"]
    fc["DryerMULTIFORM_Peff"]["Cl"]        = fc["MULTIFORM_Dryer"]["Cl"]
    fc["DryerMULTIFORM_Peff"]["Ca_ion"]    = fc["MULTIFORM_Dryer"]["Ca_ion"]
    fc["DryerMULTIFORM_Peff"]["K_ion"]     = fc["MULTIFORM_Dryer"]["K_ion"]
    
    
    fc["MULTIFORM_LiqEff"]["P-PO4"]      = fc["Src1_MULTIFORM"]["P-PO4"]-fc["MULTIFORM_Dryer"]["P-PO4"]-fc["MULTIFORM_Dryer"]["Struvite"]/MW['Struvite']*MW['P']-3*fc["MULTIFORM_Dryer"]["HAP"]/MW['HAP']*MW['P']
    fc["MULTIFORM_LiqEff"]["N-NH3"]      = fc["Src1_MULTIFORM"]["N-NH3"]-fc["MULTIFORM_Dryer"]["N-NH3"]-fc["MULTIFORM_Dryer"]["Struvite"]/MW['Struvite']*MW['N']
    fc["MULTIFORM_LiqEff"]["Mg"]         = ((fc["Src1_MULTIFORM"]["MgCl2"]/MW['MgCl2']-fc["MULTIFORM_Dryer"]["Struvite"]/MW['Struvite'])*MW['Mg'])-fc["MULTIFORM_Dryer"]["Mg"]
    fc["MULTIFORM_LiqEff"]["Cl"]         = fc["Src1_MULTIFORM"]["Cl"]-fc["MULTIFORM_Dryer"]["Cl"]
    fc["MULTIFORM_LiqEff"]["Ca_ion"]     = fc["Src1_MULTIFORM"]["Ca_ion"]-fc["MULTIFORM_Dryer"]["Ca_ion"]-5*fc["MULTIFORM_Dryer"]["HAP"]/MW['HAP']*MW['Ca']-fc["MULTIFORM_Dryer"]["CaCO3"]/MW['CaCO3']*MW['Ca']
    fc["MULTIFORM_LiqEff"]["K_ion"]      = fc["Src1_MULTIFORM"]["K_ion"]-fc["MULTIFORM_Dryer"]["K_ion"]
    
    fc["MULTIFORM_LiqEff"]["Rest"]       = fc["Src1_MULTIFORM"]["Rest"]-fc["MULTIFORM_Dryer"]["Rest"]
    fc["MULTIFORM_LiqEff"]["C"]          = fc["Src1_MULTIFORM"]["C"]-fc["MULTIFORM_Dryer"]["C"]
    fc["MULTIFORM_LiqEff"]["K"]          = fc["Src1_MULTIFORM"]["K"]-fc["MULTIFORM_Dryer"]["K"]
    fc["MULTIFORM_LiqEff"]["N"]          = fc["Src1_MULTIFORM"]["N"]-fc["MULTIFORM_Dryer"]["N"]
    fc["MULTIFORM_LiqEff"]["P"]          = fc["Src1_MULTIFORM"]["P"]-fc["MULTIFORM_Dryer"]["P"]
    fc["MULTIFORM_LiqEff"]["Ca"]         = fc["Src1_MULTIFORM"]["Ca"]-fc["MULTIFORM_Dryer"]["Ca"]
    fc["MULTIFORM_LiqEff"]["Struvite"]   = fc["MULTIFORM_Dryer"]["Struvite"]*(fines_struvite/100)
    fc["MULTIFORM_LiqEff"]["HAP"]        = fc["MULTIFORM_Dryer"]["HAP"]*(fines_struvite/100)
    fc["MULTIFORM_LiqEff"]["CaCO3"]      = fc["MULTIFORM_Dryer"]["CaCO3"]*(fines_struvite/100)
    
    
    
    #for i in total_elements.keys():
        ##if i!="Struvite":
        #fc["MULTIFORM_LiqEff"][i] =fc["MULTIFORM_LiqEff"][i]
    
    #fc["HydrocycloneSink1"]["Struvite"] = fc["MULTIFORM_LiqEff"]["Struvite"]
    
    Struvite_recovered = fc["MULTIFORM_Dryer"]["Struvite"]#+fc["MULTIFORM_Dryer"]["HAP"]# + fc["HydrocycloneSink1"]["Struvite"]
    
    #Struvite_benefits = Struvite_recovered*3600*24*334*price_struvite
    
    
    for i in nodes_list:
        F[i] = sum(fc[i][ii] for ii in total_elements)
        
    for i in nodes_list:
        if i!="Src1_MULTIFORM":
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
    
    if abs(F["Src1_MULTIFORM"] + F["Src2_MULTIFORM"] + F["Src3_MULTIFORM"] - (F["MULTIFORM_LiqEff"] + F["DryerMULTIFORM_VaporEff"] + F["DryerMULTIFORM_Peff"])) <= 0.005: #F["HydrocycloneSink1"] +
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
    #pH adjust
    pH_ini=feedstock_parameters['pH']
    pOH_ini=14 - pH_ini
    conc_OH_molar_ini = 10**(-pOH_ini)
    
    pH_fin=9
    pOH_fin=14 - pH_fin
    conc_OH_molar_fin = 10**(-pOH_fin)
    
    delta_conc_OH_molar = conc_OH_molar_fin - conc_OH_molar_ini
    
    NaOH_molar = delta_conc_OH_molar
    
    NaOH_kg = NaOH_molar*MW['NaOH']*F_ini/(feedstock_parameters['digestate_density']*1000)*3600*24*365 #kg/year
     
            
    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================
    
    # ECONOMICS/DESIGN
    MULTIFORM_results                 = MULTIFORM_cost_module(fc["Src1_MULTIFORM"]["P-PO4"])
    #MULTIFORMA_size                   = MULTIFORM_results['FBR_size']
    n_MULTIFORM                       = MULTIFORM_results['n_MULTIFORM']
    MULTIFORM_equipment_cost          = MULTIFORM_results['MULTIFORM_equipment_cost']
    MULTIFORM_operating_cost_partial  = MULTIFORM_results['MULTIFORM_operating_cost_partial']
    
    equipment_cost              = MULTIFORM_equipment_cost
    #    physical_plant_cost_2016    = 3.15*equipment_cost
    #    fixed_capital_cost_2016     = 1.4*physical_plant_cost_2016
    fixed_capital_cost_2016     = equipment_cost #FINAL REPORT Struvite or traditional chemical phosphorus precipitation. What option rocks?
    
    #chemicals_cost_2016 = price["MgCl2"]*fc["Src2_MULTIFORM"]["MgCl2"]*3600*24*334 + NaOH_kg*price["NaOH"]
    
    #operation_cost_2016 = (chemicals_cost_2016+0.3*fixed_capital_cost_2016+FBR_operating_cost_partial) #+1.5*labour_cost_2016
    operation_cost_2016_amortized = ((fixed_capital_cost_2016/ec_param['plant_lifetime'])+MULTIFORM_operating_cost_partial) #+1.5*labour_cost_2016
    operation_cost_2016_non_amortized = (MULTIFORM_operating_cost_partial)
        
    Struvite_benefits = Struvite_recovered*3600*24*334*price['struvite']
        
    #Benefits = Struvite_benefits-operation_cost_2016
    
    recovered_P = fc["DryerMULTIFORM_Peff"]["P"]
    recovered_PO4 = fc["DryerMULTIFORM_Peff"]["Struvite"]/MW['Struvite']*MW['P']+3*fc["DryerMULTIFORM_Peff"]["HAP"]/MW['HAP']*MW['P']
    released_P = fc["MULTIFORM_LiqEff"]["P"]
    released_PO4 = fc["MULTIFORM_LiqEff"]["P-PO4"]+fc["DryerMULTIFORM_Peff"]["P-PO4"]
    released_NH3 = (fc["MULTIFORM_LiqEff"]["N-NH3"]+fc["DryerMULTIFORM_Peff"]["N-NH3"]) #kgN as NH4/s
    released_N = (fc["MULTIFORM_LiqEff"]["N"]) #kg rest of N/s
    
    fraction_recoved_P = recovered_P/(fc["Src1_MULTIFORM"]["P"])
    fraction_recoved_PO4 = recovered_PO4/fc["Src1_MULTIFORM"]["P-PO4"]
    fraction_recoved_TP = (recovered_P+recovered_PO4)/(fc["Src1_MULTIFORM"]["P-PO4"]+fc["Src1_MULTIFORM"]["P"])
    fraction_released_P = released_P/(fc["Src1_MULTIFORM"]["P"])
    fraction_released_PO4 = released_PO4/(fc["Src1_MULTIFORM"]["P-PO4"])
    fraction_released_TP = (released_P+released_PO4)/(fc["Src1_MULTIFORM"]["P-PO4"]+fc["Src1_MULTIFORM"]["P"])
    PO4_conc_released = (fc["MULTIFORM_LiqEff"]["P-PO4"]+fc["DryerMULTIFORM_Peff"]["P-PO4"])/MW['P']*MW['PO4']*UnitConv['K_to_mili']/(F["MULTIFORM_LiqEff"]/feedstock_parameters['Wa_density'])
    Eutrophication_potential = (released_N*0.42+released_NH3/MW['N']*MW['NH4']*0.33+released_PO4/MW['P']*MW['PO4']*1+released_P*3.06)*365*3600*24/1000 #ton phosphate equivalent / year , potency factors from 'The sustainability metrics. Sustainable Development Progess Metrics recommended for use in the Process Idustries' IChenE
    
    return {'fc':fc,
            'x':x,
            'F':F,
            'check_F':checks_F,
            'checks_x':checks_x,
            'tech':'MULTIFORM',
            'fixed_capital_cost_2016':fixed_capital_cost_2016,
            'MULTIFORM_size':'Unique', 'n_MULTIFORM':n_MULTIFORM,
            'Struvite_benefits':Struvite_benefits, 
            'investment_cost':fixed_capital_cost_2016,
            'equipment_cost':MULTIFORM_equipment_cost,
            'operation_cost_2016_amortized':operation_cost_2016_amortized,
            'operation_cost_2016_non_amortized':operation_cost_2016_non_amortized,
            'Struvite_benefits':Struvite_benefits,
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



