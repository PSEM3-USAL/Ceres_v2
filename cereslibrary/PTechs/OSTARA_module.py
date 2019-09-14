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
def OSTARA_module(F_ini,fc_ini,x_ini):
    import pandas as pd
    import numpy as np
    
    # IMPORT MODULES
    from global_parameters_module import UnitConv, MW, c_p_liq_sol, dH_vap_0, Tc, Tb, dH_f, dH_c, c_p_v_1, c_p_v_2, c_p_v_3, c_p_v_4, coef_vapor_pressure_1, coef_vapor_pressure_2, coef_vapor_pressure_3, CEI, price, nu_p, k_p, n_watson, epsilon, T_amb, P_ref, density
    from feedstock_input_module import elements_wet, elements_dry, nutrients, feedstock_parameters, elements_dry_comp
    from equipment_costs.OSTARA_cost_module import OSTARA_cost_module
    from economic_parameters_module import ec_param
    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================
    # SETS, PARAMETERS AND NODES DATA ACQUISITION
    nodes_matrix            = pd.read_csv('nodes/nodes_OSTARA.csv', sep=",", header=None)
    #nodes_matrix            = pd.read_csv('cereslibrary/techmodels/nodes/nodes_FBR.csv', sep=",", header=None)
    nodes       = np.array(nodes_matrix[0].dropna())   
    process_elements_matrix = pd.read_csv('process_elements/process_elements_OSTARA.csv', sep=",", header=0)
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
#        x["Src1Mixer"][i]  = total_elements[i]/100
#        fc["Src1Mixer"][i] = x["Src1Mixer"][i]*F_ini
        x["Src1Mixer"][i]  = x_ini[i]
        fc["Src1Mixer"][i] = fc_ini[i]
        

    fc["Src2Mixer"]["MgCl2"]         = (fc["Src1Mixer"]["P-PO4"]/MW["P"])*Mg_P_ratio*MW["MgCl2"]
    
    #fc["Src3Mixer"]["NaOH"]=0 #TEMPORAL
    
    for i in total_elements.keys():
        fc["MixerFBR"][i] = fc["Src1Mixer"][i] + fc["Src2Mixer"][i] #+ fc["Src3Mixer"][i]
    
    # It is considered that the materia volume is contributed by the digestate. Other contributions are negligible
    # It is considered that volume refered in seeds_amount is the reactor volume
    fc["Src3FBR"]["Struvite_seeds"] = (sum(fc["MixerFBR"][i] for i in total_elements)/feedstock_parameters['digestate_density'])*struvite_seeds_dose
    

    fc["FBRDryer"]["Struvite_seeds"] = fc["Src3FBR"]["Struvite_seeds"]
    # It is consider that it is dry struvite (dry base)
    #fc["FBRSink2"]["Struvite"] = (fc["FBRSink2"]["PO4"]/MW["PO4"])*MW["Struvite"]+fc["FBRSink2"]["Struvite_seeds"]
    
    #Total suspended solids penalization in precipitates formation
    DM_MixerFBR_percentage      = 100*(fc["MixerFBR"]["P"]+fc["MixerFBR"]["N"]+fc["MixerFBR"]["C"]+fc["MixerFBR"]["Rest"]+fc["MixerFBR"]["Ca"]+fc["MixerFBR"]["K"])/sum(fc["MixerFBR"][ii] for ii in total_elements)
    TSS_MixerFBR_percentage     = 0.078*DM_MixerFBR_percentage + 0.051
    penalization_index = 1
    #if TSS_MixerFBR_percentage <= 0.215:
        #penalization_index = 1
    #else:
        #penalization_index = -7.335E-01*TSS_MixerFBR_percentage + 1.158
    
    # Calcium correlations
    Ca_PO4_ratio=(fc["MixerFBR"]["Ca_ion"]/MW["Ca"])/(fc["MixerFBR"]["P-PO4"]/MW["PO4"])
    fc["FBRDryer"]["Struvite"]   = penalization_index*(fc["MixerFBR"]["P-PO4"]*(1-(fines_struvite/100)))*(0.798/(1+(Ca_PO4_ratio*0.576)**2.113))/MW['P']*MW['Struvite']
    fc["FBRDryer"]["CaCO3"]      = penalization_index*(fc["MixerFBR"]["Ca_ion"])*(1.020/(1+((Ca_PO4_ratio)*4.097E-01)**(1.029)))/MW['Ca']*MW['CaCO3']
    fc["FBRDryer"]["HAP"]        = penalization_index*(fc["MixerFBR"]["Ca_ion"])*(-4.321E-02*Ca_PO4_ratio*Ca_PO4_ratio + 3.128E-01*Ca_PO4_ratio - 3.619E-02)/MW['Ca']/5*MW['HAP']
    
    fc["FBRDryer"]["Rest"]   = fc["MixerFBR"]["Rest"]
    fc["FBRDryer"]["C"]      = fc["MixerFBR"]["C"]
    fc["FBRDryer"]["N"]      = fc["MixerFBR"]["N"]
    fc["FBRDryer"]["P"]      = fc["MixerFBR"]["P"]
    fc["FBRDryer"]["Ca"]     = fc["MixerFBR"]["Ca"]
    fc["FBRDryer"]["K"]      = fc["MixerFBR"]["K"]
    
    fc["FBRDryer"]["Wa"] = ((fc["FBRDryer"]["Struvite"]+fc["FBRDryer"]["CaCO3"] +fc["FBRDryer"]["HAP"]+fc["FBRDryer"]["Rest"]+fc["FBRDryer"]["C"] +fc["FBRDryer"]["N"] +fc["FBRDryer"]["P"]+fc["FBRDryer"]["Ca"] +fc["FBRDryer"]["K"])/(1-(humidity_FBR/100)))*(humidity_FBR/100)
    
    fc["FBRHydrocyclone"]["Wa"]     = fc["MixerFBR"]["Wa"]-fc["FBRDryer"]["Wa"]
    Relative_water                  = fc["FBRDryer"]["Wa"]/fc["MixerFBR"]["Wa"]
    
    
    fc["FBRDryer"]["P-PO4"]    = Relative_water*(fc["MixerFBR"]["P-PO4"]/MW['P']-fc["FBRDryer"]["Struvite"]/MW['Struvite']-3*fc["FBRDryer"]["HAP"]/MW['HAP'])*MW['P']
    if fc["FBRDryer"]["P-PO4"] >= 0:
        pass
    else:
        print('ERROR negative fc["FBRDryer"]["P-PO4"] variable')
    
    fc["FBRDryer"]["N-NH3"]    = Relative_water*(fc["MixerFBR"]["N-NH3"]/MW['N']-fc["FBRDryer"]["Struvite"]/MW['Struvite'])*MW['N']
    if fc["FBRDryer"]["N-NH3"] >= 0:
        pass
    else:
        print('ERROR negative fc["FBRDryer"]["N-NH3"] variable')
    
    fc["FBRDryer"]["Mg"]    = Relative_water*(fc["MixerFBR"]["MgCl2"]/MW['MgCl2']-fc["FBRDryer"]["Struvite"]/MW['Struvite'])*MW['Mg']
    if fc["FBRDryer"]["Mg"] >= 0:
        pass
    else:
        print('ERROR negative fc["FBRDryer"]["Mg"] variable')
        
    fc["FBRDryer"]["Cl"]     = Relative_water*(fc["MixerFBR"]["MgCl2"]/MW["MgCl2"])*2*MW["Cl"]
    
    fc["FBRDryer"]["Ca_ion"] = Relative_water*(fc["MixerFBR"]["Ca_ion"]-fc["FBRDryer"]["CaCO3"]/MW['CaCO3']*MW['Ca']-5*fc["FBRDryer"]["HAP"]/MW['HAP']*MW['Ca'])
    if fc["FBRDryer"]["Ca_ion"] >= 0:
        pass
    else:
        print('ERROR negative fc["FBRDryer"]["Ca_ion"] variable')
        
    fc["FBRDryer"]["K_ion"]  = Relative_water*fc["MixerFBR"]["K_ion"]
    
    fc["DryerSink3"]["Wa"]          = fc["FBRDryer"]["Wa"] 
    
    fc["DryerSink4"]["Struvite"]  = fc["FBRDryer"]["Struvite"]
    fc["DryerSink4"]["HAP"]       = fc["FBRDryer"]["HAP"]
    fc["DryerSink4"]["CaCO3"]     = fc["FBRDryer"]["CaCO3"]
    fc["DryerSink4"]["Rest"]      = fc["FBRDryer"]["Rest"]
    fc["DryerSink4"]["C"]         = fc["FBRDryer"]["C"]
    fc["DryerSink4"]["P"]         = fc["FBRDryer"]["P"]
    fc["DryerSink4"]["N"]         = fc["FBRDryer"]["N"]
    fc["DryerSink4"]["Ca"]        = fc["FBRDryer"]["Ca"]
    fc["DryerSink4"]["K"]         = fc["FBRDryer"]["K"]
    fc["DryerSink4"]["P-PO4"]     = fc["FBRDryer"]["P-PO4"]
    fc["DryerSink4"]["N-NH3"]     = fc["FBRDryer"]["N-NH3"]
    fc["DryerSink4"]["Mg"]        = fc["FBRDryer"]["Mg"]
    fc["DryerSink4"]["Cl"]        = fc["FBRDryer"]["Cl"]
    fc["DryerSink4"]["Ca_ion"]    = fc["FBRDryer"]["Ca_ion"]
    fc["DryerSink4"]["K_ion"]     = fc["FBRDryer"]["K_ion"]
    
    
    fc["FBRHydrocyclone"]["P-PO4"]      = fc["MixerFBR"]["P-PO4"]-fc["FBRDryer"]["P-PO4"]-fc["FBRDryer"]["Struvite"]/MW['Struvite']*MW['P']-3*fc["FBRDryer"]["HAP"]/MW['HAP']*MW['P']
    fc["FBRHydrocyclone"]["N-NH3"]      = fc["MixerFBR"]["N-NH3"]-fc["FBRDryer"]["N-NH3"]-fc["FBRDryer"]["Struvite"]/MW['Struvite']*MW['N']
    fc["FBRHydrocyclone"]["Mg"]         = ((fc["MixerFBR"]["MgCl2"]/MW['MgCl2']-fc["FBRDryer"]["Struvite"]/MW['Struvite'])*MW['Mg'])-fc["FBRDryer"]["Mg"]
    fc["FBRHydrocyclone"]["Cl"]         = fc["MixerFBR"]["Cl"]-fc["FBRDryer"]["Cl"]
    fc["FBRHydrocyclone"]["Ca_ion"]     = fc["MixerFBR"]["Ca_ion"]-fc["FBRDryer"]["Ca_ion"]-5*fc["FBRDryer"]["HAP"]/MW['HAP']*MW['Ca']-fc["FBRDryer"]["CaCO3"]/MW['CaCO3']*MW['Ca']
    fc["FBRHydrocyclone"]["K_ion"]      = fc["MixerFBR"]["K_ion"]-fc["FBRDryer"]["K_ion"]
    
    fc["FBRHydrocyclone"]["Rest"]       = fc["MixerFBR"]["Rest"]-fc["FBRDryer"]["Rest"]
    fc["FBRHydrocyclone"]["C"]          = fc["MixerFBR"]["C"]-fc["FBRDryer"]["C"]
    fc["FBRHydrocyclone"]["K"]          = fc["MixerFBR"]["K"]-fc["FBRDryer"]["K"]
    fc["FBRHydrocyclone"]["N"]          = fc["MixerFBR"]["N"]-fc["FBRDryer"]["N"]
    fc["FBRHydrocyclone"]["P"]          = fc["MixerFBR"]["P"]-fc["FBRDryer"]["P"]
    fc["FBRHydrocyclone"]["Ca"]         = fc["MixerFBR"]["Ca"]-fc["FBRDryer"]["Ca"]
    fc["FBRHydrocyclone"]["Struvite"]   = fc["FBRDryer"]["Struvite"]*(fines_struvite/100)
    fc["FBRHydrocyclone"]["HAP"]        = fc["FBRDryer"]["HAP"]*(fines_struvite/100)
    fc["FBRHydrocyclone"]["CaCO3"]      = fc["FBRDryer"]["CaCO3"]*(fines_struvite/100)
    
    
    
    #fc["FBRHydrocyclone"]["Cl"] = 2*(fc["MixerFBR"]["MgCl2"]/MW["MgCl2"])*MW["Cl"]
    #fc["FBRHydrocyclone"]["Mg"] = (fc["Src2Mixer"]["MgCl2"]/MW["MgCl2"])*MW["Mg"]-fc["FBRSink2"]["Mg"]
    #fc["FBRHydrocyclone"]["Rest"] = fc["MixerFBR"]["Rest"]
    #fc["FBRHydrocyclone"]["C"] = fc["MixerFBR"]["C"]
    #fc["FBRHydrocyclone"]["K"] = fc["MixerFBR"]["K"]
    #fc["FBRHydrocyclone"]["K_ion"] = fc["MixerFBR"]["K_ion"]
    #fc["FBRHydrocyclone"]["N"] = fc["MixerFBR"]["N"]
    #fc["FBRHydrocyclone"]["P"] = fc["MixerFBR"]["P"]
    #fc["FBRHydrocyclone"]["Ca"] = fc["MixerFBR"]["Ca"]
    #fc["FBRHydrocyclone"]["Ca_ion"] = fc["MixerFBR"]["Ca_ion"]
    #fc["FBRHydrocyclone"]["NaOH"] = fc["MixerFBR"]["NaOH"]
    #fc["FBRHydrocyclone"]["Wa"] = fc["MixerFBR"]["Wa"] - fc["FBRSink2"]["Wa"] 
    #fc["FBRHydrocyclone"]["NH3"] = fc["MixerFBR"]["NH3"] - fc["FBRSink2"]["NH3"] 
    #fc["FBRHydrocyclone"]["PO4"] = fc["MixerFBR"]["PO4"] - fc["FBRSink2"]["PO4"] 
    #fc["FBRHydrocyclone"]["Struvite"] = fc["FBRDryer"]["Struvite"]*(fines_struvite/100)
    
    for i in total_elements.keys():
        #if i!="Struvite":
        fc["HydrocycloneSink2"][i] =fc["FBRHydrocyclone"][i]
    
    #fc["HydrocycloneSink1"]["Struvite"] = fc["FBRHydrocyclone"]["Struvite"]
    
    Struvite_recovered = fc["FBRDryer"]["Struvite"]#+fc["FBRDryer"]["HAP"]# + fc["HydrocycloneSink1"]["Struvite"]
    
    #Struvite_benefits = Struvite_recovered*3600*24*334*price_struvite
    
    
    for i in nodes_list:
        F[i] = sum(fc[i][ii] for ii in total_elements)
        
    for i in nodes_list:
        if i!="Src1Mixer":
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
    
    if abs(F["Src1Mixer"] + F["Src2Mixer"] + F["Src3FBR"] - (F["HydrocycloneSink2"] + F["DryerSink3"] + F["DryerSink4"])) <= 0.005: #F["HydrocycloneSink1"] +
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
    OSTARA_results                 = OSTARA_cost_module(fc["FBRDryer"]["Struvite"])
    OSTARA_size                    = OSTARA_results['FBR_size']
    n_OSTARA                       = OSTARA_results['n_FBR']
    OSTARA_equipment_cost          = OSTARA_results['FBR_equipment_cost']
    OSTARA_operating_cost_partial  = OSTARA_results['FBR_operating_cost_partial']
    
    equipment_cost              = OSTARA_equipment_cost
    #    physical_plant_cost_2016    = 3.15*equipment_cost
    #    fixed_capital_cost_2016     = 1.4*physical_plant_cost_2016
    fixed_capital_cost_2016     = equipment_cost #1.9*equipment_cost #FINAL REPORT Struvite or traditional chemical phosphorus precipitation. What option rocks?
    
    #chemicals_cost_2016 = price["MgCl2"]*fc["Src2Mixer"]["MgCl2"]*3600*24*334 + NaOH_kg*price["NaOH"]
    
    #operation_cost_2016 = (chemicals_cost_2016+0.3*fixed_capital_cost_2016+FBR_operating_cost_partial) #+1.5*labour_cost_2016
    operation_cost_2016_amortized = ((fixed_capital_cost_2016/ec_param['plant_lifetime'])+FBR_operating_cost_partial) #+1.5*labour_cost_2016
    operation_cost_2016_non_amortized = (FBR_operating_cost_partial)
        
    Struvite_benefits = Struvite_recovered*3600*24*334*price['struvite']
        
    #Benefits = Struvite_benefits-operation_cost_2016
    
    recovered_P = fc["DryerSink4"]["P"]
    recovered_PO4 = fc["DryerSink4"]["Struvite"]/MW['Struvite']*MW['P']+3*fc["DryerSink4"]["HAP"]/MW['HAP']*MW['P']
    released_P = fc["HydrocycloneSink2"]["P"]
    released_PO4 = fc["HydrocycloneSink2"]["P-PO4"]+fc["DryerSink4"]["P-PO4"]
    
    fraction_recoved_P = recovered_P/(fc["Src1Mixer"]["P"])
    fraction_recoved_PO4 = recovered_PO4/fc["Src1Mixer"]["P-PO4"]
    fraction_recoved_TP = (recovered_P+recovered_PO4)/(fc["Src1Mixer"]["P-PO4"]+fc["Src1Mixer"]["P"])
    fraction_released_P = released_P/(fc["Src1Mixer"]["P"])
    fraction_released_PO4 = released_PO4/(fc["Src1Mixer"]["P-PO4"])
    fraction_released_TP = (released_P+released_PO4)/(fc["Src1Mixer"]["P-PO4"]+fc["Src1Mixer"]["P"])
    PO4_conc_released = (fc["HydrocycloneSink2"]["P-PO4"]+fc["DryerSink4"]["P-PO4"])*UnitConv['K_to_mili']/(F["HydrocycloneSink2"]/feedstock_parameters['Wa_density'])
    #(sum(fc["HydrocycloneSink2"][ii] for ii in total_elements)/feedstock_parameters['Wa_density'])
    
    return {'fc':fc,
            'x':x,
            'F':F,
            'check_F':checks_F,
            'checks_x':checks_x,
            'tech':'OSTARA',
            'fixed_capital_cost_2016':fixed_capital_cost_2016,
            'OSTARA_size':OSTARA_size, 'n_OSTARA':n_OSTARA,
            'Struvite_benefits':Struvite_benefits, 
            'investment_cost':fixed_capital_cost_2016,
            'equipment_cost':FBR_equipment_cost,
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



