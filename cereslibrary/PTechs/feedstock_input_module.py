#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 09:15:46 2018

@author: emh
"""

#feedstock_input_module
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
def feedstock_input_module (N_animals_DairyCow, N_animals_DairyHeifer, N_animals_DairyCalf, N_animals_BeefCow, N_animals_BeefCalf):
    
    import numpy as np
    import pandas as pd
        
    pd.options.display.max_columns = 50

    # MANURE DATA ACQUISITION AND COMPOSITION SETTLEMENT
    #waste_data          = pd.read_csv('feedstock_data/waste_data_input.csv', sep=",", header=0)
    ##waste_data          = pd.read_csv('feedstock_data/waste_data_input.csv', sep=",", header=0)
    #input_comp_raw_comp = np.array(waste_data["Component"].dropna()) #NaN removed
    #input_comp_raw_conc = np.array(waste_data["Concentration"].dropna()) #NaN removed

    #input_comp_raw  = dict(zip(input_comp_raw_comp, input_comp_raw_conc))
    #comp_ready      = ["Wa", "OM", "C", "Rest"]
    #input_comp      = {key:input_comp_raw[key] for key in comp_ready}
    #input_comp.update({"N":input_comp_raw["N"]*(1-input_comp_raw["NH3_N_ratio"]), 
                    #"P":input_comp_raw["P"]*(1-input_comp_raw["PO4_P_ratio"]),
                    #"Ca":input_comp_raw["Ca"]*(1-input_comp_raw["Caion_Ca_ratio"]), 
                    #"K":input_comp_raw["K"]*(1-input_comp_raw["Kion_K_ratio"]),
                    #"N-NH3":input_comp_raw["N"]*input_comp_raw["NH3_N_ratio"], 
                    #"P-PO4":input_comp_raw["P"]*input_comp_raw["PO4_P_ratio"], 
                    #"Ca_ion":input_comp_raw["Ca"]*input_comp_raw["Caion_Ca_ratio"], 
                    #"K_ion":input_comp_raw["K"]*input_comp_raw["Kion_K_ratio"]})

    #elements_wet_comp   = ["Wa", "C", "N", "P", "Ca", "K", "N-NH3", "P-PO4", "Ca_ion", "K_ion", "Rest"]
    #elements_dry_comp   = ["C", "N", "P", "Ca", "K", "N-NH3", "P-PO4", "Ca_ion", "K_ion", "Rest"]
    #nutrients_comp      = ["N-NH3", "P-PO4"]
    #elements_wet        = {key:input_comp[key] for key in elements_wet_comp} #SUBSET wet digestate elements
    #elements_dry        = {key:input_comp[key] for key in elements_dry_comp} #SUBSET dry digestate elements
    #nutrients           = {key:input_comp[key] for key in nutrients_comp} #SUBSET nutrients

    livestock_list = ['Dairy Cow','Dairy Heifer','Dairy Calf','Beef Cow','Beef Calf']
    waste_data          = pd.read_csv('feedstock_data/waste_data_input.csv', sep=",", header=0)
    waste_data = waste_data.set_index('Livestock type', drop=False)

    waste_data['N'] = waste_data['N_tot']*(1-waste_data["NH3_N_ratio"])
    waste_data['P'] = waste_data['P_tot']*(1-waste_data["PO4_P_ratio"])
    waste_data['Ca'] = waste_data['Ca_tot']*(1-waste_data["Caion_Ca_ratio"])
    waste_data['K'] = waste_data['K_tot']*(1-waste_data["Kion_K_ratio"])
    waste_data['N-NH3'] = waste_data['N_tot']*(waste_data["NH3_N_ratio"])
    waste_data['P-PO4'] = waste_data['P_tot']*(waste_data["PO4_P_ratio"])
    waste_data['Ca_ion'] = waste_data['Ca_tot']*(waste_data["Caion_Ca_ratio"])
    waste_data['K_ion'] = waste_data['K_tot']*(waste_data["Kion_K_ratio"])
    waste_data['C'] = waste_data['OM']*(waste_data["C_OM_ratio"])

    F_animals_DairyCow = N_animals_DairyCow/waste_data.loc['Dairy Cow','AU_equivalence']*waste_data.loc['Dairy Cow','Weight']
    F_animals_DairyHeifer = N_animals_DairyHeifer/waste_data.loc['Dairy Heifer','AU_equivalence']*waste_data.loc['Dairy Cow','Weight']
    F_animals_DairyCalf = N_animals_DairyCalf/waste_data.loc['Dairy Calf','AU_equivalence']*waste_data.loc['Dairy Calf','Weight']
    F_animals_BeefCow = N_animals_BeefCow/waste_data.loc['Beef Cow','AU_equivalence']*waste_data.loc['Beef Cow','Weight']
    F_animals_BeefCalf = N_animals_BeefCalf/waste_data.loc['Beef Calf','AU_equivalence']*waste_data.loc['Beef Calf','Weight']

    Total_animals = F_animals_DairyCow+F_animals_DairyHeifer+F_animals_DairyCalf+F_animals_BeefCow+F_animals_BeefCalf

    Manure_inlet = pd.DataFrame({'N_animals_DairyCow':[N_animals_DairyCow],
                                'N_animals_DairyHeifer':[N_animals_DairyHeifer],
                                'N_animals_DairyCalf':[N_animals_DairyCalf],
                                'N_animals_BeefCow':[N_animals_BeefCow],
                                'N_animals_BeefCalf':[N_animals_BeefCalf],
                                'Dairy Cow':[F_animals_DairyCow],
                                'Dairy Heifer':[F_animals_DairyHeifer],
                                'Dairy Calf':[F_animals_DairyCalf],
                                'Beef Cow':[F_animals_BeefCow],
                                'Beef Calf':[F_animals_BeefCalf],
                                'Total_animals':[Total_animals],
                                })

    elements_wet_comp   = ["Wa", 'C', "N", "P", "Ca", "K", "N-NH3", "P-PO4", "Ca_ion", "K_ion", "Rest"]
    elements_dry_comp   = ["C", "N", "P", "Ca", "K", "N-NH3", "P-PO4", "Ca_ion", "K_ion", "Rest"]
    nutrients_comp      = ["N-NH3", "P-PO4"]

    elements_wet = dict()
    for i in elements_wet_comp:
        elements_wet[i] = sum(Manure_inlet[ii].item()/Manure_inlet['Total_animals'].item()*waste_data.loc[ii,i].item() for ii in livestock_list)

    elements_dry = dict()
    for i in elements_dry_comp:
        elements_dry[i] = sum(Manure_inlet[ii].item()/Manure_inlet['Total_animals'].item()*waste_data.loc[ii,i].item() for ii in livestock_list)
        
    nutrients = dict()
    for i in nutrients_comp:
        nutrients[i] = sum(Manure_inlet[ii].item()/Manure_inlet['Total_animals'].item()*waste_data.loc[ii,i].item() for ii in livestock_list)

    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================
    # FEDSTOCK PARAMETERS
    #feedstock_parameters_matrix          = pd.read_csv('feedstock_data/feedstock_parameters.csv', sep=",", header=0)
    feedstock_parameters_matrix          = pd.read_csv('feedstock_data/feedstock_parameters.csv', sep=",", header=0)
    feedstock_parameters                 = dict(zip(np.array(feedstock_parameters_matrix["Name"].dropna()), np.array(feedstock_parameters_matrix["Value"].dropna())))
    
    return {'elements_wet':elements_wet,
            'elements_dry':elements_dry,
            'nutrients':nutrients,
            'feedstock_parameters':feedstock_parameters,
            'elements_dry_comp':elements_dry_comp,
            'Total_animals':Total_animals,
            }
