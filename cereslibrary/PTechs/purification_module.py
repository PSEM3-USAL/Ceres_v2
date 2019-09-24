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
def purification_module(F_ini,fc_ini,x_ini, elements_wet, elements_dry, nutrients, feedstock_parameters):
    import numpy as np
    import pandas as pd

    # SETS, PARAMETERS AND NODES DATA ACQUISITION
    nodes_matrix_purif            = pd.read_csv('nodes/nodes_purif.csv', sep=",", header=None)
    nodes_purif                   = np.array(nodes_matrix_purif[0].dropna())
    process_elements_matrix_purif = pd.read_csv('process_elements/process_elements_digestor.csv', sep=",", header=0)

    from global_parameters_module import UnitConv, MW, c_p_liq_sol, dH_vap_0, Tc, Tb, dH_f, dH_c, c_p_v_1, c_p_v_2, c_p_v_3, c_p_v_4, coef_vapor_pressure_1, coef_vapor_pressure_2, coef_vapor_pressure_3, CEI, price, nu_p, k_p, n_watson, epsilon, T_amb, P_ref, density, latent_heat_evap, nat_gas_heat_value
    #from feedstock_input_module import elements_wet, elements_dry, nutrients, feedstock_parameters
    from economic_parameters_module import ec_param


    # MANURE DATA ACQUISITION AND COMPOSITION SETTLEMENT   
    gases_comp_purif   = np.array(process_elements_matrix_purif["Component"].dropna())
    gases_conc_purif   = np.full((len(gases_comp_purif)), 0)
    gases_purif        = dict(zip(gases_comp_purif, gases_conc_purif))

    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================  
    # TOTAL ELEMENTS
    total_elements_purif = {**elements_wet,**gases_purif} # Merge the two dictionaries

    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================
    # AD PARTICULAR PARAMETERS
    T_biogas_out = 55
    rendcom = 0.85 #*Rendimiento compresores/0.85/
    Pcomp1 = 2.5 # *Presión de salida del compresor para el reactor de eliminación del H2S /2.5/
    R = 8.314 #*Constante de los gases: kJ/kmol·K
    H_Wa = 2254.62 #*Calor estándar de vaporización: kJ/kg (100ºC)
    Tc_Wa = 374.15 #*Temperatura crítica del agua: ºC
    Tb_Wa = 100 #*Temperatura punto de ebullición del agua: ºC
    Nwatson = 0.38 #*Número de watson
    Pcomp2 = 4.5 #*Presión de salida del compresor 2:
    RelacPcomp2 = 1.1 #*Relación de presiones Pout/Pin, debe encontrarse en el reactor 3 a 4,5 bar, y aunquese alcanza esa presión en el anterior, hay pérdida presión, por lo que va a baja*se "supone" que hay una relación entre ellas de 1,1:

        
    # VARIABLES DEFINITION (IN NESTED DICTIONARIES) (INITIALIZATION)
    nodes_list_purif              = nodes_purif.tolist()
    initialization_comp_purif     = total_elements_purif #["Wa", "C", "NH3", "PO4", "Ca_ion", "K_ion"]
    initialization_nan_purif      = np.full((len(initialization_comp_purif)), 0.00)

    fc_purif = {key: dict(zip(initialization_comp_purif,initialization_nan_purif)) for key in nodes_list_purif}

    x_purif = {key: dict(zip(initialization_comp_purif,initialization_nan_purif)) for key in nodes_list_purif}

    y_purif  = {key: dict(zip(initialization_comp_purif,initialization_nan_purif)) for key in nodes_list_purif}

    F_purif = {key: np.nan for key in nodes_list_purif}

    T_purif = {key: np.nan for key in nodes_list_purif}


    # ===============================================================================================================================================================================================================================
    # ###############################################################################################################################################################################################################################
    # ===============================================================================================================================================================================================================================     
    # MASS BALANCE
    for i in total_elements_purif.keys():
        x_purif["Bioreactor_Comp1"][i]  = x_ini[i]
        fc_purif["Bioreactor_Comp1"][i] = fc_ini[i]
        #x_purif["Bioreactor_Comp1"][i]  = x["ADEval_Biogas"][i]
        #fc_purif["Bioreactor_Comp1"][i] = fc_kgs["ADEval_Biogas"][i]
        
    for i in total_elements_purif.keys():
    #     x_purif["Comp1_HX1"][i]  = x_purif["Bioreactor_Comp1"][i]
        fc_purif["Comp1_HX1"][i] = fc_purif["Bioreactor_Comp1"][i]
        
    for i in total_elements_purif.keys():
    #     x_purif["HX1_Sep1"][i]  = x_purif["Comp1_HX1"][i]
        fc_purif["HX1_Sep1"][i] = fc_purif["Comp1_HX1"][i]
        
    T_purif["Comp1_HX1"] = ((T_biogas_out+273)+(1/rendcom)*(T_biogas_out+273)*((Pcomp1+0.001)**(0.4/1.4)-1))-273

    # *Paso la fracción masica a fracción molar y obtengo el valor del nuevo peso molecular
    for i in gases_purif:
        y_purif["Bioreactor_Comp1"][i] = fc_purif["Bioreactor_Comp1"][i]/MW[i]/(sum(fc_purif["Bioreactor_Comp1"][ii]/MW[ii] for ii in gases_purif)+ fc_purif["Bioreactor_Comp1"]['Wa']/MW['Wa'])

    y_purif["Bioreactor_Comp1"]['Wa'] = fc_purif["Bioreactor_Comp1"]['Wa']/MW['Wa']/(sum(fc_purif["Bioreactor_Comp1"][ii]/MW[ii] for ii in gases_purif)+ fc_purif["Bioreactor_Comp1"]['Wa']/MW['Wa'])

    # *Obtengo el peso molecular del biogas con el agua
    MW_biogas_comp1 = sum(y_purif["Bioreactor_Comp1"][ii]*MW[ii] for ii in gases_purif)+y_purif["Bioreactor_Comp1"]['Wa']*MW['Wa']

    # *Trabajo del compresor
    # *Aquí la temperatura debe estar en Kelvin porque R está en K
    W_Comp1 = ((1/rendcom)*(R*(T_biogas_out+273)*1.4*sum(fc_purif["Bioreactor_Comp1"][ii] for ii in total_elements_purif))*(((Pcomp1)**(0.4/1.4))-1))/(0.4*(MW_biogas_comp1+0.001))

    # *Una vez salen los gases del compresor pasan al intercambiador de calor para reducir su temperatura
    # *hasta la de trabajo del reactor de eliminación de H2S: 25ºC
    T_purif["HX1_Sep1"]=25

    #Se calcula la cantidad de agua que condensa: ((Pv/(P-Pv))*(PMvapor/PMgas))=(fcvapor/fcgas) vapor:agua, gas:biogas
    # **Presión de saturación del agua(gas) que sale del IC11
    # Purif8..  Psat_HX1=E=10**(coef_An('Wa','1')-(coef_An('Wa','2')/(eps1+T('IC11','Sep7')+coef_An('Wa','3'))));
    Psat_HX1 = 10**(coef_vapor_pressure_1['Wa']-(coef_vapor_pressure_2['Wa']/(T_purif["HX1_Sep1"]+coef_vapor_pressure_3['Wa'])))

    # *Cálculo del nuevo peso molecular del biogas*
    for i in gases_purif:
        y_purif["HX1_Sep1"][i] = fc_purif["HX1_Sep1"][i]/MW[i]/(sum(fc_purif["HX1_Sep1"][ii]/MW[ii] for ii in gases_purif))
    MW_biogas_HX1 = sum(y_purif["HX1_Sep1"][ii]*MW[ii] for ii in gases_purif)

    # *Cálculo de la humedad:
    Wa_HX1 = MW['Wa']*(Psat_HX1/750)/(MW_biogas_HX1*(Pcomp1-(Psat_HX1/750)))

    # *Con la humedad que se obtenga, se tiene el valor del flujo del vapo de agua en esa corriente (IC11->Sep7)
    Wa_nocondensada1 = Wa_HX1*sum(fc_purif["HX1_Sep1"][ii] for ii in total_elements_purif)
                                                                                                                            
    # *Se añade ahora la ecuación para el calor eliminado del gas:
    # *(Puesto que no está el set de biogasWa en ninguna variable fuera de los sumatorios, no se especifica)
    # *¡¡¡¡CALOR SENSIBLE Y CALOR LATENTE!!!!!(tener en cuenta la regla de Watson)


    Q_HX1 = sum(fc_purif["Comp1_HX1"][ii]/MW[ii]*(c_p_v_1[ii]*(T_purif["HX1_Sep1"]-T_purif["Comp1_HX1"] )+
                                                    1/2*c_p_v_2[ii]*((T_purif["HX1_Sep1"]+273)**2-(T_purif["Comp1_HX1"]+273)**2 )+
                                                    1/3*c_p_v_3[ii]*((T_purif["HX1_Sep1"]+273)**3-(T_purif["Comp1_HX1"]+273)**3 )+
                                                    1/4*c_p_v_4[ii]*((T_purif["HX1_Sep1"]+273)**4-(T_purif["Comp1_HX1"]+273)**4)) for ii in gases_purif)+fc_purif["Comp1_HX1"]['Wa']/MW['Wa']*(c_p_v_1['Wa']*(T_purif["HX1_Sep1"]-T_purif["Comp1_HX1"] )+
                                                    1/2*c_p_v_2['Wa']*((T_purif["HX1_Sep1"]+273)**2-(T_purif["Comp1_HX1"]+273)**2 )+
                                                    1/3*c_p_v_3['Wa']*((T_purif["HX1_Sep1"]+273)**3-(T_purif["Comp1_HX1"]+273)**3 )+
                                                    1/4*c_p_v_4['Wa']*((T_purif["HX1_Sep1"]+273)**4-(T_purif["Comp1_HX1"]+273)**4))-((fc_purif["Comp1_HX1"]['Wa']-Wa_nocondensada1)*H_Wa*((Tc_Wa-T_purif["HX1_Sep1"])/((Tc_Wa-Tb_Wa)**Nwatson)))

    for i in gases_purif.keys():
        x_purif["Sep1_React1"][i]  = x_purif["HX1_Sep1"][i]
        fc_purif["Sep1_React1"][i] = fc_purif["HX1_Sep1"][i]
        
    fc_purif["Sep1_React1"]['Wa'] = Wa_nocondensada1
    fc_purif["Sep1_Sink1Wa"]['Wa'] = fc_purif["HX1_Sep1"]['Wa']-Wa_nocondensada1

    # *Se fija que la temperatura se debe manterner
    T_purif["Sep1_React1"]=25

    # *Balances al reactor 2( se elimina por completo el H2S y se forma agua)
    fc_purif["React1_Comp2"]['CH4'] = fc_purif["Sep1_React1"]['CH4']
    fc_purif["React1_Comp2"]['CO2'] = fc_purif["Sep1_React1"]['CO2']
    fc_purif["React1_Comp2"]['O2'] = fc_purif["Sep1_React1"]['O2']
    fc_purif["React1_Comp2"]['N2'] = fc_purif["Sep1_React1"]['N2']
    fc_purif["React1_Comp2"]['NH3'] = fc_purif["Sep1_React1"]['NH3']
    fc_purif["React1_Comp2"]['H2S'] = 0
    fc_purif["React1_Comp2"]['Wa'] = fc_purif["Sep1_React1"]['Wa']
    fc_purif["React1_Sink2H2S"]['H2S'] = fc_purif["Sep1_React1"]['H2S']

    for i in total_elements_purif.keys():
    #     x_purif["Comp2_HX2"][i]  = x_purif["React1_Comp2"][i]
        fc_purif["Comp2_HX2"][i] = fc_purif["React1_Comp2"][i]

    T_purif["React1_Comp2"]=25
    T_purif["Comp2_HX2"] = ((T_purif["React1_Comp2"]+273)+((T_purif["React1_Comp2"]+273)*(((RelacPcomp2 )**(0.4/1.4))-1))*(1/rendcom))-273

    # *Paso la fracción masica a fracción molar y obtengo el valor del nuevo peso molecular
    for i in gases_purif:
        y_purif["React1_Comp2"][i] = fc_purif["React1_Comp2"][i]/MW[i]/(sum(fc_purif["React1_Comp2"][ii]/MW[ii] for ii in gases_purif)+ fc_purif["React1_Comp2"]['Wa']/MW['Wa'])

    y_purif["React1_Comp2"]['Wa'] = fc_purif["React1_Comp2"]['Wa']/MW['Wa']/(sum(fc_purif["React1_Comp2"][ii]/MW[ii] for ii in gases_purif)+ fc_purif["React1_Comp2"]['Wa']/MW['Wa'])

    # *Obtengo el peso molecular del biogas con el agua
    MW_biogas_comp2 = sum(y_purif["React1_Comp2"][ii]*MW[ii] for ii in gases_purif)+y_purif["React1_Comp2"]['Wa']*MW['Wa']

    # *Trabajo del compresor
    W_Comp2 = ((1/rendcom)*(R*(T_purif["Comp2_HX2"]+273)*1.4*sum(fc_purif["React1_Comp2"][ii] for ii in gases_purif))*(((RelacPcomp2)**(0.4/1.4))-1))/(0.4*(MW_biogas_comp2))

    # *Flujo total de los componentes a la salida del reactor 2:
    F_purif["React1_Comp2"] = sum(fc_purif["React1_Comp2"][ii] for ii in gases_purif) + fc_purif["React1_Comp2"]['Wa']

    for i in total_elements_purif.keys():
    #     x_purif["Comp2_HX2"][i]  = x_purif["React1_Comp2"][i]
        fc_purif["HX2_Sep2"][i] = fc_purif["Comp2_HX2"][i]

    # *Hasta aquí los pasos del compresor
    # *Una vez salen los gases del compresor pasan al intercambiador de calor para reducir su temperatura
    # *hasta la de trabajo del reactor de eliminación de NH3, H2O y CO2(los dos primeros se eliminan por completo): 25ºC
    T_purif["HX2_Sep2"]=25

    # *Ahora se calcula la cantidad de agua que condensa: ((Pv/(P-Pv))*(PMvapor/PMgas))=(fcvapor/fcgas) vapor:agua, gas:biogas
    # **Presión de saturación del agua(gas) que sale del IC11
    Psat_HX2 = 10**(coef_vapor_pressure_1['Wa']-(coef_vapor_pressure_2['Wa']/(T_purif["HX2_Sep2"]+coef_vapor_pressure_3['Wa'])))

    # *Cálculo del nuevo peso molecular del biogas*
    for i in gases_purif:
        y_purif["HX2_Sep2"][i] = fc_purif["HX2_Sep2"][i]/MW[i]/(sum(fc_purif["HX2_Sep2"][ii]/MW[ii] for ii in gases_purif))
    MW_biogas_HX2 = sum(y_purif["HX2_Sep2"][ii]*MW[ii] for ii in gases_purif)

    # *Cálculo de la humedad:
    Wa_HX2 = MW['Wa']*(Psat_HX2/750)/(MW_biogas_HX2*(Pcomp2-(Psat_HX2/750)))

    # *Con la humedad que se obtenga, se tiene el valor del flujo del vapo de agua en esa corriente (IC11->Sep7)
    Wa_nocondensada2 = Wa_HX2*sum(fc_purif["HX2_Sep2"][ii] for ii in gases_purif)

    Q_HX2 = sum(fc_purif["Comp2_HX2"][ii]/MW[ii]*(c_p_v_1[ii]*(T_purif["HX2_Sep2"]-T_purif["Comp2_HX2"] )+
                                                    1/2*c_p_v_2[ii]*((T_purif["HX2_Sep2"]+273)**2-(T_purif["Comp2_HX2"]+273)**2 )+
                                                    1/3*c_p_v_3[ii]*((T_purif["HX2_Sep2"]+273)**3-(T_purif["Comp2_HX2"]+273)**3 )+
                                                    1/4*c_p_v_4[ii]*((T_purif["HX2_Sep2"]+273)**4-(T_purif["Comp2_HX2"]+273)**4)) for ii in gases_purif)+fc_purif["Comp2_HX2"]['Wa']/MW['Wa']*(c_p_v_1['Wa']*(T_purif["HX2_Sep2"]-T_purif["Comp2_HX2"] )+
                                                    1/2*c_p_v_2['Wa']*((T_purif["HX2_Sep2"]+273)**2-(T_purif["Comp2_HX2"]+273)**2 )+
                                                    1/3*c_p_v_3['Wa']*((T_purif["HX2_Sep2"]+273)**3-(T_purif["Comp2_HX2"]+273)**3 )+
                                                    1/4*c_p_v_4['Wa']*((T_purif["HX2_Sep2"]+273)**4-(T_purif["Comp2_HX2"]+273)**4))-((fc_purif["Comp2_HX2"]['Wa']-Wa_nocondensada2)*H_Wa*((Tc_Wa-T_purif["HX2_Sep2"])/((Tc_Wa-Tb_Wa)**Nwatson)))

    fc_purif["Sep2_Sink3Wa"]['Wa'] = fc_purif["HX2_Sep2"]['Wa']-Wa_nocondensada2
    fc_purif["Sep2_React2"]['Wa'] = Wa_nocondensada2

    for i in total_elements_purif.keys():
    #     x_purif["Comp2_HX2"][i]  = x_purif["React1_Comp2"][i]
        fc_purif["Sep2_React2"][i] = fc_purif["HX2_Sep2"][i]
        
    T_purif["Sep2_React2"] = 25

    # *Ya tengo todos los compuestos que entran al reactor3
    fc_purif["React2_Sink4Biogas"]['CH4'] = fc_purif["Sep2_React2"]['CH4']
    fc_purif["React2_Sink4Biogas"]['CO2'] = fc_purif["Sep2_React2"]['CO2']
    fc_purif["React2_Sink5exhaust"]['O2'] = fc_purif["Sep2_React2"]['O2']
    fc_purif["React2_Sink5exhaust"]['N2'] = fc_purif["Sep2_React2"]['N2']
    fc_purif["React2_Sink5exhaust"]['NH3'] = fc_purif["Sep2_React2"]['NH3']
    fc_purif["React2_Sink5exhaust"]['Wa'] = fc_purif["Sep2_React2"]['Wa']

    T_purif["React2_Sink4Biogas"] = 25
    T_purif["React2_Sink5exhaust"] = 25
    
    
    for i in nodes_list_purif:
        F_purif[i] = sum(fc_purif[i][ii] for ii in total_elements_purif)
        
    for i in nodes_list_purif:
        if i!="Bioreactor_Comp1":
            for ii in total_elements_purif.keys():
                x_purif[i][ii] = fc_purif[i][ii]/F_purif[i]
                
    #Checks
    checks_store = ['OK', 'FAIL']
    checks_F = []
    checks_x = []
    
    #    for i in total_elements.keys():
    #        if abs(fc["Src1_MULTIFORM"][i] + fc["Src2Mixer"][i] + fc["Src3FBR"][i] - (fc["HydrocycloneSink1"][i] + fc["MULTIFORM_LiqEff"][i] + fc["DryerMULTIFORM_VaporEff"][i] + fc["DryerMULTIFORM_Peff"][i])) <= 0.005:
    #            print("fc check", i, "OK")
    #        else:
    #            print("fc check", i, "FAIL")
    
    if abs(F_purif["Bioreactor_Comp1"] - (F_purif["Sep1_Sink1Wa"] + F_purif["React1_Sink2H2S"] + F_purif["Sep2_Sink3Wa"] + F_purif["React2_Sink4Biogas"] + F_purif["React2_Sink5exhaust"])) <= 0.005: #F["HydrocycloneSink1"] +
    #        print("F check OK")
        checks_F = checks_store[0]
    else:
    #        print("F check FAIL")
        checks_F = checks_store[1]
    
    for i in nodes_list_purif:
        if abs(1-sum(x_purif[i][ii] for ii in total_elements_purif)) <= 0.005:
    #            print("x check", i, "OK")
            checks_x = checks_store[0]
        else:
    #            print("x check", i, "FAIL")
            checks_x = checks_store[1]
    print("\n\n")

    #ECONOMIC
    Pelect = 0.06 #           por kwh      /0.06/
    #Molecular sieves
    Molecular_Sieve_HBC_time = 1800 #s
    Molecular_Sieve_HBC_volume = (1.1*(1.077)*Molecular_Sieve_HBC_time/(45*0.454*(1/0.3048**3)))  #m3
    Molecular_Sieve_HBC_diameter = (Molecular_Sieve_HBC_volume*6/(7*3.14))**(1/3)  #m
    Molecular_Sieve_HBC_widht = 0.0023+0.003*Molecular_Sieve_HBC_diameter #m
    Molecular_Sieve_HBC_weight = 8039*(3.14*((Molecular_Sieve_HBC_diameter/2+Molecular_Sieve_HBC_widht)**2-(Molecular_Sieve_HBC_diameter/2)**2)*4*Molecular_Sieve_HBC_diameter+(4/3)*3.14*((Molecular_Sieve_HBC_diameter/2+Molecular_Sieve_HBC_widht)**3-(Molecular_Sieve_HBC_diameter/2)**3)) #kg
    n_Molecular_Sieve_HBC = 4 #total
    Molecular_Sieve_HBC_4Asieve = (1.219*Molecular_Sieve_HBC_time*1.1*1)
    Molecular_Sieve_HBC_cost = n_Molecular_Sieve_HBC*(209.04*Molecular_Sieve_HBC_weight**0.7163+Molecular_Sieve_HBC_4Asieve)

    #Compresors
    Cost_Comp1 = 335.27*W_Comp1+36211
    Cost_Comp2 = 335.27*W_Comp2+36211

    #HXs
    T_in_HOT_HX1 = T_purif["Comp1_HX1"]+273
    T_out_HOT_HX1 = T_purif["HX1_Sep1"]+273
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

    T_in_HOT_HX2 = T_purif["Comp2_HX2"]+273
    T_out_HOT_HX2 = T_purif["HX2_Sep2"]+273
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
        
    Purification_equipment_cost = Molecular_Sieve_HBC_cost + Cost_Comp1 + Cost_Comp2 + Cost_HX1 + Cost_HX2
    Purification_OperatingCost = (W_Comp1+W_Comp2)*3600*24*365/3600*Pelect
    Purification_OperatingCost_amortized = ((Purification_equipment_cost/ec_param['plant_lifetime'])+Purification_OperatingCost)
    
    return {'fc':fc_purif,
            'x':x_purif,
            'F':F_purif,
            'check_F':checks_F,
            'checks_x':checks_x,
            'tech':'Biogas purification',
            'fixed_capital_cost_2016':Purification_equipment_cost,
            'investment_cost':Purification_equipment_cost,
            'equipment_cost':Purification_equipment_cost,
            'operation_cost_2016_amortized':Purification_OperatingCost_amortized,
            'operation_cost_2016_non_amortized':Purification_OperatingCost,
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
            'MW_biogas_HX2':MW_biogas_HX2
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



