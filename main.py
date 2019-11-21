#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 16:35:13 2019

@author: noahigra
"""
import opioidmodule as op
from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np
import csv



def main():
    '''Program will take in a specific table of values and output data to
    be visualized in R'''
    
    print('\nWelcome to The Prescription Opioid Epidemic Model')
    print('Brought to you by Simon and the boiz')
    
    whattodo , whattoprint = menu()#Select Desired Output
    
    run = whattodo != 4
    
    if run: #If not quiting
        
        ###DATA DUMP###
        input_csv = input('Please select your dataset: ')
        param_dict = parameter_plug(input_csv)
        #Depending on selection, will define our variable parameters accordingly
        
        #Populations from 2010 - 2017
        pop_2010_2017 = {2010:param_dict['pop_2010'],2011:param_dict['pop_2011'],
                         2012:param_dict['pop_2012'],2013:param_dict['pop_2013'],
                         2014:param_dict['pop_2014'],2015:param_dict['pop_2015'],
                         2016:param_dict['pop_2016'],2017:param_dict['pop_2017']}
        #Estimated proportion of prescription opioid overdoses attributed to addicts
        EST_ADDICT_OD = param_dict['est_addict_od']
        #Actual overdoses from 2011-2017 and compiling those values into a list
        opioid_od_2011 = param_dict['opioid_od_2011']*EST_ADDICT_OD
        opioid_od_2012 = param_dict['opioid_od_2012']*EST_ADDICT_OD
        opioid_od_2013 = param_dict['opioid_od_2013']*EST_ADDICT_OD
        opioid_od_2014 = param_dict['opioid_od_2014']*EST_ADDICT_OD
        opioid_od_2015 = param_dict['opioid_od_2015']*EST_ADDICT_OD
        opioid_od_2016 = param_dict['opioid_od_2016']*EST_ADDICT_OD
        opioid_od_2017 = param_dict['opioid_od_2017']*EST_ADDICT_OD
        addict_od_list = [opioid_od_2011,opioid_od_2012,opioid_od_2013,
                          opioid_od_2014,opioid_od_2015,opioid_od_2016,
                          opioid_od_2017]#for modeling
        od_dict = {2011:opioid_od_2011, 2012:opioid_od_2012, 2013:opioid_od_2013,
                   2014:opioid_od_2014, 2015:opioid_od_2015, 2016:opioid_od_2016,
                   2017:opioid_od_2017}

        if whattodo !=3:#For all that isn't numerical sensitivity
                        
            addict_rate = param_dict['addict_rate']#Rate of addiction from prescription
            addict_fromaddict_rate = param_dict['addict_fromaddict_rate']#Rate of addiction from illicit drugs from addicts 
            addict_frompres_rate = param_dict['addict_frompres_rate']#Rate of addiction from illicit drugs from prescribed users 
            relapse_rate = param_dict['relapse_rate']#Rate of relapse

            if whattodo == 1:
                pres_rate_list = [param_dict['pres_rate']]#Prescription rate
                no_addict_rate_list = [0.005,0.5,1.0,1.5,2.0,2.5,3]#Rate of ending prescription without addiction
                treat_entry_rate_list = np.linspace(0.001,1,
                                                    300,endpoint = True).tolist()#Treatment entry rate
            else:#Correlation plot
                pres_rate_list = np.linspace(param_dict['pres_rate']-0.1,param_dict['pres_rate']+0.05,
                                        5,endpoint = True).tolist() 
                no_addict_rate_list = np.linspace(0.008,8,
                                             300,endpoint = True).tolist()
                treat_entry_rate_list = np.linspace(0.005,5,
                                               300,endpoint = True).tolist()
            


        #Values needed to determine initial conditions    
        death_rate = param_dict['natural_death_rate']#Natural death rate
        actual_od_death_2015 = param_dict['actual_od_death_2015']#Rate of actual overdose in 2015
        presopioid_death_rate_2010 = param_dict['presopioid_death_rate_2010']#Rate of actual overdose in 2010
        #Death rate attributed to the addicted class
        addict_od_death_rate = op.addict_death_rate(EST_ADDICT_OD,actual_od_death_2015,
                                                     death_rate,pop_2010_2017)
        
        
        if whattodo == 3:#Parameter choice for numerical sensitivity
            parameter_names = ['pres_rate', 'no_addict_rate','treat_entry_rate',
                               'addict_rate','illictaddict/pres','illictpres/pres',
                               'relapse_rate','illictaddict/illicitpres',
                               'death_rate','addict_od_death_rate']
            
            parameter_ranges = {'num_vars':10,'names': parameter_names,
                                'bounds':[[(param_dict['pres_rate']+0.1)/10,
                                           param_dict['pres_rate']+0.1],[0.3,3],
                                           [0.2*param_dict['relative_treat'],2*param_dict['relative_treat']],
                                           [0.00235,0.0235],
                                           [(param_dict['addict_fromaddict_rate']/param_dict['pres_rate'])*3/10,
                                            (param_dict['addict_fromaddict_rate']*3/param_dict['pres_rate'])],
                                           [(param_dict['addict_frompres_rate']/param_dict['pres_rate'])*3/10,
                                            (param_dict['addict_frompres_rate']*3/param_dict['pres_rate'])],
                                           [0.001,1],[0.01,3],
                                           [(death_rate*3)/10,death_rate*3],
                                           [(addict_od_death_rate*3)/10,addict_od_death_rate*3]]}
                        
  
        
        t_0 = 2010
        span = 8#end at 2017
        p_0 = param_dict['p_0']#Proportion of population prescribed opioids per week
        a_0 = op.addict_initial(addict_od_death_rate,EST_ADDICT_OD,death_rate,
                             presopioid_death_rate_2010,pop_2010_2017)
        r_0 = 0.1*a_0
        s_0 = 1 - (p_0 + a_0 + r_0)
        
        
        ###END OF DATA DUMP###
        
        #Some opening notes, show user initial conditions
        print('\nModeling epidemic starting at 2010:')
        print('s = {:.4f}'.format(s_0))
        print('p = {:.4f}'.format(p_0))
        print('a = {:.4f}'.format(a_0))
        print('r = {:.4f}'.format(r_0))
        print('\n### MODEL IN PROCESS, PLEASE HOLD ###\n')
        
              
        ### ACTUAL MODELING ###
        
        #Initiate list for results
        results = []
        #Iterate through range of variables
        if whattodo == 1 or whattodo == 2 :#For options 1 and 2, same procedure
            step = 0 #initiate step_count
            step_total = len(pres_rate_list)*len(treat_entry_rate_list)*len(no_addict_rate_list)*7 # Total iterations, not accounting for internal manipulations
            for pres_rate in pres_rate_list: 
                for treat_entry_rate in treat_entry_rate_list:
                    for no_addict_rate in no_addict_rate_list:
                        #Reset proportions to t = t_0
                        s,p,a,r = s_0, p_0, a_0, r_0
                        
                        #Compile all of our parameters into a list
                        paramlist = [pres_rate,no_addict_rate,treat_entry_rate,
                                     addict_rate,addict_fromaddict_rate,
                                     addict_frompres_rate,relapse_rate,
                                     death_rate,addict_od_death_rate]
                        
                        #Calculate RK4 Steps for fixed set of parameters across span of time
                        for year in range(t_0+1,(t_0 + span)):#For results in span from 2011
                            #Setup lists to be used in rk4 solver
                            pop_list = [s,p,a,r]
                            #Find proportions for t0+1
                            s_1,p_1,a_1,r_1 = op.rk4_odesolver(pop_list,paramlist)
                            
                            #Calculate S-ORD(t)
                            od_death = (pop_2010_2017[year])*(addict_od_death_rate - death_rate)*a_1
                            od_dif = abs(od_death - od_dict[year])#Compare to A-ORD(2017)
                            results.append([year,
                                            pres_rate,no_addict_rate,treat_entry_rate,
                                            s_1,p_1,a_1,r_1,od_death,od_dif])
                  
        
                            #Return to vogue, assign initial conditions for next year
                            s,p,a,r = s_1, p_1, a_1, r_1

                            step = step + 1
                            op.progress_report(step,step_total,10)
        
        else:#Output for Numerical Sensitivity
            
            og_stat_tuple = t_0,s_0,p_0,a_0,r_0#Initial conditions tuple
            #Fixed 1000X10 (Number of Runs x Number of Variables) Matrix of parameters 
            param_matrix = saltelli.sample(parameter_ranges, 1000)
            
            sense_results = []#Initiate list of results
            step = 0
            pop_index = (1,2,3,4)#Index to be used to choose sensitivity output
            for pop in pop_index:
                
                #Solve RK4 using saltelli amtrix of parameters for 1 population of 4
                pop_results = np.array(op.sense_RK4(param_matrix,pop,og_stat_tuple))
                #Find S1 and ST for 1 population out of 4
                pop_sense = sobol.analyze(parameter_ranges,pop_results,print_to_console = False)
                s1_list = list(pop_sense['S1'])
                st_list = list(pop_sense['ST'])
                
                if pop == 1:
                    value = ['s']
                elif pop == 2:
                    value = ['p']
                elif pop == 3:
                    value = ['a']
                else:
                    value = ['r']
                
                first_order_row = value + ['1'] + s1_list 
                total_order_row = value + ['2'] + st_list

                sense_results.append(first_order_row)#Row of 1 sensitivities for give pop
                sense_results.append(total_order_row)#Row of T sensitivities for given pop
                
                if whattoprint != 1 and whattoprint != 3:#for statistical comparisons
                    s1_conf_list = list(pop_sense['S1_conf'])
                    st_conf_list = list(pop_sense['ST_conf'])
                    first_order_conf_row = value + ['3'] + s1_conf_list
                    total_order_conf_row = value + ['4'] + st_conf_list
                    sense_results.append(first_order_conf_row)
                    sense_results.append(total_order_conf_row)
                
                step = step + 1
                op.progress_report(step,4,4)
        
        
        print('\n### MODEL COMPLETE, MANIPULATING RESULTS NOW ### ')
        
        ### MANIPULATIONS ###

        
        if whattodo == 2:#correlation plot
            results = np.array(results)
            #initialize list of close results
            close_od = np.zeros([1,10])
            for pres_rate in pres_rate_list:
                results_by_rate = results[results[:,1] == pres_rate]
                #Find statistically significant cutoff                
                cutoff = op.find_cutoff(results_by_rate,opioid_od_2017,2017,70)
                print('\nCutting off results for prescription rate {:.2f}'.format(pres_rate))
                print('@ Difference < {:.2f}'.format(cutoff))
                results_cutoff = results_by_rate[results_by_rate[:,0] == 2017]
                results_cutoff = results_cutoff[results_cutoff[:,9] <= cutoff]

                close_od = np.vstack((close_od,results_cutoff))
        
            
        ### OUTPUT ###
        
        #Output results as csv files to be plotted and analyzed in R
        #also tell user how many rows of results have been calculated
        if whattodo == 1:
            #Output all results
            if whattoprint == 1 or whattoprint == 3:
                csv_output(results,'modeloutput',1,addict_od_list)
                print('\ncreated modeloutput.csv with {} rows'.format(len(results)))
            if whattoprint == 2 or whattoprint == 3:
                csv_output(results,'modeldataset',1,addict_od_list)
                print('\ncreated modeldataset.csv with {} rows'.format(len(results)))
        elif whattodo == 2:
            close_od = np.delete(close_od,0,0)
            #Perform ANOVA on list within neighborhood
            anova_result = op.anova(close_od,pres_rate_list)
            close_od = list(close_od)
            if not anova_result:
                print('FAILED TO FIND CLOSE RESULTS')
            if whattoprint == 1 or whattoprint == 3:
                csv_output(close_od,'modeloutput',2,addict_od_list)
                print('\ncreated modeloutput.csv with {} rows'.format(len(close_od)))
            if whattoprint == 2 or whattoprint == 3:
                 csv_output(close_od,'modeldataset',2,addict_od_list)
                 print('\ncreated modeldataset.csv with {} rows'.format(len(close_od)))
        else:#Option 3
            if whattoprint == 1 or whattoprint == 3:
                csv_output(sense_results,'modeloutput',3,addict_od_list)
                print('\ncreated modeloutput.csv with {} rows'.format(len(sense_results)))
            if whattoprint == 2 or whattoprint == 3:
                csv_output(sense_results,'modeldataset',3,addict_od_list)
                print('\ncreated modeldataset.csv with {} rows'.format(len(sense_results)))
        
        print('\n### COMPLETE ###\n')
        
            
    else:#QUIT
        
        print('\nThanks for Playing!')
        print('\n### END OF PROGRAM ###\n')
    

###END OF MAIN###
    

    
#Following functions are used but not directly involved in any computations

def parameter_plug(mycsv):
    '''Converts csv table of coefficients into a dictionary'''
    
    all_param = {}
    
    with open(mycsv,'r') as csv_file:
        data = csv.reader(csv_file, delimiter = ',')
        
        for line in data:
            #Make dictionary of values s.t the keys match the table key
            #inputformat is the proper format for such a table, but any order
            #of these specifically named parameters would work
            param = line[0]
            all_param[param] = float(line[1])
        
    return all_param

    
def csv_output(mylist,title,outputtype,carryoverlist):#Used for all
    '''Outputs results as a csv'''
    
    csv_name = str(title) + '.csv'

    with open(csv_name,'w') as results_csv:
        writer = csv.writer(results_csv)
        lencarry = len(carryoverlist)
        firstrow = (len(mylist[0]) - lencarry)*[outputtype] + carryoverlist
        writer.writerow(firstrow)

        for row in mylist:
            writer.writerow(row)
            




def menu():
    '''User chooses desired activity'''
    
    choice_dict = {1: 'Overdose per Year', 2: 'Parameter Correlation', 3:'Sensitivity Analysis'}
    output_dict = {1: 'Plot', 2: 'Dataset', 3:'Plot and Dataset'}
    print('\n### MENU ###\n')
    print('Please Choose a Action:')
    print('1) Calculate Overdose per Year 2011-2017')
    print('2) Find Parameter Correlations for Accurate Models ')
    print('3) Run First and Total Order Sensitivity Analysis')
    print('4) Quit Program')
    choice_bad = True
    good_pick = 'n'
    while choice_bad or good_pick == 'n':
        choice = int(input('Choice: '))
        if choice in [1,2,3,4]:
            choice_bad = False
            if choice != 4:
                print('\nPlease Choose a Desired Output:')
                print('1) Plot')
                print('2) CSV Dataset')
                print('3) Both')
                output_choice = int(input('Choice: '))
                if output_choice in [1,2,3]:
                    print('\nSelection: {}\nOutput: {}'.format(choice_dict[choice],output_dict[output_choice]))
                    good_pick = input('Would you like to continue?(y/n) ').strip()
                else:
                    print('{} is not an appropriate option'.format(output_choice))
                    print('Please pick an item from the menu above!')
            else:
                print('\nSelection: Quit')
                good_pick = input('Would you like to continue?(y/n) ').strip()
                output_choice = 'NA'
                
        else:
            print('{} is not an appropriate option'.format(choice))
            print('Please pick an item from the menu above!')

    return choice , output_choice   

    
main()
