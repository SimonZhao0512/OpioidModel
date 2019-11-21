#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 11:43:07 2019

@author: noahigra
"""
from __future__ import division
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import ttest_1samp as ttest
from scipy.stats import f_oneway

### FUNCTIONS USED IN MAIN ###

def progress_report(current_step,totalstep,howmanyreports):
    '''Returns a progress report a specified amount of times'''
    
    steps_per_report = totalstep//howmanyreports
    if current_step % steps_per_report == 0:
        report_percent = math.ceil((current_step / totalstep)*100)
        print('Progress: {}%'.format(report_percent))   



def addict_initial(addictdeathrate,estaddictod,naturaldeath,presodrate2010,poprange):
    '''Calculates the proportion of the sample attributed to the addicted
    population at 2010'''
    
    #Based on Battista(2019) evaluation of addict pop at t0
    actual_pres_od_2010 = presodrate2010*poprange[2010]
    initialaddict = estaddictod*actual_pres_od_2010 / poprange[2010]
    initialaddict = initialaddict / (addictdeathrate - naturaldeath)
    
    return initialaddict


def addict_death_rate(estaddictod,actualod2015,naturaldeath,poprange):
    '''Calculates the death rate for the addicted population'''
    
    #Based on Battista(2019) evaluation of yearly addict death rate
    OP_SUD = 1.9/100 # northeast Opioid SUD proportion (Han 2015)
    addict_pop_2015 = poprange[2015] * OP_SUD
    addict_death = estaddictod * actualod2015 * poprange[2015]#2015
    addict_death = addict_death / addict_pop_2015
    addict_death = addict_death + (naturaldeath - estaddictod*actualod2015)
    rate = -1*math.log(1 - addict_death)
    
    return rate


### The following functions are involved in our ODE solver ###
    
def pop_assign(listofpop):
    '''Takes a list of population proportions and assigns them to individual 
    variables, returns a tupple'''
    
    news = listofpop[0]
    newp = listofpop[1]
    newa = listofpop[2]
    newr = listofpop[3]
    
    #Make sure that new proportions satisfy s + p + a + r = 1
    newpop = news + newp + newa + newr
    mys,myp,mya,myr = news / newpop, newp / newpop, newa / newpop, newr / newpop
    
    return mys,myp,mya,myr



#The following func_s-r are used in the RK4 solver

def g_s(t,pres_rate,s,p,a,r,k,g,
           no_addict_rate,relapse_rate,death_rate,addict_od_death_rate):
    '''Returns RHS of ds/dt'''
    
    return -(pres_rate+a+p)*s + (no_addict_rate+death_rate)*p/k + pres_rate*((1/relapse_rate)-1)*r + addict_od_death_rate*a/g


def g_p(t,pres_rate,s,p,k,no_addict_rate,addict_rate,death_rate):
    '''Returns RHS of dp/dt'''
    
    return pres_rate*k*s - (no_addict_rate + addict_rate + death_rate)*p


def g_a(t,addict_rate,s,p,a,r,pres_rate,relapse_rate,k,
           g,h,treat_entry_rate,addict_od_death_rate):
    '''Returns RHS of da/dt'''
    
    return pres_rate*g*r + (h*addict_rate + g*s)*p + (g*s - (treat_entry_rate + addict_od_death_rate))*a


def g_r(t,treat_entry_rate,a,r,pres_rate,g,relapse_rate,death_rate):
    '''Returns RHS of dr/dt'''
    
    return treat_entry_rate*relapse_rate*a/(pres_rate*g) - r


def dg_dt(t,poplist,alpha,epsilon,iota,gamma,k,g,h,sigma,
              addict_od_death_rate,death_rate):
    '''Returns RHS of dg/dt for 4x1 vector of ODE'S'''
    s,p,a,r = poplist[0], poplist[1], poplist[2],poplist[3]
    
    gs = g_s(t,alpha,s,p,a,r,k,g,epsilon,sigma,death_rate,
             addict_od_death_rate)
    gp = g_p(t,alpha,s,p,k,epsilon,gamma,death_rate)
    ga = g_a(t,gamma,s,p,a,r,alpha,sigma,k,
           g,h,iota,addict_od_death_rate)
    gr = g_r(t,iota,a,r,alpha,g,sigma,death_rate)
    
    return [gs,gp,ga,gr]


    
    
    
def rk4_odesolver(poplist,listofparameters):
    '''solves the 4-D ODE for the opioid epidemic for a step size of 1
    using the RK4 method using proportions at some t_n
    returns a tupple of our proportions at t_n+1
    approach options:
    'iteration' -- solve's the 4 ODE's 1 by 1, keeping 3/4 pop parameters
    constant
    'allatonce' -- solves a 4x1 ODE '''
    
    #Local parameter assignement
    alpha = listofparameters[0]
    epsilon = listofparameters[1]
    iota = listofparameters[2]
    gamma = listofparameters[3]
    beta_a = listofparameters[4]
    beta_p= listofparameters[5]
    sigma = listofparameters[6]
    death_rate = listofparameters[7]
    addict_od_death_rate = listofparameters[8]
    
    k = beta_p / alpha
    h = beta_a / beta_p
    g = h*k
    
    #Some values for solveivp to be used for all populations
    method = 'RK45'
    atol = 1.e-10
    rtol = 1.e-6
    t_min = 0.
    t_max = 1
    t_span= np.array([t_min, t_max])
    
    #scaling
    poplist[0] = poplist[0]*alpha
    poplist[1] = poplist[1]*alpha*k
    poplist[2] = poplist[2]*alpha*g
    poplist[3] = poplist[3]*sigma
    
    
    
    
    pop_n_1 = solve_ivp(lambda t,poplist : 
        dg_dt(t,poplist,alpha,epsilon,iota,gamma,k,g,h,sigma,
              addict_od_death_rate,death_rate),
        t_span,np.array(poplist),
        method = method, rtol=rtol, atol=atol,vectorized = True)
        
    values_n_1 = pop_n_1.y[:,-1]
    sfin,pfin,afin,rfin = values_n_1[0], values_n_1[1], values_n_1[2], values_n_1[3]
    
    sfin = sfin/alpha
    pfin = pfin/(alpha*k)
    afin = afin/(alpha*g)
    rfin = rfin/sigma
            
    sfin,pfin,afin,rfin = pop_assign([sfin,pfin,afin,rfin])
            
    return sfin,pfin,afin,rfin
            


    
      

def find_cutoff(resultsinstance,realmean,cutoffyear,desired_size = 30):#Used for 2
    '''Based on model output, choosing best neighborhood around the result
    to pick smallest option that fails to reject a null hypothesis of
    no difference @ p = 0.1'''
    
    not_good = True
    myresults = resultsinstance
    cutoffrange = np.linspace(0.03,300,1000,endpoint = True).tolist()
    t = 0
    cutoff = cutoffrange[t]
    

    while not_good:

        result_subset = myresults[myresults[:,0] == cutoffyear]
        result_subset = result_subset[result_subset[:,9] < cutoff]
        deathinstance = result_subset[:,8]
        
        #One Sample Student T test
        samplesize = deathinstance.shape[0]
        
        if samplesize > 1:
            
            calc_result = ttest(deathinstance,realmean)
            if calc_result.pvalue > 0.1 and samplesize >= desired_size:
                good_cutoff = cutoffrange[t]
                not_good = False
            elif t == 999:
                good_cutoff = 25
                not_good = False
                print('\nFAILED TO FIND CUTOFF FOR {}'.format(cutoffyear))
            else:
                t = t + 1
                cutoff = cutoffrange[t]
        else:
            t = t + 1
            cutoff = cutoffrange[t]
            
            
    return good_cutoff


def anova(allresults,categories):
    '''Evaluates ANOVA of no difference between subgroups (by category)
    returns boolean failing to reject null @ 0.05'''
    
    group_1 = allresults[allresults[:,1] == categories[0]][:,8]
    group_2 = allresults[allresults[:,1] == categories[1]][:,8]
    group_3 = allresults[allresults[:,1] == categories[2]][:,8]
    group_4 = allresults[allresults[:,1] == categories[3]][:,8]
    group_5 = allresults[allresults[:,1] == categories[4]][:,8]
    
    instant_anova = f_oneway(group_1,group_2,group_3,group_4,group_5)
    p_value = instant_anova.pvalue
    
    
    if p_value > 0.05:
        return True
    else:
        return False
    
    






#This function is the RK4 solver method applied in main() using numpy arrays
#as the input and output

def sense_RK4(param_values,currentpop,initialconditions):
    '''For given population, find population values for saltelli matrix of 
    N parameters for 10 years, output as list of values'''
    
    t_0,s_0,p_0,a_0,r_0 = initialconditions #Initial conditons
    
    #Setup N x D matrix to receive results
    s_1 = np.zeros([param_values.shape[0]])
    p_1 = np.zeros([param_values.shape[0]])
    a_1 = np.zeros([param_values.shape[0]])
    r_1 = np.zeros([param_values.shape[0]])

    #Reset to t_0
    s = s_0
    p = p_0
    a = a_0
    r = r_0

    og_pop_set = (s,p,a,r)
    s,p,a,r = og_pop_set 
    
    A = [] #Initiate output matrix
    
    for i, x in enumerate(param_values):
        paramlist = [x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9]]#For row of parameters
        for m in range(t_0 + 1,t_0 + 11):#Span of 10 years
          #For results in span from 2011

          pop_list = [s,p,a,r]
          #Find proportions for t0+1
          s_1,p_1,a_1,r_1 = rk4_odesolver(pop_list,paramlist,'irregular')
          s,p,a,r = s_1,p_1,a_1,r_1
        
        #Using our index, choose which population type to output for  
        if currentpop == 1:
             A.append(s)
        elif currentpop == 2:
            A.append(p)
        elif currentpop == 3:
            A.append(a)
        else:
            A.append(r)

    return A #A is an Nx1 vector of a given population, to be used in sobol analysis

        

        


