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

### FUNCTIONS USED IN MAIN ###

def stepint(mini,maxi,step):
    '''Outputs interval of values between a specific elements for a desired
    number of elements'''
    
    span = maxi - mini
    step_size = span / step
    step_int_out = list()
    step_int_out.append(mini)
    laststep = mini
    for i in range(0,step):
        step_int_out.append(laststep + step_size)
        laststep = laststep + step_size
    return step_int_out
        


def addict_initial(addictdeathrate,estaddictod,naturaldeath,presodrate2010,poprange):
    '''Calculates the proportion of the sample attributed to the addicted
    population at 2010'''
    
    #Based on Battista(2019) evaluation of addict pop at t0
    actual_pres_od_2010 = presodrate2010*poprange[0]
    initialaddict = estaddictod*actual_pres_od_2010 / poprange[0]
    initialaddict = initialaddict / (addictdeathrate - naturaldeath)
    
    return initialaddict


def addict_death_rate(estaddictod,actualod2015,naturaldeath,poprange):
    '''Calculates the death rate for the addicted population'''
    
    #Based on Battista(2019) evaluation of yearly addict death rate
    OP_SUD = 1.9/100 # northeast Opioid SUD proportion (Han 2015)
    addict_pop_2015 = poprange[5] * OP_SUD
    addict_death = estaddictod * actualod2015 * poprange[5]#2015
    addict_death = addict_death / addict_pop_2015
    addict_death = addict_death + (naturaldeath - estaddictod*actualod2015)
    rate = -1*math.log(1 - addict_death)
    
    return rate

### The following 6 functions are involved in our ODE solver ###

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
def func_s(t,x,s,addict_fromaddict_rate,a,addict_frompres_rate,p,j,rehab_success_rate,r,death_rate,addict_od_death_rate):
    '''returns RHS of ds/dt'''
    
    return -x*s -addict_fromaddict_rate*s*a -addict_frompres_rate*s*p + j*p + rehab_success_rate*r + death_rate*(p + r) + addict_od_death_rate*a
    
    
def func_p(t,x,s,j,addict_rate,death_rate,p):
    '''returns RHS of dp/dt'''
    
    return x*s -(j + addict_rate + death_rate)*p 
    

def func_a(t,addict_rate,p,relapse_rate,r,addict_fromaddict_rate,s,a,addict_frompres_rate,k,addict_od_death_rate):
    '''returns RHS of da/dt'''
    
    return addict_rate*p + relapse_rate*r + addict_fromaddict_rate*s*a + addict_frompres_rate*s*p -(k + addict_od_death_rate)*a
    


def func_r(t,k,a,rehab_success_rate,relapse_rate,death_rate,r):
    '''returns RHS of dr/dt'''
    
    return k*a -(rehab_success_rate + relapse_rate + death_rate)*r
    


def rk4_odesolver(poplist,listofparameters):
    '''solves the four ode's for the opioid epidemic for a step size of 1
    using the RK4 method using proportions at some t_n
    returns a tupple of our proportions at t_n+1'''
    
    #Local parameter assignement
    s,p,a,r = tuple(poplist) #my props at t_n
    x = listofparameters[0]
    j = listofparameters[1]
    k = listofparameters[2]
    addict_rate = listofparameters[3]
    addict_fromaddict_rate = listofparameters[4]
    addict_frompres_rate = listofparameters[5]
    relapse_rate = listofparameters[6]
    rehab_success_rate = listofparameters[7]
    death_rate = listofparameters[8]
    addict_od_death_rate = listofparameters[9]
    #Some values for solveivp to be used for all populations
    method = 'RK45'
    atol = 1.e-10
    rtol = 1.e-6
    t_min = 0.
    t_max = 1
    t_span= np.array([t_min, t_max])
    
    #Solve ODE's iterably:
    #1. Set initial condition for IV
    #2. Solve ODE using related variables, rates, and IC's(!!!)
    #3. Pull out solution for t = 1 (i.e. t = t_n + 1)
    for i in range(0,4):
        if i == 0:#S step
            s_0 = np.array([s])
            sol_s = solve_ivp(lambda t,s : func_s(t,x,s,addict_fromaddict_rate,a,addict_frompres_rate,p,j,rehab_success_rate,r,death_rate,addict_od_death_rate), t_span, s_0, method = method, rtol=rtol, atol=atol)  
            snew = sol_s.y[0,-1]
        elif i == 1:
            p_0 = np.array([p])
            sol_p = solve_ivp(lambda t,p : func_p(t,x,s,j,addict_rate,death_rate,p),t_span, p_0, method = method, rtol=rtol, atol=atol)
            pnew = sol_p.y[0,-1]
        elif i == 2:
            a_0 = np.array([a])
            sol_a = solve_ivp(lambda t,a : func_a(t,addict_rate,p,relapse_rate,r,addict_fromaddict_rate,s,a,addict_frompres_rate,k,addict_od_death_rate),t_span, a_0, method = method, rtol=rtol, atol=atol)
            anew = sol_a.y[0,-1]
        else:
            r_0 = np.array([r])
            sol_r = solve_ivp(lambda t,r : func_r(t,k,a,rehab_success_rate,relapse_rate,death_rate,r), t_span, r_0, method = method, rtol=rtol, atol=atol)
            rnew = sol_r.y[0,-1]
    
    #Ensure total = 1       
    new_pop = [snew,pnew,anew,rnew]
    sfin,pfin,afin,rfin = pop_assign(new_pop)
        
    return sfin,pfin,afin,rfin
        

def find_cutoff(resultsinstance,realmean,cutoffyear):#Used for 2
    '''Based on model output, choosing best neighborhood around the result
    to pick smallest option that fails to reject a null hypothesis of
    no difference @ p = 0.1 for over 48 degrees of freedom'''
    
    not_good = True
    myresults = resultsinstance
    cutoffrange = stepint(0.01,100,1000)
    t = 0
    cutoff = cutoffrange[t]
    if cutoffyear != 2017:#In case I go back to Option 3 with Means
        t_crit = 1.98
        desired_size = 100
    else:
        t_crit = 1.677#for significance level of 0.1
        desired_size = 50
    while not_good:
        deathinstance = []
        for i in range(0,len(myresults)):
            if (myresults[i][9] != 'NA') and (myresults[i][9] < cutoff) and myresults[i][0] == cutoffyear:
                deathinstance.append(myresults[i][8])
        
        #One Sample Student T test
        samplesize = len(deathinstance)
        if samplesize > 1:
            
            samplesum = 0
            for value in deathinstance:
                samplesum = samplesum + value
            samplemean = samplesum / samplesize
            
            diffsum = 0
            for value in deathinstance:
                diffsum = diffsum + (value - samplemean)**2   
            samplesd = (diffsum / (samplesize - 1))**0.5
            
            samplestderror = samplesd / (samplesize)**0.5
            
            calc_t = (samplemean - realmean) / samplestderror
            
            if abs(calc_t) < t_crit and samplesize >= desired_size:
                good_cutoff = cutoffrange[t]
                not_good = False
            elif t == 999:
                good_cutoff = 5
                not_good = False
                print('\nFAILED TO FIND CUTOFF FOR {}'.format(cutoffyear))
            else:
                t = t + 1
                cutoff = cutoffrange[t]
        else:
            t = t + 1
            cutoff = cutoffrange[t]
            
    return good_cutoff


def colmean(mylist,interestcol):#Used for option 3
    '''Calculates mean for column of interest in specific matrix'''
    
    col = []
    for row in mylist:
        col.append(row[interestcol])
    col_len = len(col)
    col_sum = 0
    for value in col:
        col_sum = col_sum + value
    col_mean = col_sum / col_len
    
    return col_mean


def find_next_value(mylist,interestcol):#Used for option 3
    '''Forms a linear trend-line based on a sample size of 6 (2011-2016)
    for column values of a parameter of interest
    NOTE: this method is NOT accurate, and will be changed'''
    
    dv = []
    for row in mylist:
        dv.append(row[interestcol])
    iv = [1,2,3,4,5,6]
    iv_sum = 0
    iv_sqr_sum = 0
    for value in iv:
        iv_sum = iv_sum + value
        iv_sqr_sum = iv_sqr_sum + value**2
    dv_sum = 0
    for value in dv:
        dv_sum = dv_sum + value
    iv_sqr_sum = 0
    iv_dv_sum = 0
    for i in range(0,len(iv)):
        iv_dv_sum = iv_dv_sum + iv[i]*dv[i]
    
    intercept = (dv_sum*iv_sqr_sum - iv_sum*iv_dv_sum) / (6*iv_sqr_sum - iv_sum**2)
    slope = (6*iv_dv_sum - iv_sum*dv_sum) / (6*iv_sqr_sum - iv_sum**2)
    next_val = intercept + slope*7
    
    return next_val

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
        paramlist = [x[0],x[1],x[5],x[4],x[3],x[2],x[7],x[6],x[8],x[9]]#For row of parameters
        for m in range(t_0 + 1,t_0 + 11):#Span of 10 years
          #For results in span from 2011

          pop_list = [s,p,a,r]
          #Find proportions for t0+1
          s_1,p_1,a_1,r_1 = rk4_odesolver(pop_list,paramlist)
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

        

        


