# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 11:13:47 2017

@author: rwest
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 11:04:16 2017

@author: rwest
"""


#import tkinter
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import csv
import tkinter as Tk
from scipy.stats import chisquare


#USEFUL DEFINITIONS
def indexfind(array,val):
    diff = np.abs(np.subtract(array,val))
    minimum = np.min(diff)
    s=0
    n=0
    while n <= len(diff) and s == 0:
        if minimum == diff[n]:
            ind = n
            s=1
            return ind
            break
        n=n+1

#LOADING DATA
data = np.loadtxt('acetonitrile.dat',delimiter='\t',unpack=True)
t = data[0,:]  #t = np.linspace(0,10,100)
S1ICTPPdat = data[1,:]
S1ICTPDPdat = data[2,:]
S1ICTPDPPPdat = data[3,:]

#MODEL CALCULATIONS
def PPeqns(N,t):
    #nS1 = N[0], nICT = N[1], nS0 = N[2]
    #K = S1-->ICT, ICT-->S1, S1-->S0, ICT-->S0
    dS1dt = K[1]*N[1] - K[0]*N[0] - K[2]*N[0]
    dICTdt = K[0]*N[0] - K[1]*N[1] - K[3]*N[1]
    dS0dt = K[2]*N[0] + K[3]*N[1]
    return [dS1dt,dICTdt,dS0dt]

def PDPeqns(ND,t):
     #nS1 = N[0], nICT = N[1], nS0 = N[2]
     #K = S1-->ICT, ICT-->S1, S1-->S0, ICT-->S0   
    dS1Ddt = K[1]*ND[1] - K[0]*ND[0] - K[2]*ND[0]
    GaussD = np.exp(np.dot(-2.77,(np.square(np.divide((t-tD),DPW)))))
    dICTDdt = K[0]*ND[0] - K[1]*ND[1] - K[3]*ND[1] - Dpow*GaussD*kD*ND[1]
    dS0Ddt = K[2]*ND[0] + K[3]*ND[1] + Dpow*GaussD*kD*ND[1]
    return [dS1Ddt,dICTDdt,dS0Ddt]

#definition for removing zeroes for divided arrays  #new
def remove0(array):
    for n in range(0,len(array)):
        if array[n] == 0:
            array[n] == 0.0000000001

#USER INPUT
#USER OPTIONS
#apply common and equilibrium parameters
applyLSQrat = 0         #U include LSQ fit of S1ICTPP and S1ICTPDP in total LSQ
applyLSQratrat = 1      #U include LSQ fit of S1ICTPDPPP in total LSQ
#fix kinetic parameters
fixS1ICT = 0
fixICTS1 = 0
fixS1S0 = 1
fixICTS0 = 0
fixnS1ICT0 = 1
fixcoeq = 0             #U apply two-level model strictly
fixDpow = 0 #new
#USER ENTRY
#dump pulse properties
Dpow = 500               #U
tD = 1.91               #U
DPW = 0.065             #U
kD = 0.12   #new
#common and equilibration decay
Tco = 50.461            #U
Teq = 7.927             #U
#-----------------------------
TS1ICT = 25         #U
TICTS1 = 20          #U
TS1S0 = 68              #U
TICTS0 = 15             #U

#initial populations
S10 = 100
RnS1ICT0 = 1.0          #U
S00 = 0

#strictly apply two-level system model
C = 0.25*((1/Tco)**2 + (1/Teq)**2 - 2*(1/Tco)*(1/Teq) - 1/Tco - 1/Teq)
if fixcoeq == 1:
    fixICTS0 = 1
    fixnS1ICT0 = 1
    #at least begin with two-level system model
    TS1ICT = -(1/TICTS0-1/TS1S0) / (C + (1/TS1S0)**2 - (1/Tco+1/Teq)*(1/TS1S0))
    TICTS1 = (1/TS1S0-1/TICTS0) / (C + (1/TICTS0)**2 - (1/Tco+1/Teq)*(1/TICTS0))
elif fixcoeq == 0 and fixS1ICT == 1:
    TICTS1 = (1/TS1S0-1/TICTS0) / (C + (1/TICTS0)**2 - (1/Tco+1/Teq)*(1/TICTS0))
elif fixcoeq == 0 and fixICTS1 == 1:
    TS1ICT = -(1/TICTS0-1/TS1S0) / (C + (1/TS1S0)**2 - (1/Tco+1/Teq)*(1/TS1S0))
    

#     T[0]   T[1]   T[2]  T[3]
T = [TS1ICT,TICTS1,TS1S0,TICTS0]
    
#FITTING ALGORITHM
#percent of 'lifetimes' to test above and below
Tper = 60
Tres = 60
#points which to average to the right for Dpow fitting  #new
avgDpow = 4    #new
#prepare ranges of values for fitting loop
T0range = np.linspace(T[0]-np.round(T[0]*Tper,2),T[0]+np.round(T[0]*Tper,2),2*Tres)
T1range = np.linspace(T[1]-np.round(T[1]*Tper,2),T[1]+np.round(T[1]*Tper,2),2*Tres)
T2range = np.linspace(T[2]-np.round(T[2]*Tper,2),T[2]+np.round(T[2]*Tper,2),2*Tres)
T3range = np.linspace(1,80,1)
RnS1ICT0range = np.linspace(0,5,10)
#conditions for fixed rate parameters
if fixS1ICT == 1:
    T0range = [TS1ICT]
if fixICTS1 == 1:
    T1range = [TICTS1]
if fixS1S0 == 1:
    T2range = [TS1S0]
if fixICTS0 == 1:
    T3range = [TICTS0]
if fixnS1ICT0 == 1:
    RnS1ICT0range = [RnS1ICT0]
#total number of computations  #new
totcount = len(T0range)*len(T1range)*len(T2range)*len(T3range)*len(RnS1ICT0range)
#for normalization of Dpow to max of S1ICTPDPPPdat
spikeind = indexfind(S1ICTPDPPPdat,np.max(S1ICTPDPPPdat)) + 1
#setting high value for LSQ to begin finding lowest LSQ value
LSQprev = 99999
#fitting loop
#current computation # count
count = 0    #new
for RnS1ICT0 in RnS1ICT0range:
    for T0 in T0range:
        T[0] = T0
        for T1 in T1range:
            T[1] = T1
            for T2 in T2range:
                T[2] = T2
                for T3 in T3range:
                    T[3] = T3
                    #increase computation number count  #new
                    count=count+1
                    if totcount % count == 0:
                        print(count,' / ',totcount)
                    #varied initial ICT population
                    ICT0 = S10*RnS1ICT0
                    #     S1  ICT  S0
                    N0 = [S10, ICT0, S00]
                    if fixcoeq == 1:
                        C = 0.25*((1/Tco)**2 + (1/Teq)**2 - 2*(1/Tco)*(1/Teq) - 1/Tco - 1/Teq)
                        T[0] = -(1/TICTS0-1/TS1S0) / (C + (1/TS1S0)**2 - (1/Tco+1/Teq)*(1/TS1S0))
                        T[1] = (1/TS1S0-1/TICTS0) / (C + (1/TICTS0)**2 - (1/Tco+1/Teq)*(1/TICTS0))
                    #rate array
                    K = [1/T[0], 1/T[1], 1/T[2], 1/T[3]]
                    #Solve Differential equations
                    PP = odeint(PPeqns,N0,t)
                    PDP = odeint(PDPeqns,N0,t)
                    
                    #MANIPULATING MODEL RESULTS
                    #calculate ratios
                    S1ICTPP = np.divide(PP[:,0],PP[:,1])
                    S1ICTPDP = np.divide(PDP[:,0],PDP[:,1])
                    S1ICTPDPPP = np.divide(S1ICTPDP,S1ICTPP)
                    #normalizing of Dpow
                    if fixDpow == 0:   #new
                        D = np.mean(S1ICTPDPPPdat[spikeind:spikeind+avgDpow])
                        F = np.mean(S1ICTPDPPP[spikeind:spikeind+avgDpow])
                        if D-F >= 0:
                            Dpow = ((D-F)/F+1)*Dpow
                        else:
                            Dpow = (1-(F-D)/F)*Dpow
                    #Least square values of fit (check out module on internet) {always display this}
                    if applyLSQrat == 1 or (applyLSQrat == 0 and applyLSQratrat == 0):
                        LSQS1ICTPP = np.sqrt(np.sum(np.square(np.subtract(S1ICTPP,S1ICTPPdat))))
                        LSQS1ICTPDP = np.sqrt(np.sum(np.square(np.subtract(S1ICTPDP,S1ICTPDPdat))))
                        prefS1ICT = 'S1ICT'
                    else:
                        LSQS1ICTPP = 0
                        LSQS1ICTPDP = 0
                        prefS1ICT = ''
                    if applyLSQratrat == 1 or (applyLSQrat == 0 and applyLSQratrat == 0):
                        LSQS1ICTPDPPP = np.sqrt(np.sum(np.square(np.subtract(S1ICTPDPPP,S1ICTPDPPPdat))))
                        prefS1ICTPDPPP = 'S1ICTPDPPP'
                    else:
                        LSQS1ICTPDPPP = 0
                        prefS1ICTPDPPP = ''
                    LSQtotal = LSQS1ICTPP + LSQS1ICTPDP + LSQS1ICTPDPPP
                    LSQ = [LSQS1ICTPP,LSQS1ICTPDP,LSQS1ICTPDPPP,LSQtotal]
                    if LSQtotal < LSQprev:
                        LSQfit = LSQ
                        LSQprev = LSQtotal
                        PPfit = PP
                        PDPfit = PDP
                        S1ICTPPfit = S1ICTPP
                        S1ICTPDPfit = S1ICTPDP
                        S1ICTPDPPPfit = S1ICTPDPPP
                        RnS1ICT0fit = RnS1ICT0
                        Tfit = T
                        Dpowfit = Dpow

#Saving best fit values 
PP = PPfit
PDP = PDPfit
S1ICTPP = S1ICTPPfit
S1ICTPDP = S1ICTPDPfit
S1ICTPDPPP = S1ICTPDPPPfit

#OUTPUT
#QUALITY OF FIT NUMERICS
X2S1ICTPP = chisquare(S1ICTPP,S1ICTPPdat)
X2S1ICTPDP = chisquare(S1ICTPDP,S1ICTPDPdat)
X2S1ICTPDPPP = chisquare(S1ICTPDPPP,S1ICTPDPPPdat)





print('Chi-squared values \n  S1ICTPP: ',X2S1ICTPP,'\n  S1ICTPDP: ',X2S1ICTPDP,'\n  S1ICTPPPDP: ',X2S1ICTPDPPP)
#PRINTING NUMBERS
print('LSQ: [S1ICTPP,S1ICTPDP,S1ICTPDPPP,total]')
print('LSQ:',np.round(LSQfit,4))
print('initial nICT % of S1',RnS1ICT0fit)
print('lifetimes:',Tfit)
print('Dpow = ',Dpowfit)
print('Fit to graph(s): ',prefS1ICT,' ',prefS1ICTPDPPP)







#PLOTTING
#plot solutions
popfig = plt.figure(figsize=(6,6))
plt.plot(t,PP[:,0],color='b',label='S1')
plt.plot(t,PP[:,1],color='r',label='ICT')
plt.plot(t,PP[:,2],color='k',label='S0')
plt.plot(t,PDP[:,0],':',color='b')
plt.plot(t,PDP[:,1],':',color='r')
plt.plot(t,PDP[:,2],':',color='k')
plt.xlabel('time (ps)')
plt.ylabel('population')
plt.title('Populations')
plt.legend(loc=4)

ratfig = plt.figure(figsize=(6,6))
S1ICTfig = ratfig.add_subplot(2,1,1)
plt.plot(t,S1ICTPP,color='r',label='PP')
plt.plot(t,S1ICTPPdat,color='k')
plt.plot(t,S1ICTPDP,':',color='r',label='PDP')
plt.plot(t,S1ICTPDPdat,':',color='k')
plt.xlabel('time (ps)')
plt.ylabel('nS1 / nICT')
plt.title('S1/ICT ratio of pump-probe and pump-dump-probe')
plt.legend(loc=4)

PDPPPfig = ratfig.add_subplot(2,1,2)
plt.plot(t,S1ICTPDPPP,color='r')
plt.plot(t,S1ICTPDPPPdat,color='k')
plt.xlabel('time (ps)')
plt.ylabel('PDP(nS1/nICT) / PP(nS1/nICT)')
plt.title('Ratio of S1/ICT ratios (PDP/PP)')