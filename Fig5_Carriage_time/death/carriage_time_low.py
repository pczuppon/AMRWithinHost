import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.integrate as integrate
from numpy import genfromtxt
from scipy.special import hyp2f1
from scipy.special import expi

################################
################################ data
################################

os.chdir(os.path.realpath(''))

################################
################################ Parameter definitions
################################

################################ Logistic growth dynamics
bS = 2.5
bR = 2.25
dS = 0.5
dR = 0.5

gamma = 1.

K = 1000.

N0 = (bS-dS)*K/gamma

TT = 7

################################ Antimicrobial response / Pharmacodynamics
c = np.arange(0,1,0.0001)

def a(c,X):                                 # X = 0 for sensitive growth rate, X = 1 for resistant growth rate
    if (X==0):
        psimax = bS-dS
        mic = 0.017
    if (X==1):
        psimax = bR-dR
        mic = 0.017*5
    psimin = np.log(10)*(-6.5)*24
    kappa = 1.1
    
    return((psimax-psimin)*(c/mic)**kappa / ((c/mic)**kappa -psimin/psimax))


################################
################################ Extinction time of the resistant
################################

################################ Preliminary definitions
################################ survival probability during treatment
################################ biostatic
surv_bs = []

for i in range(len(c)):
    rhoS = max(bS-a(c[i],0),0)-dS
    rhoR = max(bR-a(c[i],1),0)-dR
    s = rhoR-rhoS
    
    if (bR > a(c[i],1)):
        F = K*s*rhoS*rhoR / (K*s*rhoS*rhoR + gamma*rhoR*N0*(dR+rhoS)*(1-np.exp(-s*TT)) - dR*s*(gamma*N0 - K*rhoS)*(1-np.exp(-rhoR*TT)) )
    else:
        F = np.exp(-dR*TT + np.log(-K*rhoS) - np.log(gamma*N0-gamma*np.exp(TT*rhoS)*N0-K*rhoS) )
    surv_bs.append(max(0,F))

i = surv_bs.index(0.0)
surv_bs[i] = (surv_bs[i-1] + surv_bs[i+1])/2
i = 850
surv_bs[i] = (surv_bs[i-1] + surv_bs[i+1])/2


################################ biocidic
surv_bc = []
for i in range(len(c)):
    rhoS = bS-a(c[i],0)-dS
    rhoR = bR-a(c[i],1)-dR
    s = rhoR-rhoS
    
    F = K*s*rhoS*rhoR / (K*s*rhoS*rhoR + gamma*rhoR*N0*(dR+a(c[i],1) + rhoS)*(1-np.exp(-s*TT)) - (dR+a(c[i],1))*s*(gamma*N0 - K*rhoS)*(1-np.exp(-rhoR*TT)) )
    surv_bc.append(max(0,F))


i = surv_bc.index(0.0)
surv_bc[i] = (surv_bc[i-1] + surv_bc[i+1])/2
i = 850
surv_bc[i] = (surv_bc[i-1] + surv_bc[i+1])/2

################################
################################ Size at time t
################################
dt = 0.01
xS_bs = []
xR_bs = []
xS_bc = []
xR_bc = []


tfin_int = 50
t_aux = np.arange(0,tfin_int,dt)

################################ biostatic
for i in range(len(c)):
    print(i)
    xc = 1/surv_bs[i]
    ### a integral
    aInt = np.zeros(len(t_aux))
    xS_int_aux = np.zeros(len(t_aux))
    #
    xS_int_aux[0] = N0
    aInt[0] = dt * (max(0,bR-a(c[i],1))-xS_int_aux[0]*gamma/K-dR)
    expTime = 0
    #
    for j in range(len(aInt)-1):  
        xSold = xS_int_aux[j]
        xS_int_aux[j+1] = xSold + (max(0,bS-a(c[i],0) )- gamma*xSold/K-dS)*xSold*dt
        aInt[j+1] = aInt[j] + dt*max(0,( max(0,bR-a(c[i],1) )- gamma*xS_int_aux[j+1]/K  - dR))
        expTime += dt * t_aux[j+1] * np.exp(-surv_bs[i]*xc*np.exp(-aInt[j+1])) * surv_bs[i] * xc * np.exp(-aInt[j+1]) * (max(0,max(0,bR-a(c[i],1))-gamma*xS_int_aux[j+1]/K-dR))
    #
    #
    t = min(expTime,TT)
    #
    if (t<TT):
        xR = xc
        xS = xS_int_aux[int(t/dt)]    
        while (t <= TT):
            xS_old = xS
            xR_old = xR
            xR = xR_old + dt * xR_old * (max(0, bR- a(c[i],1))- gamma*(xR_old + xS_old)/K -  dR)
            xS = xS_old + dt * xS_old * (max(0, bS -a(c[i],0))- gamma*(xR_old + xS_old)/K - dS)
            t += dt
    #
    else:
        xR = xc
        xS = N0
        t=0
        while (t <= TT):
            xS_old = xS
            xR_old = xR
            xR = xR_old + dt * xR_old * (max(0, bR- a(c[i],1))- gamma*(xR_old + xS_old)/K -  dR)
            xS = xS_old + dt * xS_old * (max(0,bS-a(c[i],0))- gamma*(xR_old+xS_old)/K - dS)
            t += dt
    #
    xR_surv = 1/surv_bs[i]
    xS_surv = N0
    #
    t = 0
    while (t <= TT):
        xR_surv_old = xR_surv
        xS_surv_old = xS_surv
        xR_surv = xR_surv_old + dt * xR_surv_old * (max(0,bR - a(c[i],1) )- gamma * (xR_surv_old + xS_surv_old)/K - dR)
        xS_surv = xS_surv_old + dt * xS_surv_old * (max(0,bS - a(c[i],0) )- gamma * (xR_surv_old + xS_surv_old)/K - dS)
        t += dt
    #
    if (xR_surv < xR):
        xR_bs.append(max(0,xR_surv))
        xS_bs.append(max(0,xS_surv))
    else:
        xR_bs.append(max(0,xR))
        xS_bs.append(max(0,xS))


################################ biocidal
for i in range(len(c)):
    print(i)
    xc = 1/surv_bc[i]
    ### a integral
    aInt = np.zeros(len(t_aux))
    xS_int_aux = np.zeros(len(t_aux))
    #
    aInt[0] = dt * (max(0,bR-N0*gamma/K)-a(c[i],1)-dR)
    xS_int_aux[0] = N0
    expTime = 0
    #
    for j in range(len(aInt)-1):  
        xSold = xS_int_aux[j]
        xS_int_aux[j+1] = xSold + (bS-gamma*xSold/K-a(c[i],0)-dS)*xSold*dt
        aInt[j+1] = aInt[j] + dt*( max(0,bR- gamma*xS_int_aux[j+1]/K -a(c[i],1) - dR))
        expTime += dt * t_aux[j+1] * np.exp(-surv_bc[i]*xc*np.exp(-aInt[j+1])) * surv_bc[i] * xc * np.exp(-aInt[j+1]) * (max(0,bR-gamma*xS_int_aux[j+1]/K-a(c[i],1)-dR))
    #
    t = min(expTime,TT)
    #print(expTime)
    #
    if (t<TT): 
        xR = xc
        xS = xS_int_aux[int(t/dt)]
        while (t <= TT):
            xS_old = xS
            xR_old = xR
            xR = xR_old + dt * xR_old * (bR- gamma*(xR_old+xS_old)/K - dR-a(c[i],1) )
            xS = xS_old + dt * xS_old * (bS- gamma*(xR_old+xS_old)/K - dS-a(c[i],0) )  
            t += dt
    else:
        xR = xc
        xS = N0
        t=0
        while (t <= TT):
            xS_old = xS
            xR_old = xR
            xR = xR_old + dt * xR_old * (bR- gamma*(xR_old+xS_old)/K - dR-a(c[i],1) )
            xS = xS_old + dt * xS_old * (bS- gamma*(xR_old+xS_old)/K - dS-a(c[i],0) )  
            t += dt   
    #
    xR_surv = 1/surv_bc[i]
    xS_surv = N0
    #
    t = 0
    while (t <= TT):
        xR_surv_old = xR_surv
        xS_surv_old = xS_surv
        xR_surv = xR_surv_old + dt * xR_surv_old * (bR - gamma * (xR_surv_old + xS_surv_old)/K - a(c[i],1) - dR)
        xS_surv = xS_surv_old + dt * xS_surv_old * (bS - gamma * (xR_surv_old + xS_surv_old)/K - a(c[i],0) - dS)
        t += dt
    #
    if (xR_surv < xR):
        xR_bc.append(max(0,xR_surv))
        xS_bc.append(max(0,xS_surv))
    else:
        xR_bc.append(max(0,xR))
        xS_bc.append(max(0,xS))


################################ Defining the integrals
def integrand1(p,initfreq):
    rhoS = bS-dS
    rhoR = bR-dR
    C = rhoS*K*(rhoS-rhoR)/(bS*gamma)
    return(np.exp(C*p) * (np.exp(-C*initfreq) - np.exp(-C)) * (1-np.exp(-C*p))**2 * rhoS / ((1-np.exp(-C)) * (1-np.exp(-C*initfreq)) * (rhoS-rhoR) * p * (1-p) ) ) 

def integrand2(p,initfreq):
    rhoS = bS-dS
    rhoR = bR-dR
    C = rhoS*K*(rhoS-rhoR)/(bS*gamma)
    return(np.exp(C*p) * (np.exp(-C*p) - np.exp(-C)) * (1-np.exp(-C*p)) * rhoS / ((1-np.exp(-C)) * p * (1-p) * (rhoS-rhoR) ) )

################################ Calculating the stochastic estimates
time_bs = []
time_bc = []

for i in range(len(c)):
    rhoS = bS-dS
    
    ############# biostatic
    if (xS_bs[i]>=1):
        #if (xR_bs[i] < 1):
        #    xR_bs[i] = 1
        #
        initfreq_bs = xS_bs[i]/(xS_bs[i]+xR_bs[i])
        #    
        if (initfreq_bs <= 0.0001):
            t1 = 0.
        else:
            t1 = integrate.quad(lambda p: integrand1(p,initfreq_bs),0.0001,initfreq_bs)[0]
        #    
        t2 = integrate.quad(lambda p: integrand2(p,initfreq_bs),initfreq_bs,1.-0.0001)[0]
        #  
        time_bs.append((t1+t2)/rhoS + TT)
    
    elif (xS_bs[i] > 0 and xS_bs[i] < 1 and xR_bs[i] > 0):
        initfreq_bs = 1/(1+xR_bs[i])
        
        if (initfreq_bs <= 0.0001):
            t1 = 0.
        else:
            t1 = integrate.quad(lambda p: integrand1(p,initfreq_bs),0.0001,initfreq_bs)[0]
        
        t2 = integrate.quad(lambda p: integrand2(p,initfreq_bs),initfreq_bs,1.-0.0001)[0]
        
        time_bs.append((t1+t2)/rhoS + TT)
        
    else:
        time_bs.append(0)
            
    ############ biocidic
    if (xS_bc[i]>1):
        if (xR_bc[i] < 1):
            xR_bc[i] = 1
        #
        initfreq_bc = xS_bc[i]/(xS_bc[i]+xR_bc[i])
        
        if (initfreq_bc <= 0.0001):
            t1 = 0.
        else:
            t1 = integrate.quad(lambda p: integrand1(p,initfreq_bc),0.0001,initfreq_bc)[0]
        
        t2 = integrate.quad(lambda p: integrand2(p,initfreq_bc),initfreq_bc,1.-0.0001)[0]
        time_bc.append((t1+t2)/rhoS + TT)
    
    elif (xS_bc[i] > 0 and xS_bc[i] < 1):
        initfreq_bc = 1/(1+xR_bc[i])
        
        if (initfreq_bc <= 0.0001):
            t1 = 0.
        else:
            t1 = integrate.quad(lambda p: integrand1(p,initfreq_bc),0.0001,initfreq_bc)[0]
        
        t2 = integrate.quad(lambda p: integrand2(p,initfreq_bc),initfreq_bc,1.-0.0001)[0]
        time_bc.append((t1+t2)/rhoS + TT)
    
    else:
        time_bc.append(0)

time_bc[169] = (time_bc[168]+time_bc[170])/2




################################
################################ Import Data
################################

################################ biostatic
c_plot = [0.0001, 0.0002, 0.0004, 0.0006, 0.001, 0.002, 0.004, 0.006, 0.01, 0.02, 0.04, 0.06, 0.1,0.2,0.4,0.6]
mean_bs = []
low_bs = []
high_bs = []
for j in range(len(c_plot)):
    my_data = genfromtxt('ext_time_death_c_%f_fac_5_sc_0.txt' % c_plot[j], delimiter=',')
    data = []
    for i in range(len(my_data)):
        data.append(float(my_data[i]))
    data = np.asarray(data)
    mean_bs.append(np.mean(data))
    #low_bs.append(np.percentile(data,5))
    #high_bs.append(np.percentile(data,95))

################################ biocidic
mean_bc = []
low_bc = []
high_bc =[]
for j in range(len(c_plot)):
    my_data = genfromtxt('ext_time_death_c_%f_fac_5_sc_1.txt' % c_plot[j], delimiter=',')
    data = []
    for i in range(len(my_data)):
        data.append(float(my_data[i]))
    data = np.asarray(data)
    mean_bc.append(np.mean(data))


################################
################################ Plot popsizes
################################
plt.semilogx(c,time_bs,color='C0',linewidth=3)
plt.plot(c,time_bc,color='C1',linewidth=3)
plt.plot(c_plot,mean_bs,linestyle='None',marker='s',markersize = '10',color='C0')
plt.plot(c_plot,mean_bc,linestyle='None',marker='s',markersize = '10',color='C1')
#plt.axvspan(c[low], c[up], facecolor='0.2', alpha=0.15)
plt.vlines(c[170],0,65,linestyle='dashed',color='black',lw=2)
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.ylim((0,70))
plt.xlim((0,2))
plt.show()


