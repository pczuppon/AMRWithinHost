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
bS = 11.
bR = 10.6
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
        mic = 0.017*10
    psimin = np.log(10)*(-6.5)*24.
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
    
    dt = 0.01
    ############## sensitive strain exponentially declining for all t
    if (bS < a(c[i],0)):
        # time until growth rate of the resistant type becomes positive
        if (bR-a(c[i],1) > 0):
            TT1 = min(TT,max(0,-np.log(((bR-a(c[i],1))*K/(gamma))/N0)/dS))
            # sensitive abundance at TT1
            N0alt = max(0,N0*np.exp(-dS*TT1))
            # compute integral of survival probability numerically
            integral_time = np.arange(TT1,TT,dt)
            integral = 0
            for j in range(len(integral_time)):
                integral += dt*np.exp((integral_time[j]-TT1)*(a(c[i],1)-bR+dR) + gamma*(1-np.exp(-(integral_time[j]-TT1)*dS))*N0alt/(K*dS))*dR
            
            # multiply this survival probability with probability of no death until TT1
            surv_bs.append(np.exp(-dR*TT1)/(1+integral)) #, np.exp(-dR*TT)) )         ##### Maximum because of finite population but TT1 defined on continuity with unreasonably small values (< 1/1000) possible that are not "visible" in the simulations
        else:
            surv_bs.append(np.exp(-dR*TT))
    
    ############## sensitive strain not exponentially declining for all t, just until TTS
    else:
        TTS = min(TT,max(0,-np.log(((bS-a(c[i],0))*K/(gamma))/N0)/dS))          ## time until sensitive strain is exponentially declining
        TTR = min(TT,max(0,-np.log(((bR-a(c[i],1))*K/(gamma))/N0)/dS))          ## time until resistant strain is not reproducing (always smaller than TTS!)
        # sensitive abundance at TTR
        N0atR = max(0,N0*np.exp(-dS*TTR))
        # sensitive abundance at TTS
        N0atS = max(0,N0*np.exp(-dS*TTS))
        
        ### if sensitive strain initially declines exponentially!
        if (TTS > 0):
            # compute integral of survival probability numerically between TTR and TTS
            integral_time = np.arange(TTR,TTS,dt)
            integral = 0
            for j in range(len(integral_time)):
                integral += dt*np.exp((integral_time[j]-TTR)*(a(c[i],1)-bR+dR) + gamma*(1-np.exp(-(integral_time[j]-TTR)*dS))*N0atR/(K*dS))
            
            # Survival between TTR and TTS
            # Notation as in Uecker & Hermisson (Eq. 15)       
            denom = (1+integral*dR)
            A = np.exp( gamma*(1-np.exp(-dS*(TTS-TTR)))*N0atR / (dS*K) + (TTS-TTR)*(dR+a(c[i],1)-bR)) 
            B = 1 - A + integral*dR
            
            offspring = 10
            j=1
            survival = 0
            countprob = 0
            while (j<offspring):
                # probability of extinction of j offspring
                F = (dR*(rhoR*gamma*N0atS*(1-np.exp(-s*(TT-TTS))) - s*(gamma*N0atS - K*rhoS) * (1-np.exp(-rhoR*(TT-TTS)))) / (K*rhoR*rhoS*s + dR*(rhoR*gamma*N0atS*(1-np.exp(-s*(TT-TTS))) - s*(gamma*N0atS - K*rhoS) * (1-np.exp(-rhoR*(TT-TTS))))))**j
                # deriving the probability of j offspring until time TTS (using generating functions) -- see also Mathematica notebook for the derivatives of the probability generating function
                prob_of_j_offspring = A * B**(j-1) / (denom**(j+1))
                # probability of between 1 and 10 offspring
                countprob += prob_of_j_offspring
                # adding the survival probability of j offspring
                survival += (1-F)*prob_of_j_offspring
                j += 1
            
            # total survival probability = surviving until TTR * (survival of 10 offspring until tau + more than 10 offspring (1 - countprob - extinction) )
            surv_bs.append(np.exp(-dR*TTR)* (survival+1-countprob-1-1/(-A-B))) 
        
        ##### if sensititve strain does not decline exponentially
        else:
            F = K*rhoR*rhoS*s / (K*rhoR*rhoS*s + dR*(rhoR*gamma*N0atS*(1-np.exp(-s*(TT))) - s*(gamma*N0atS - K*rhoS) * (1-np.exp(-rhoR*(TT)))))
            surv_bs.append(F)


################################ biocidic
surv_bc = []
for i in range(len(c)):
    rhoS = bS-a(c[i],0)-dS
    rhoR = bR-a(c[i],1)-dR
    s = rhoR-rhoS
    
    F = K*s*rhoS*rhoR / (K*s*rhoS*rhoR + (dR + a(c[i],1)) * (N0 * gamma * rhoR * (1-np.exp(-s*TT)) + s * (rhoS*K - N0*gamma)*(1-np.exp(-rhoR*TT))))
    surv_bc.append(max(0,F))

i = surv_bc.index(0.0)
surv_bc[i] = (surv_bc[i-1] + surv_bc[i+1])/2
i = 1700
surv_bc[i] = (surv_bc[i-1] + surv_bc[i+1])/2

################################
################################ Size at time t
################################
dt = 0.01
xS_bs = []
xR_bs = []
xS_bc = []
xR_bc = []
xS_bs_surv = []
xR_bs_surv = []
xS_bc_surv = []
xR_bc_surv = []

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
    aInt[0] = dt * (max(0,bR-a(c[i],1)-xS_int_aux[0]*gamma/K)-dR)
    expTime = 0
    #
    for j in range(len(aInt)-1):  
        xSold = xS_int_aux[j]
        xS_int_aux[j+1] = xSold + (max(0,bS-a(c[i],0) - gamma*xSold/K)-dS)*xSold*dt
        aInt[j+1] = aInt[j] + dt*max(0,( max(0,bR-a(c[i],1) - gamma*xS_int_aux[j+1]/K)  - dR))
        expTime += dt * t_aux[j+1] * np.exp(-surv_bs[i]*xc*np.exp(-aInt[j+1])) * surv_bs[i] * xc * np.exp(-aInt[j+1]) * (max(0,max(0,bR-a(c[i],1)-gamma*xS_int_aux[j+1]/K)-dR))
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
            xR = xR_old + dt * xR_old * (max(0, bR- a(c[i],1)- gamma*(xR_old + xS_old)/K) -  dR)
            xS = xS_old + dt * xS_old * (max(0, bS -a(c[i],0)- gamma*(xR_old + xS_old)/K) - dS)
            t += dt
    #
    else:
        xR = xc
        xS = N0
        t=0
        while (t <= TT):
            xS_old = xS
            xR_old = xR
            xR = xR_old + dt * xR_old * (max(0, bR- a(c[i],1)- gamma*(xR_old + xS_old)/K) -  dR)
            xS = xS_old + dt * xS_old * (max(0,bS-a(c[i],0)- gamma*(xR_old+xS_old)/K) - dS)
            t += dt
    #
    xR_surv = xc
    xS_surv = N0
    #
    t = 0
    while (t <= TT):
        xR_surv_old = xR_surv
        xS_surv_old = xS_surv
        xR_surv = xR_surv_old + dt * xR_surv_old * (max(0,bR - a(c[i],1) - gamma * (xR_surv_old + xS_surv_old)/K) - dR)
        xS_surv = xS_surv_old + dt * xS_surv_old * (max(0,bS - a(c[i],0) - gamma * (xR_surv_old + xS_surv_old)/K) - dS)
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
    if (surv_bc[i] <= 10**(-10)):
        xR_bc.append(0)
        xS_bc.append(1)
    else:
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
            xS_int_aux[j+1] = xSold + (max(0,bS-gamma*xSold/K)-a(c[i],0)-dS)*xSold*dt
            aInt[j+1] = aInt[j] + dt*( max(0,max(0,bR- gamma*xS_int_aux[j+1]/K) -a(c[i],1) - dR))
            expTime += dt * t_aux[j+1] * np.exp(-surv_bc[i]*xc*np.exp(-aInt[j+1])) * surv_bc[i] * xc * np.exp(-aInt[j+1]) * (max(0,max(0,bR-gamma*xS_int_aux[j+1]/K)-a(c[i],1)-dR))
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
                xR = xR_old + dt * xR_old * (max(0,bR- gamma*(xR_old+xS_old)/K) - dR-a(c[i],1) )
                xS = xS_old + dt * xS_old * (max(0,bS- gamma*(xR_old+xS_old)/K) - dS-a(c[i],0) )  
                t += dt
        else:
            xR = xc
            xS = N0
            t=0
            while (t <= TT):
                xS_old = xS
                xR_old = xR
                xR = xR_old + dt * xR_old * (max(0,bR- gamma*(xR_old+xS_old)/K) - dR-a(c[i],1) )
                xS = xS_old + dt * xS_old * (max(0,bS- gamma*(xR_old+xS_old)/K) - dS-a(c[i],0) )  
                t += dt   
        #
        xR_surv = xc
        xS_surv = N0
        #
        t = 0
        while (t <= TT):
            xR_surv_old = xR_surv
            xS_surv_old = xS_surv
            xR_surv = xR_surv_old + dt * xR_surv_old * (max(0,bR - gamma * (xR_surv_old + xS_surv_old)/K) - dR - a(c[i],1) )
            xS_surv = xS_surv_old + dt * xS_surv_old * (max(0,bS  - gamma * (xR_surv_old + xS_surv_old)/K) - dS - a(c[i],0))
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
    C = rhoS*K*(rhoS-rhoR)/(dS*gamma)
    return((np.exp(-C*(initfreq-p)) - np.exp(-C*(1-p))) * (1-np.exp(-C*p))**2 * rhoS / ((1-np.exp(-C)) * (1-np.exp(-C*initfreq)) * (rhoS-rhoR) * p * (1-p) ) ) 

def integrand2(p,initfreq):
    rhoS = bS-dS
    rhoR = bR-dR
    C = rhoS*K*(rhoS-rhoR)/(dS*gamma)
    return((1 - np.exp(-C*(1-p))) * (1-np.exp(-C*p)) * rhoS / ((1-np.exp(-C)) * p * (1-p) * (rhoS-rhoR) ) )

################################ Calculating the estimates
time_bs = []
time_bc = []

for i in range(len(c)):
    rhoS = bS-dS
    if (xS_bs[i]>=1):
        initfreq_bs = xS_bs[i]/(xS_bs[i]+xR_bs[i])
        #    
        if (initfreq_bs <= 0.0001):
            t1 = 0.
        else:
            t1 = integrate.quad(lambda p: integrand1(p,initfreq_bs),0.0001,initfreq_bs)[0]
            t2 = integrate.quad(lambda p: integrand2(p,initfreq_bs),initfreq_bs,1.-0.0001)[0]
            time_bs.append((t1+t2)/rhoS + TT)
    #
    elif (xS_bs[i] > 0 and xS_bs[i] < 1):
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
    if (xS_bc[i]>=1):
        initfreq_bc = xS_bc[i]/(xS_bc[i]+xR_bc[i])
        #    
        if (initfreq_bc <= 0.0001):
            t1 = 0.
        else:
            t1 = integrate.quad(lambda p: integrand1(p,initfreq_bc),0.0001,initfreq_bc)[0]
            #    
        t2 = integrate.quad(lambda p: integrand2(p,initfreq_bc),initfreq_bc,1.-0.0001)[0]
        time_bc.append((t1+t2)/rhoS+TT)
    #
    elif (xS_bc[i] >= 0 and xS_bc[i] < 1):
        initfreq_bc = 1/(1+xR_bc[i])
        #    
        if (initfreq_bc <= 0.0001):
            t1 = 0.
        else:
            t1 = integrate.quad(lambda p: integrand1(p,initfreq_bc),0.0001,initfreq_bc)[0]
        #
        t2 = integrate.quad(lambda p: integrand2(p,initfreq_bc),initfreq_bc,1.-0.0001)[0]
        time_bc.append((t1+t2)/rhoS+TT)
    #
    else:   # in this case due to numerical inaccurracies the population size goes below zero, which is unrealistic, therefore this code does not apply here, and we change the upper condition to xS_bc[i]>=0
        time_bc.append(0)


time_bs[170] = (time_bs[169]+time_bs[171])/2


################################
################################ Import Data
################################

################################ biostatic
c_plot = [0.0001, 0.0002, 0.0004, 0.0006, 0.001, 0.002, 0.004, 0.006, 0.01, 0.02, 0.04, 0.06, 0.1,0.2,0.4,0.6]
mean_bs = []
low_bs = []
high_bs = []
for j in range(len(c_plot)):
    my_data = genfromtxt('ext_time_birth_c_%f_fac_10_sc_0.txt' % c_plot[j], delimiter=',')
    data = []
    for i in range(len(my_data)):
        data.append(float(my_data[i]))
    data = np.asarray(data)
    mean_bs.append(np.mean(data))
    #low_bs.append(np.percentile(data,5))
    #high_bs.append(np.percentile(data,95))

################################ biocidic
c_plot2 = [0.0001, 0.0002, 0.0004, 0.0006, 0.001, 0.002, 0.004, 0.006, 0.01, 0.02, 0.04, 0.06, 0.1]
mean_bc = []
low_bc = []
high_bc =[]
for j in range(len(c_plot2)):
    my_data = genfromtxt('ext_time_birth_c_%f_fac_10_sc_1.txt' % c_plot[j], delimiter=',')
    data = []
    for i in range(len(my_data)):
        data.append(float(my_data[i]))
    data = np.asarray(data)
    mean_bc.append(np.mean(data))
    #low_bc.append(np.percentile(data,5))
    #high_bc.append(np.percentile(data,95))


low = np.min(np.where(bS-dS-a(c,0)-(bR-dR-a(c,1))<=0))
up = np.min(np.where(bR-dR-a(c,1)<=0))

################################
################################ Plot popsizes
################################
plt.semilogx(c,time_bs,color='C0',linewidth=3)
plt.plot(c,time_bc,color='C1',linewidth=3)
plt.semilogx(c_plot,mean_bs,linestyle='None',marker='^',markersize = '10',color='C0')
plt.plot(c_plot2,mean_bc,linestyle='None',marker='^',markersize = '10',color='C1')
plt.vlines(c[170],0,65,linestyle='dashed',color='black',linewidth=2)
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.ylim((0,70))
plt.show()



