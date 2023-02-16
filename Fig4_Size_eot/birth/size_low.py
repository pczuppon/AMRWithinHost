import numpy as np
import matplotlib.pyplot as plt
import os
from numpy import genfromtxt

################################
################################ data
################################

os.chdir(os.path.realpath(''))

################################
################################ Import Data - popsizes
################################
################################ 
c_plot = [0.000100, 0.000200, 0.000400, 0.000600, 0.001000, 0.002000, 0.004000, 0.006000, 0.010000, 0.020000, 0.040000, 0.060000, 0.100000, 0.200000, 0.400000, 0.600000]
mean_bs = []
low_bs = []
high_bs = []
for j in range(len(c_plot)):
    my_data = genfromtxt('popsize_birth_c_%f_fac_5_sc_0.txt' % c_plot[j], delimiter=',')
    data = []
    for i in range(len(my_data)):
        data.append(float(my_data[i][1]))
    data = np.asarray(data)
    mean_bs.append(np.mean(data))
    low_bs.append(np.percentile(data,5))
    high_bs.append(np.percentile(data,95))

c_plot2 = [0.000100, 0.000200, 0.000400, 0.000600, 0.001000, 0.002000, 0.004000, 0.006000, 0.010000, 0.020000, 0.040000, 0.060000, 0.100000]
mean_bc = []
low_bc = []
high_bc = []
for j in range(len(c_plot2)):
    my_data = genfromtxt('popsize_birth_c_%f_fac_5_sc_1.txt' % c_plot2[j], delimiter=',')
    data = []
    for i in range(len(my_data)):
        data.append(float(my_data[i][1]))
    data = np.asarray(data)
    mean_bc.append(np.mean(data))
    low_bc.append(np.percentile(data,5))
    high_bc.append(np.percentile(data,95))


################################
################################ Parameter definitions
################################

################################ Logistic growth dynamics
bS = 2.5
bR = 2.25
dS = 0.5
dR = 0.5

gamma = 1.

K = 1000

N0 = (bS-dS)*K/gamma

################################ Antimicrobial response / Pharmacodynamics
c = np.arange(0,1,0.0001)

def a(c,X):
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
################################ Survival probability
################################
TT = 7
TTinf = 100

################################ biostatic
surv_bs = []
survinf_bs = []
N0 = (bS-dS)*K/gamma
for i in range(len(c)):
    print(i)
    rhoS = max(bS-a(c[i],0),0)-dS   
    rhoR = max(bR-a(c[i],1),0)-dR   
    
    s = rhoR-rhoS
    
    dt = 0.01
    ############## sensitive strain exponentially declining for all t
    if (bS < a(c[i],0)):
        # time until growth rate of the resistant type becomes positive
        if (bR-a(c[i],1) > 0):
            TT1 = min(TT,max(0,-np.log(((bR-a(c[i],1))*K/(gamma))/N0)/dS))
            TT1inf = max(0,-np.log(((bR-a(c[i],1))*K/(gamma))/N0)/dS)
            # sensitive abundance at TT1
            N0alt = max(0,N0*np.exp(-dS*TT1))
            N0altinf = max(0,N0*np.exp(-dS*TT1inf))
            # compute integral of survival probability numerically
            integral_time = np.arange(TT1,TT,dt)
            integral = 0
            for j in range(len(integral_time)):
                integral += dt*np.exp((integral_time[j]-TT1)*(a(c[i],1)-bR+dR) + gamma*(1-np.exp(-(integral_time[j]-TT1)*dS))*N0alt/(K*dS))*dR
            
            # multiply this survival probability with probability of no death until TT1
            surv_bs.append( np.exp(-dR*TT1)/(1+integral)) 
            
            # compute integral of survival probability numerically
            integral_time = np.arange(TT1inf,TTinf,dt)
            integral = 0
            for j in range(len(integral_time)):
                integral += dt*np.exp((integral_time[j]-TT1inf)*(a(c[i],1)-bR+dR) + gamma*(1-np.exp(-(integral_time[j]-TT1inf)*dS))*N0altinf/(K*dS))*dR
            survinf_bs.append( max(np.exp(-dR*TT1inf)/(1+integral),0))
        
        else:
            surv_bs.append(np.exp(-dR*TT))
            survinf_bs.append(0)
    
    ############## sensitive strain not exponentially declining for all t, just until TTS
    else:
        TTS = min(TT,max(0,-np.log(((bS-a(c[i],0))*K/(gamma))/N0)/dS))          ## time until sensitive strain is exponentially declining
        TTR = min(TT,max(0,-np.log(((bR-a(c[i],1))*K/(gamma))/N0)/dS))          ## time until resistant strain is not reproducing (always smaller than TTS!)
        # sensitive abundance at TTR
        N0atR = max(0,N0*np.exp(-dS*TTR))
        # sensitive abundance at TTS
        N0atS = max(0,N0*np.exp(-dS*TTS))
        
        # infinite treatment scenario:
        TTSinf = max(0,-np.log(((bS-a(c[i],0))*K/(gamma))/N0)/dS)
        TTRinf = max(0,-np.log(((bR-a(c[i],1))*K/(gamma))/N0)/dS)
        N0atRinf = max(0,N0*np.exp(-dS*TTRinf))
        N0atSinf = max(0,N0*np.exp(-dS*TTSinf))
        
        ### if sensitive strain initially declines exponentially!
        if (TTS > 0):
            
            # compute integral of survival probability numerically between TTR and TTS
            integral_time = np.arange(TTR,TTS,dt)
            integral = 0
            for j in range(len(integral_time)):
                integral += dt*np.exp((integral_time[j]-TTR)*(a(c[i],1)-bR+dR) + gamma*(1-np.exp(-(integral_time[j]-TTR)*dS))*N0atRinf/(K*dS))
                   
            denom = (1+integral*dR)
            A = np.exp( gamma*(1-np.exp(-dS*(TTS-TTR)))*N0atRinf / (dS*K) + (TTS-TTR)*(dR+a(c[i],1)-bR)) 
            B = 1 - A + integral*dR
            
            
            offspring = 10
            j=1
            survival = 0
            countprob = 0
            while (j<offspring):
                F = F = (dR*(rhoR*gamma*N0atSinf*(1-np.exp(-s*(TTinf-TTS))) - s*(gamma*N0atSinf - K*rhoS) * (1-np.exp(-rhoR*(TTinf-TTS)))) / (K*rhoR*rhoS*s + dR*(rhoR*gamma*N0atSinf*(1-np.exp(-s*(TTinf-TTS))) - s*(gamma*N0atSinf - K*rhoS) * (1-np.exp(-rhoR*(TTinf-TTS))))))**j
                prob = A * B**(j-1) / (denom**(j+1))
                countprob += prob
                survival += (1-F)*prob
                j += 1
            
            survinf_bs.append(max(0,min(1,np.exp(-dR*TTR)* (survival+1-countprob-1-1/(-A-B)))))
            
            
            if (TTS < TT):
                # compute integral of survival probability numerically between TTR and TTS
                integral_time = np.arange(TTR,TTS,dt)
                integral = 0
                for j in range(len(integral_time)):
                    integral += dt*np.exp((integral_time[j]-TTR)*(a(c[i],1)-bR+dR) + gamma*(1-np.exp(-(integral_time[j]-TTR)*dS))*N0atR/(K*dS))
                       
                denom = (1+integral*dR)
                A = np.exp( gamma*(1-np.exp(-dS*(TTS-TTR)))*N0atR / (dS*K) + (TTS-TTR)*(dR+a(c[i],1)-bR)) 
                B = 1 - A + integral*dR
                
                offspring = 10
                j=1
                survival = 0
                countprob = 0
                while (j<offspring):
                    F = (dR*(rhoR*gamma*N0atS*(1-np.exp(-s*(TT-TTS))) - s*(gamma*N0atS - K*rhoS) * (1-np.exp(-rhoR*(TT-TTS)))) / (K*rhoR*rhoS*s + dR*(rhoR*gamma*N0atS*(1-np.exp(-s*(TT-TTS))) - s*(gamma*N0atS - K*rhoS) * (1-np.exp(-rhoR*(TT-TTS))))))**j
                    prob = A * B**(j-1) / (denom**(j+1))
                    countprob += prob
                    survival += (1-F)*prob
                    j += 1
                
                surv_bs.append(np.exp(-dR*TTR)* (survival+1-countprob-1-1/(-A-B))) 
                
            else:
                integral_time = np.arange(TTR,TT,dt)
                integral = 0
                for j in range(len(integral_time)):
                    integral += dt*np.exp((integral_time[j]-TTR)*(a(c[i],1)-bR+dR) + gamma*(1-np.exp(-(integral_time[j]-TTR)*dS))*N0atR/(K*dS))
                
                surv_bs.append(np.exp(-dR*TTR)/(1+integral*dR))
        
        ##### if sensititve strain does not decline exponentially
        else:
            F = K*rhoR*rhoS*s / (K*rhoR*rhoS*s + dR*(rhoR*gamma*N0atS*(1-np.exp(-s*(TT))) - s*(gamma*N0atS - K*rhoS) * (1-np.exp(-rhoR*(TT)))))
            F2 = K*rhoR*rhoS*s / (K*rhoR*rhoS*s + dR*(rhoR*gamma*N0atS*(1-np.exp(-s*(TTinf))) - s*(gamma*N0atS - K*rhoS) * (1-np.exp(-rhoR*(TTinf)))))
            surv_bs.append(F)
            survinf_bs.append(max(0,F2))

surv_bs[170] = (surv_bs[169] + surv_bs[171])/2

############# bacteriocidic
surv_bc = []
survinf_bc = []
N0 = (bS-dS)*K/gamma
for i in range(len(c)):
    rhoS = bS-a(c[i],0)-dS
    rhoR = bR-a(c[i],1)-dR
    s = rhoR-rhoS
    
    F = K*s*rhoS*rhoR / (K*s*rhoS*rhoR + (dR + a(c[i],1)) * (N0 * gamma * rhoR * (1-np.exp(-s*TT)) + s * (rhoS*K - N0*gamma)*(1-np.exp(-rhoR*TT))))
    F2 = K*s*rhoS*rhoR / (K*s*rhoS*rhoR + (dR + a(c[i],1)) * (N0 * gamma * rhoR * (1-np.exp(-s*TTinf)) + s * (rhoS*K - N0*gamma)*(1-np.exp(-rhoR*TTinf))))
    surv_bc.append(max(0,F))
    survinf_bc.append(max(0,F2))

i = 170
surv_bc[i] = (surv_bc[i-1] + surv_bc[i+1])/2
i = 850
surv_bc[i] = (surv_bc[i-1] + surv_bc[i+1])/2

################################
################################ Size at time t
################################
dt = 0.01
size_bs = []
size_bs_det = []
size_bc = []
size_bc_det = []

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
    aInt[0] = dt * (max(0,max(0,bR-a(c[i],1)-N0*gamma/K)-dR))
    xS_int_aux[0] = N0
    expTime = 0
    #
    for j in range(len(aInt)-1):  
        xSold = xS_int_aux[j]
        xS_int_aux[j+1] = xSold + (max(0,bS-a(c[i],0)-gamma*xSold/K)-dS)*xSold*dt
        aInt[j+1] = aInt[j] + dt*( max(0,max(0,bR-a(c[i],1) - gamma*xS_int_aux[j+1]/K)  - dR))
        expTime += dt * t_aux[j+1] * np.exp(-surv_bs[i]*xc*np.exp(-aInt[j+1])) * surv_bs[i] * xc * np.exp(-aInt[j+1]) * (max(0,max(0,bR-a(c[i],1)-gamma*xS_int_aux[j+1]/K)-dR))
    #
    t = min(expTime,TT)
    #
    xR = xc
    xS = xS_int_aux[int(t/dt)]
    xR_det = 1
    xS_det = N0
    xR_surv = xc
    xS_surv = N0
    #
    while (t <= TT):
        xS_old = xS
        xR_old = xR
        xR = xR_old + dt * xR_old * (max(0, bR- a(c[i],1)- gamma*(xR_old + xS_old)/K) -  dR)
        xS = xS_old + dt * xS_old * (max(0,bS-a(c[i],0)- gamma*(xR_old+xS_old)/K) - dS)
        t += dt
    #
    t = 0
    while (t <= TT):
        xS_old = xS_det
        xR_old = xR_det
        xR_det = xR_old + dt * xR_old * (max(0,bR-a(c[i],1)- gamma*(xR_old+xS_old)/K) - dR)
        xS_det = xS_old + dt * xS_old * (max(0,bS-a(c[i],0)- gamma*(xR_old+xS_old)/K) - dS)   
        xR_surv_old = xR_surv
        xS_surv_old = xS_surv
        xR_surv = xR_surv_old + dt * xR_surv_old * (max(0,bR - a(c[i],1) - gamma * (xR_surv_old + xS_surv_old)/K) - dR)
        xS_surv = xS_surv_old + dt * xS_surv_old * (max(0,bS - a(c[i],0) - gamma * (xR_surv_old + xS_surv_old)/K) - dS)
        t += dt
    #
    size_bs.append(min(xR_surv,xR))
    size_bs_det.append(xR_det)

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
        xS_int_aux[j+1] = xSold + (max(0,bS-gamma*xSold/K)-a(c[i],0)-dS)*xSold*dt
        aInt[j+1] = aInt[j] + dt*( max(0,max(0,bR- gamma*xS_int_aux[j+1]/K) -a(c[i],1) - dR))
        expTime += dt * t_aux[j+1] * np.exp(-surv_bc[i]*xc*np.exp(-aInt[j+1])) * surv_bc[i] * xc * np.exp(-aInt[j+1]) * (max(0,max(0,bR-gamma*xS_int_aux[j+1]/K)-a(c[i],1)-dR))
    #
    t = min(expTime,TT)
    #
    xR = xc    
    xS = xS_int_aux[int(t/dt)]
    xR_det = 1
    xS_det = N0
    xR_surv = xc
    xS_surv = N0
    #
    while (t <= TT):
        xS_old = xS
        xR_old = xR
        xR = xR_old + dt * xR_old * (max(0,bR- gamma*(xR_old+xS_old)/K) -dR-a(c[i],1) )
        xS = xS_old + dt * xS_old * (max(0,bS- gamma*(xR_old+xS_old)/K) - dS-a(c[i],0) )  
        t += dt
    #
    t = 0
    while (t <= TT):
        xS_old = xS_det
        xR_old = xR_det
        xR_det = xR_old + dt * xR_old * (max(0,bR- gamma*(xR_old+xS_old)/K) -dR-a(c[i],1) )
        xS_det = xS_old + dt * xS_old * (max(0,bS- gamma*(xR_old+xS_old)/K) - dS-a(c[i],0) )  
        xR_surv_old = xR_surv
        xS_surv_old = xS_surv
        xR_surv = xR_surv_old + dt * xR_surv_old * (max(0,bR  - gamma * (xR_surv_old + xS_surv_old)/K) - a(c[i],1) - dR)
        xS_surv = xS_surv_old + dt * xS_surv_old * (max(0,bS  - gamma * (xR_surv_old + xS_surv_old)/K) - a(c[i],0) - dS)
        t += dt
    #
    size_bc.append(min(xR_surv,xR))
    size_bc_det.append(xR_det)


low = np.min(np.where(bS-dS-a(c,0)-(bR-dR-a(c,1))<=0))
up = np.min(np.where(bR-dR-a(c,1)<=0))
#mid = np.min(np.where(bS-a(c,0) <= 0))


################################
################################ Plot popsizes
################################
plt.semilogx(c,size_bs,color='C0',linewidth=3)
plt.semilogx(c,size_bc,color='C1',linewidth=3)
plt.plot(c,size_bs_det,color='C0',linewidth=3,linestyle='dotted',alpha=0.5)
plt.plot(c,size_bc_det,color='C1',linewidth=3,linestyle='dotted',alpha=0.5)
plt.plot(c_plot,mean_bs,linestyle='None',marker='s',markersize = '10',color='C0')
plt.plot(c_plot2,mean_bc,linestyle='None',marker='s',markersize = '10',color='C1')
plt.vlines(c[170],0,1550,linestyle='dashed',color='black',lw=2)
#plt.vlines(c[mid],0,1550,linestyle='dashed',color='red',lw=2)
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.ylim((0,1600))
plt.show()

