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

################################ long
c_plot = [0.0001, 0.0002, 0.0004, 0.0006, 0.001, 0.002, 0.004, 0.006, 0.01, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.6]
mean_bs = []
low_bs = []
high_bs = []
for j in range(len(c_plot)):
    my_data = genfromtxt('popsize_death_c_%f_fac_10_sc_0.txt' % c_plot[j], delimiter=',')
    data = []
    for i in range(len(my_data)):
        data.append(float(my_data[i][1]))
    data = np.asarray(data)
    mean_bs.append(np.mean(data))
    low_bs.append(np.percentile(data,5))
    high_bs.append(np.percentile(data,95))

c_plot2 = [0.0001, 0.0002, 0.0004, 0.0006, 0.001, 0.002, 0.004, 0.006, 0.01, 0.02, 0.04, 0.06, 0.1, 0.2]
mean_bc = []
low_bc = []
high_bc = []
for j in range(len(c_plot2)):
    my_data = genfromtxt('popsize_death_c_%f_fac_10_sc_1.txt' % c_plot2[j], delimiter=',')
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
        mic = 0.017*10
    psimin = np.log(10)*(-6.5)
    kappa = 1.1
    
    return((psimax-psimin)*(c/mic)**kappa / ((c/mic)**kappa -psimin/psimax))

################################
################################ Survival probability
################################
TT = 7

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
    aInt[0] = dt * (max(0,max(0,bR-a(c[i],1))-N0*gamma/K-dR))
    xS_int_aux[0] = N0
    expTime = 0
    #
    for j in range(len(aInt)-1):  
        xSold = xS_int_aux[j]
        xS_int_aux[j+1] = xSold + (max(0,bS-a(c[i],0))-gamma*xSold/K-dS)*xSold*dt
        aInt[j+1] = aInt[j] + dt*( max(0,max(0,bR-a(c[i],1))- gamma*xS_int_aux[j+1]/K  - dR))
        expTime += dt * t_aux[j+1] * np.exp(-surv_bs[i]*xc*np.exp(-aInt[j+1])) * surv_bs[i] * xc * np.exp(-aInt[j+1]) * (max(0,max(0,bR-a(c[i],1))-gamma*xS_int_aux[j+1]/K-dR))
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
        xR = max(0,xR_old + dt * xR_old * (max(0,bR-a(c[i],1)) - dR - gamma*(xR_old+xS_old)/K))
        xS = max(0,xS_old + dt * xS_old * (max(0,bS-a(c[i],0)) - dS - gamma*(xR_old+xS_old)/K))
        t += dt
    #
    t=0
    while (t<= TT):
        xS_old = xS_det
        xR_old = xR_det
        xR_det = max(0,xR_old + dt * xR_old * (max(0,bR-a(c[i],1)) - dR - gamma*(xR_old+xS_old)/K))
        xS_det = max(0,xS_old + dt * xS_old * (max(0,bS-a(c[i],0)) - dS - gamma*(xR_old+xS_old)/K))
        xR_surv_old = xR_surv
        xS_surv_old = xS_surv
        xR_surv = xR_surv_old + dt * xR_surv_old * (max(0,bR - a(c[i],1) )- gamma * (xR_surv_old + xS_surv_old)/K - dR)
        xS_surv = xS_surv_old + dt * xS_surv_old * (max(0,bS - a(c[i],0) )- gamma * (xR_surv_old + xS_surv_old)/K - dS)
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
    aInt[0] = dt * max(0,(bR-N0*gamma/K-a(c[i],1)-dR))
    xS_int_aux[0] = N0
    expTime = 0
    #
    for j in range(len(aInt)-1):  
        xSold = xS_int_aux[j]
        xS_int_aux[j+1] = xSold + (bS-gamma*xSold/K-a(c[i],0)-dS)*xSold*dt
        aInt[j+1] = aInt[j] + dt*( max(0,bR- gamma*xS_int_aux[j+1]/K -a(c[i],1) - dR))
        expTime += dt * t_aux[j+1] * np.exp(-surv_bc[i]*xc*np.exp(-aInt[j+1])) * surv_bc[i] * xc * np.exp(-aInt[j+1]) * max(0,(bR-gamma*xS_int_aux[j+1]/K-a(c[i],1)-dR))
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
        xR = max(0,xR_old + dt * xR_old * (bR-dR-a(c[i],1) - gamma*(xR_old+xS_old)/K))
        xS = max(0,xS_old + dt * xS_old * (bS-dS-a(c[i],0) - gamma*(xR_old+xS_old)/K))
        t += dt
    #
    t=0
    while (t <= TT):    
        xS_old = xS_det
        xR_old = xR_det
        xR_det = max(0,xR_old + dt * xR_old * (bR-dR-a(c[i],1) - gamma*(xR_old+xS_old)/K))
        xS_det = max(0,xS_old + dt * xS_old * (bS-dS-a(c[i],0) - gamma*(xR_old+xS_old)/K))
        xR_surv_old = xR_surv
        xS_surv_old = xS_surv
        xR_surv = xR_surv_old + dt * xR_surv_old * (bR  - gamma * (xR_surv_old + xS_surv_old)/K - a(c[i],1) - dR)
        xS_surv = xS_surv_old + dt * xS_surv_old * (bS  - gamma * (xR_surv_old + xS_surv_old)/K - a(c[i],0) - dS)
        t += dt
    #
    size_bc.append(min(xR_surv,xR))
    size_bc_det.append(xR_det)


low = np.min(np.where(bS-dS-a(c,0)-(bR-dR-a(c,1))<=0))
up = np.min(np.where(bR-dR-a(c,1)<=0))

################################
################################ Plot popsizes
################################
plt.semilogx(c,size_bs,color='C0',linewidth=3)
plt.semilogx(c,size_bc,color='C1',linewidth=3)
plt.plot(c,size_bs_det,color='C0',linewidth=3,linestyle='dotted',alpha=.5)
plt.plot(c,size_bc_det,color='C1',linewidth=3,linestyle='dotted',alpha=.5)
plt.plot(c_plot,mean_bs,linestyle='None',marker='^',markersize = '10',color='C0')
plt.plot(c_plot2,mean_bc,linestyle='None',marker='^',markersize = '10',color='C1')
plt.vlines(c[170],0,1550,linestyle='dashed',color='black',lw = 2)
#plt.axvspan(c[low], c[up], facecolor='0.2', alpha=0.15)
#plt.fill_between(c_plot,low_bs,high_bs,alpha=0.25,color='C0')
#plt.fill_between(c_plot,low_bc,high_bc,alpha=0.25,color='C1')
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.ylim((0,1600))
plt.show()




