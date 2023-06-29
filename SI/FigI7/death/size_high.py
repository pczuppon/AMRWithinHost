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

################################ Logistic growth dynamics
bS = 2.5
bR = 2.25
dS = 0.5
dR = 0.5

gamma = 1.

K = 1000

N0 = (bS-dS)*K/gamma

################################ Antimicrobial response / Pharmacodynamics
c = np.arange(0,50,0.01)

def a(c,X):
    if (X==0):
        psimax = bS-dS
        mic = 0.67
    if (X==1):
        psimax = bR-dR
        mic = 0.67*20
    psimin = np.log(10)*(-8.1)*24
    kappa = .61
    
    return((psimax-psimin)*(c/mic)**kappa / ((c/mic)**kappa -psimin/psimax))


################################
################################ Surv prob at time t - theory
################################

################################ biostatic
surv_bs = []
survinf_bs = []
TT = 7
for i in range(len(c)):
    rhoS = max(bS-a(c[i],0),0)-dS
    rhoR = max(bR-a(c[i],1),0)-dR
    s = rhoR-rhoS
    
    if (bR > a(c[i],1)):
        F = K*s*rhoS*rhoR / (K*s*rhoS*rhoR + gamma*rhoR*N0*(dR+rhoS)*(1-np.exp(-s*TT)) - dR*s*(gamma*N0 - K*rhoS)*(1-np.exp(-rhoR*TT)) )
        F2 = K*s*rhoS*rhoR / (K*s*rhoS*rhoR + gamma*rhoR*N0*(dR+rhoS)*(1) - dR*s*(gamma*N0 - K*rhoS)*(1) )
    else:
        F = np.exp(-dR*TT + np.log(-K*rhoS) - np.log(gamma*N0-gamma*np.exp(TT*rhoS)*N0-K*rhoS) )
        F2 = 0
    surv_bs.append(max(0,F))
    survinf_bs.append(max(0,F2))

i = 67
surv_bs[i] = (surv_bs[i-1] + surv_bs[i+1])/2
survinf_bs[i] = (survinf_bs[i-1] + survinf_bs[i+1])/2

i = 1340
surv_bs[i] = (surv_bs[i-1] + surv_bs[i+1])/2
survinf_bs[i] = (survinf_bs[i-1] + survinf_bs[i+1])/2


################################ biocidic
surv_bc = []
survinf_bc = []
for i in range(len(c)):
    rhoS = bS-a(c[i],0)-dS
    rhoR = bR-a(c[i],1)-dR
    s = rhoR-rhoS
    
    F = K*s*rhoS*rhoR / (K*s*rhoS*rhoR + gamma*rhoR*N0*(dR+a(c[i],1) + rhoS)*(1-np.exp(-s*TT)) - (dR+a(c[i],1))*s*(gamma*N0 - K*rhoS)*(1-np.exp(-rhoR*TT)) )
    F2 = K*s*rhoS*rhoR / (K*s*rhoS*rhoR + gamma*rhoR*N0*(dR+a(c[i],1) + rhoS)*(1) - (dR+a(c[i],1))*s*(gamma*N0 - K*rhoS)*(1) )
    surv_bc.append(max(0,F))
    survinf_bc.append(max(0,F2))

i = 67
surv_bc[i] = (surv_bc[i-1] + surv_bc[i+1])/2
survinf_bc[i] = (survinf_bc[i-1] + survinf_bc[i+1])/2

i = 1340
surv_bc[i] = (surv_bc[i-1] + surv_bc[i+1])/2
survinf_bc[i] = (survinf_bc[i-1] + survinf_bc[i+1])/2

################################
################################ Size at time t
################################
dt = 0.01
size_bs = []
size_bc = []

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
        xR_surv_old = xR_surv
        xS_surv_old = xS_surv
        xR_surv = xR_surv_old + dt * xR_surv_old * (max(0,bR - a(c[i],1) )- gamma * (xR_surv_old + xS_surv_old)/K - dR)
        xS_surv = xS_surv_old + dt * xS_surv_old * (max(0,bS - a(c[i],0) )- gamma * (xR_surv_old + xS_surv_old)/K - dS)
        t += dt
    #
    size_bs.append(min(xR_surv,xR))

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
        xR_surv_old = xR_surv
        xS_surv_old = xS_surv
        xR_surv = xR_surv_old + dt * xR_surv_old * (bR  - gamma * (xR_surv_old + xS_surv_old)/K - a(c[i],1) - dR)
        xS_surv = xS_surv_old + dt * xS_surv_old * (bS  - gamma * (xR_surv_old + xS_surv_old)/K - a(c[i],0) - dS)
        t += dt
    #
    size_bc.append(min(xR_surv,xR))


################################
################################ Plot surv prob
################################
plt.semilogx(c,size_bs,linewidth=2)
plt.semilogx(c,size_bc,linewidth=2)
plt.vlines(c[67],0,1400,linestyle='dashed',lw=2,color='black')
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.ylim((0,1450))
plt.show()




