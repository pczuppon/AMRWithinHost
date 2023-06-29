import numpy as np
import matplotlib.pyplot as plt
import os
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

K = 1000

N0 = (bS-dS)*K/gamma

################################ Antimicrobial response / Pharmacodynamics
c = np.arange(0,2000,0.1)

def a(c,X):
    if (X==0):
        psimax = bS-dS
        mic = 20
    if (X==1):
        psimax = bR-dR
        mic = 20*20
    psimin = np.log(10)*(-3)*24
    kappa = 5
    
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

i = 200
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

i = 200
surv_bc[i] = (surv_bc[i-1] + surv_bc[i+1])/2
survinf_bc[i] = (survinf_bc[i-1] + survinf_bc[i+1])/2

################################
################################ Plot surv prob
################################
plt.semilogx(c,surv_bs,linewidth=2)
plt.semilogx(c,surv_bc,linewidth=2)
plt.vlines(c[200],0,.75,linestyle='dashed',lw=2,color='black')
plt.plot(c,survinf_bs,color='C0',linewidth=2,linestyle='dashed')
plt.plot(c,survinf_bc,color='C1',linewidth=2,linestyle='dashed')
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.ylim((-0.025,0.8))
plt.show()



