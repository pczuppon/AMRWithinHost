import numpy as np
import matplotlib.pyplot as plt

################################
################################ Parameter definitions
################################

################################ Logistic growth dynamics
bS = 2.5
bR = 2.25
dS = 0.5
dR = 0.5

gamma = 1.

gamma2 = 1./dS
gamma1 = 1./bS

K = 1000


################################ Antimicrobial response / Pharmacodynamics
c = 0.001

def a(c,X):
    if (X==0):
        psimax = bS-dS
        mic = 0.017
    if (X==1):
        psimax = bR-dR
        mic = 0.017*20
    psimin = np.log(10)*(-6.5)*24
    kappa = 1.1
    
    return((psimax-psimin)*(c/mic)**kappa / ((c/mic)**kappa -psimin/psimax))

################################
################################ Survival probability
################################
TT = 7

################################
################################ Size at time t
################################
dt = 0.01
t = np.arange(0,TT+2*dt,dt)
N0 = (bS-dS)*K/(gamma)
size_bs_birth = [N0]
size_bc_birth = [N0]
size_bs_death = [N0]
size_bc_death = [N0]

t0 = 0
xS = N0
while (t0 < TT):
    xS_old = size_bs_birth[-1]
    xS = max(0,xS_old + dt * xS_old * (max(0,(bS-a(c,0))*(1- gamma1*xS_old/K)) - dS ))
    t0 += dt
    size_bs_birth.append(max(0,xS))

t0 = 0
xS = N0
while (t0 < TT):
    xS_old = size_bc_birth[-1]
    xS = max(0,xS_old + dt * xS_old * (max(0,(bS*(1- gamma1*xS_old/K))) - dS - a(c,0)))
    t0 += dt
    size_bc_birth.append(max(0,xS))

t0 = 0
xS = N0
while (t0 < TT):
    xS_old = size_bs_death[-1]
    xS = max(0,xS_old + dt * xS_old * (max(0,(bS-a(c,0)))- dS*(1+gamma2*xS_old/K) ))
    t0 += dt
    size_bs_death.append(max(0,xS))

t0 = 0
xS = N0
while (t0 < TT):
    xS_old = size_bc_death[-1]
    xS = max(0,xS_old + dt * xS_old * (bS- (1+gamma2*xS_old/K)*(dS + a(c,0))))
    t0 += dt
    size_bc_death.append(max(0,xS))

################################
################################ Plot popsizes
################################
plt.plot(t,size_bs_birth,color='C0',linewidth=3)
plt.plot(t,size_bc_birth,color='C1',linewidth=3,linestyle='dashed')
plt.plot(t,size_bs_death,color='C0',linewidth=3,linestyle='-.')
plt.plot(t,size_bc_death,color='C1',linewidth=3,linestyle='dotted')
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.ylim((0,2200))
plt.show()




