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

K = 10000

N0 = (bS-dS)*K/gamma

################################ Antimicrobial response / Pharmacodynamics -- Ampicillin
c = np.arange(0,1000,0.01)

def a(c,X):
    if (X==0):
        psimax = bS-dS
        mic = 3.4
    if (X==1):
        psimax = bR-dR
        mic = 3.4*5
    psimin = np.log(10)*(-4.)*24
    kappa = 0.75
    if (X==2):
        psimax = bR-dR
        mic = 3.4*10
    if (X==3):
        psimax = bR-dR
        mic = 3.4*20
    
    return((psimax-psimin)*(c/mic)**kappa / ((c/mic)**kappa -psimin/psimax))



plt.semilogx(c,bS-dS-a(c,0),color='black',linewidth=3)
plt.plot(c,bR-dR-a(c,1),color='C0',linewidth=3)
plt.plot(c,bR-dR-a(c,2),color='C0',linewidth=3,linestyle='dashed')
plt.plot(c,bR-dR-a(c,3),color='C0',linewidth=3,linestyle='dotted')
plt.plot(c,[0]*len(c),color='black',ls='dotted')
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.xlim((0,1000))
plt.ylim((-10,2.5))
plt.show()




