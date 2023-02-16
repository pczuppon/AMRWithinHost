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

K = 1000

N0 = (bS-dS)*K/gamma

################################ Antimicrobial response / Pharmacodynamics
c = np.arange(0,10,0.0001)

def a(c,X):
    if (X==0):
        psimax = bS-dS
        mic = 0.017
    if (X==1):
        psimax = bR-dR
        mic = 0.017*5
    psimin = np.log(10)*(-6.5) *24.
    kappa = 1.1
    if (X==2):
        psimax = bR-dR
        mic = 0.017*10
    if (X==3):
        psimax = bR-dR
        mic = 0.017*20
    
    return((psimax-psimin)*(c/mic)**kappa / ((c/mic)**kappa -psimin/psimax))


low = np.min(np.where(bS-dS-a(c,0)-(bR-dR-a(c,1))<=0))
low2 = np.min(np.where(bS-dS-a(c,0)-(bR-dR-a(c,2))<=0))
low3 = np.min(np.where(bS-dS-a(c,0)-(bR-dR-a(c,3))<=0))
up = np.min(np.where(bR-dR-a(c,1)<=0))
up2 = np.min(np.where(bR-dR-a(c,2)<=0))
up3 = np.min(np.where(bR-dR-a(c,3)<=0))

################################
################################ Plot popsizes
################################
#plt.semilogx(c,a(c,0),color='black',linewidth=3)
#plt.plot(c,a(c,1),color='C0',linewidth=3)
#plt.semilogx(c,a(c,2),color='C0',linewidth=3,linestyle='dashed')
#plt.plot(c,a(c,3),color='C0',linewidth=3,linestyle='dotted')
#plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
#plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
#plt.show()

    
plt.semilogx(c,bS-dS-a(c,0),color='black',linewidth=3)
plt.semilogx(c,np.maximum(bS-a(c,0),0)-dS,color='black',linewidth=3,linestyle='-.')
#plt.plot(c,np.maximum(bR-a(c,1),0)-dR,color='C0',linewidth=3,linestyle='dotted')
#plt.semilogx(c,np.maximum(bR-a(c,2),0)-dR,color='C0',linewidth=3,linestyle='dashed')
#plt.plot(c,np.maximum(bR-a(c,3),0)-dR,color='C0',linewidth=3)
plt.plot(c,bR-dR-a(c,1),color='C1',linewidth=3,linestyle='dotted')
plt.semilogx(c,bR-dR-a(c,2),color='C1',linewidth=3,linestyle='dashed')
plt.plot(c,bR-dR-a(c,3),color='C1',linewidth=3)
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.plot(c,[0]*len(c),color='black',linestyle='dotted')
#plt.axvspan(c[low], c[up], facecolor='0.2', alpha=0.15)
#plt.axvspan(c[low2], c[up2], facecolor='0.2', alpha=0.15)
#plt.axvspan(c[low3], c[up3], facecolor='0.2', alpha=0.15)
plt.xlim((0,1))
plt.ylim((-15,2.5))
plt.show()

#plt.semilogx(c,bS-dS-a(c,0),color='black',linewidth=3)
#plt.plot(c,bR-dR-a(c,1),color='C0',linewidth=3)
#plt.semilogx(c,bR-dR-a(c,2),color='C0',linewidth=3,linestyle='dashed')
#plt.plot(c,bR-dR-a(c,3),color='C0',linewidth=3,linestyle='dotted')
#plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
#plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
#plt.plot(c,[0]*len(c),color='black',linestyle='dotted')
#plt.xlim((0,0.02))
#plt.ylim((0,2.5))
#plt.show()  




