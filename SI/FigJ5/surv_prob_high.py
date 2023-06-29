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
################################ Import Data
################################

################################ Short - bacteriostatic
my_data = genfromtxt('survIR_day_vary_c_factor_0.6.txt', delimiter=',')

data1 = []
c_data1 = []
for i in range(len(my_data)):
    data1.append(float(my_data[i][1]))
    c_data1.append(float(my_data[i][0]))


################################ Short - bacteriostatic
my_data = genfromtxt('survIR_day_vary_c_factor_0.6_alt.txt', delimiter=',')

data2 = []
c_data2 = []
for i in range(len(my_data)):
    data2.append(float(my_data[i][1]))
    c_data2.append(float(my_data[i][0]))

################################
################################ Compute MIC
################################
gamma = 0.01
hS = 0.3
bS = 0.6
c= np.arange(0,1,0.001)

repl = bS*(1-np.tanh(15*(c-0.3)))
growth = repl-gamma
ind = np.where(growth<= 0)[0][0]

################################
################################ Plot surv prob
################################
plt.plot(c_data1,data1,'o',color='C0',markersize=10)
plt.vlines(c[ind],0,.75,linestyle='dashed',color='black')
plt.semilogy(c_data2,data2,'o',color='C1',markersize=10)
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.ylim((-0.025,1.))
plt.show()



