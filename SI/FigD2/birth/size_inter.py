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
c_plot = [0.0001, 0.0002, 0.0004, 0.0006, 0.001, 0.002, 0.004, 0.006, 0.01, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4]
mean_bs = []
low_bs = []
high_bs = []
for j in range(len(c_plot)):
    my_data = genfromtxt('popsize_freq_birth_c_%f_fac_10_sc_0.txt' % c_plot[j], delimiter=',')
    data = []
    for i in range(len(my_data)):
        data.append(float(my_data[i][1])/(float(my_data[i][0])+float(my_data[i][1])))
    data = np.asarray(data)
    mean_bs.append(np.mean(data))
    low_bs.append(np.percentile(data,5))
    high_bs.append(np.percentile(data,95))

c_plot2 = [0.0001, 0.0002, 0.0004, 0.0006, 0.001, 0.002, 0.004, 0.006, 0.01, 0.02, 0.04, 0.06, 0.1, 0.2]
mean_bc = []
low_bc = []
high_bc = []
for j in range(len(c_plot2)):
    my_data = genfromtxt('popsize_freq_birth_c_%f_fac_10_sc_1.txt' % c_plot2[j], delimiter=',')
    data = []
    d2 = []
    d3 = []
    for i in range(len(my_data)):
        data.append(float(my_data[i][1])/(float(my_data[i][0])+float(my_data[i][1])))
    data = np.asarray(data)
    mean_bc.append(np.mean(data))
    low_bc.append(np.percentile(data,5))
    high_bc.append(np.percentile(data,95))


################################
################################ Plot popsizes
################################
plt.semilogx(c_plot,mean_bs,linestyle='None',marker='^',markersize = '10',color='C0')
plt.plot(c_plot2,mean_bc,linestyle='None',marker='^',markersize = '10',color='C1')
plt.vlines(0.017,-0.01,.95,linestyle='dashed',color='black')
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
#plt.ylim((0,1700))
plt.show()




