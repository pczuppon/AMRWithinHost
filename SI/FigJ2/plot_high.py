import numpy as np
import matplotlib.pyplot as plt
import os
from numpy import genfromtxt

################################
################################ data
################################

os.chdir(os.path.realpath(''))


################################
################################ read data
################################

for j in range(10):
    data = genfromtxt('trajectory_%s_high.txt' % (j), delimiter=',')
    t, s, r, I = [], [], [], []
    #
    for i in range(len(data)):
        t.append(float(data[i][3]))
        s.append(float(data[i][0]))
        r.append(float(data[i][1]))
        I.append(float(data[i][2]))
    #
    plt.plot(t,s,lw=1,color='black',alpha=0.5)
    plt.plot(t,r,lw=1,color='C1',alpha=.5)
    plt.plot(t,I,lw=1,color='C0',alpha=.5)

s = np.asarray(s)
ind = np.where(s >= 100)[0][0]

plt.axvspan(t[ind],t[ind]+ 7, facecolor='0.2', alpha=0.15)
plt.xlim((0,20))
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.ylim((0,400))
plt.show()





