#!/usr/bin/python
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib2tikz import save as tikz_save

base_dir = '/media/Home/projects/shtns/results/polar_pfc_new'

radius = [10,20,40,60,80] #,100]
V0 = [0.1, 0.2, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.425, 0.45, 0.475, 0.5,0.6,0.7,0.8,0.9,1.0]
r = 60
v0 = 0.3

radius_= []
for r in radius:
  meanN = []
  V0_= []
  for v0 in V0:
    filename = base_dir + '/radius_'+str(r)+'_v0_'+str(v0)+'_dt_0.5/director.csv'
    if not os.path.exists(filename):
      continue
    data = np.genfromtxt(filename, delimiter=',', dtype='float', skip_header=1)
    print data.shape[0]
    n = 0
    num = 0
    for i in range(5000, data.shape[0]):
      n += np.linalg.norm(data[i,1:3])
      num += 1

    if n > 0:
      n /= float(num)
      meanN.append(n)
      V0_.append(v0)


  #plt.plot(radius, meanP, '*-', radius, meanW, '.-')
  plt.plot(V0_, meanN, '*-') #, V0_, meanW, '.-')



plt.xlabel('activity $v_0$')
plt.ylabel('polarity')
plt.legend(['R=10','R=20','R=40','R=60','R=80'], loc='upper left')
plt.show()
