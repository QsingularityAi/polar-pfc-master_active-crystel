#!/usr/bin/python
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib2tikz import save as tikz_save

base_dir = '/media/Home/projects/shtns/results/polar_pfc'

radius = [10,20,40,60,80,100]
V0 = [0.1, 0.2, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.425, 0.45, 0.475, 0.5,0.6,0.7,0.8,0.9,1.0]
r = 60
v0 = 0.3

radius_= []
for r in radius:
  meanP = [0]*len(V0)
  meanW = [0]*len(V0)
  meanP_ = [0]*len(V0)
  meanW_ = [0]*len(V0)
  meanPpsi = [0]*len(V0)
  meanE1 = [0]*len(V0)
  meanE2 = [0]*len(V0)

  nnum = [0]*len(V0)
  for run in range(0,5):
    for i,v0 in enumerate(V0):
      filename = base_dir + '/radius_'+str(r)+'_v0_'+str(v0)+'_run_'+str(run)+'/order_parameters.csv'
      if not os.path.exists(filename):
        continue
      data = np.genfromtxt(filename, delimiter=',', dtype='float')

      if data.shape[0] < 4500:
        print "skipping",filename
        continue

      p = 0
      w = 0
      p_ = 0
      w_ = 0
      ppsi = 0
      e1 = 0
      e2 = 0

      num = 0
      for j in range(500, data.shape[0]):
        p += data[j,1]
        w += data[j,2]
        p_ += np.linalg.norm(data[j,3:5])
        w_ += np.linalg.norm(data[j,6:8])
        ppsi += np.linalg.norm(data[j,9:11])
        e1 += data[j,12]
        e2 += data[j,13]
        num += 1

      if num > 0:
        meanP[i] += p / float(num)
        meanW[i] += w / float(num)
        meanP_[i] += p_ / float(num)
        meanW_[i] += w_ / float(num)
        meanPpsi[i] += ppsi / float(num)
        meanE1[i] += e1 / float(num)
        meanE2[i] += e2 / float(num)

        nnum[i] += 1
      # {end if}
    # {end for V0}
  # {end for run}

  V0_ = []
  P = []
  W = []
  P_ = []
  W_ = []
  Ppsi = []
  E1 = []
  E2 = []
  for i in range(0,len(V0)):
    if nnum[i] > 0:
      V0_.append(V0[i])
      P.append(meanP[i] / float(nnum[i]))
      W.append(meanW[i] / float(nnum[i]))
      P_.append(meanP_[i] / float(nnum[i]))
      W_.append(meanW_[i] / float(nnum[i]))
      Ppsi.append(meanPpsi[i] / float(nnum[i]))
      E1.append(meanE1[i] / float(nnum[i]))
      E2.append(meanE2[i] / float(nnum[i]))


  plt.plot(V0_, P, '*-')

  with open(base_dir + '/data_radius_'+str(r)+'.csv','w') as out_file:
    writer = csv.writer(out_file, delimiter=' ')
    for data in zip(V0_, P, W, P_, W_, Ppsi, E1, E2):
      writer.writerow(map(str, data))


plt.xlabel('activity $v_0$')
plt.ylabel('polarity')
plt.legend(['R=10','R=20','R=40','R=60','R=80','R=100'], loc='upper left')
plt.show()
