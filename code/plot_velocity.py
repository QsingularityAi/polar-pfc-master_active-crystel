#!/usr/bin/python
import sys
import os
import glob
import fnmatch
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib2tikz import save as tikz_save

# R=80
v0 = [0.1, 0.2, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.425, 0.45]
v0_ = ['R80_0.1', 'R80_0.2', 'R80_0.3', 'R80_0.31', 'R80_0.32', 'R80_0.33', 'R80_0.34', 'R80_0.35', 'R80_0.36', 'R80_0.37', 'R80_0.38', 'R80_0.39', 'R80_0.4', 'R80_0.425', 'R80_0.45']
v0_idx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

# R=100
#v0 = [0.1, 0.2, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.425, 0.45, 0.475, 0.8, 0.9]
#v0_ = ['0.1', '0.2', '0.3', '0.31', '0.32', '0.33', '0.34', '0.35', '0.36', '0.37', '0.38', '0.39', '0.4', '0.4_b', '0.425_b', '0.45', '0.45_b', '0.475_b', '0.8', '0.9']
#v0_idx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 13, 14, 14, 15, 16, 17]



mean_velocity = [0]*len(v0)
n_mean_velocity = [0]*len(v0)

for v, idx in zip(v0_, v0_idx):

  timesteps= []
  #polar = []
  #nematic = []
  velocity = []

  with open('time_series_' + str(v) + '.csv','r') as in_file:
    reader = csv.DictReader(in_file, delimiter=',')

    tol=1.0e-2

    for row in reader:
      t = float(row['t'])
      #p = float(row['polar'])
      #n = float(row['nematic'])
      v = float(row['velocity'])

      timesteps.append(t)
      #polar.append(p)
      #nematic.append(n)
      velocity.append(v)

  min_time = 500
  mean_vel = 0.0
  n_vel = 0
  for t,vel in zip(timesteps, velocity):
    if t >= min_time:
      mean_vel += vel
      n_vel += 1

  mean_velocity[idx] += mean_vel / n_vel;
  n_mean_velocity[idx] += 1

for i in range(0, len(mean_velocity)):
  mean_velocity[i] *= 1.0/n_mean_velocity[i]


# plot timeseries

plt.plot(v0, mean_velocity,'-*', v0,v0,'--k')
plt.xlabel('activity $v_0$')
plt.legend(['velocity $\\|\\bar{\\mathbf{v}}(v_0)\\|$', 'activity $v_0$'], loc='lower right')

tikz_save('mean_velocity.tikz',
           figureheight = '\\figureheight',
           figurewidth = '\\figurewidth')
plt.show()
