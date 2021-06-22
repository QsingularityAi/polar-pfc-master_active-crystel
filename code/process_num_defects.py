#!/usr/bin/python
import sys
import os
import glob
import fnmatch
import csv
import numpy as np
import cyflann
import subprocess
import matplotlib.pyplot as plt
from matplotlib2tikz import save as tikz_save


def cleanup(string):
    return string.replace('.', '_')
# {end cleanup}

out_dir = '/media/Daten/projects/polar_pfc'
setup = 'vop'

V0s = [0.3, 0.31, 0.32, 0.33, 0.35, 0.4, 0.5, 0.55]
data = []
mean_data = []

radius = 80
r = -0.98
for i,V0 in enumerate(V0s):
  postfix = cleanup('R' + str(radius) + '_r' + str(r) + '_v' + str(V0))
  filename = out_dir + '/results/' + setup + '/' + postfix + '/num_defects.csv'
  print i,V0,filename

  mean = 0.0
  n_mean = 0
  data_v0 = []
  with open(filename, 'r') as f:
    k = 0
    for row in f:
      values = row.split()
      num_defects = 0
      total = 0
      for j,value in enumerate(values):
        if j != 6:
          num_defects += int(float(value))
        else:
          total += int(float(value))
      # {end for j}
      if total > 0:
        data_v0.append(num_defects / float(total))
        if k > 1000:
          mean += num_defects / float(total)
          n_mean += 1
        # {end if k}
      # {end if total}
      k += 1

  data.append(data_v0)
  mean_data.append(mean / float(n_mean))

print V0s
print mean_data

plt.plot(V0s, mean_data, '-*')
plt.xlabel('v0')
plt.show()

data_len = min(map(len, data))

data2 = []
for j in range(0, len(V0s)):
  data2_v0 = []
  for i in xrange(0, data_len, 50):
    data2_v0.append(data[j][i])
  data2.append(data2_v0)

print len(data2)

for i,V0 in enumerate(V0s):
  plt.plot(xrange(0,data_len,50), data2[i])

plt.xlabel('time')
plt.legend(map(lambda v0: 'v0=' + str(v0), V0s), loc='upper right')

plt.show()
