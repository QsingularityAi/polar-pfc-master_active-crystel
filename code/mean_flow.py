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

if len(sys.argv) > 2:
  print "usage: ./mean_flow.py postfix"
  exit()

if len(sys.argv) > 1:
  postfix = '_' + sys.argv[1]
else:
  postfix = ''


filename_pattern = 'particles_*.csv'
base_dir_ = '/media/Home/projects/shtns/results/polar_pfc_new'

dt_write = 0.5
radius = 80
v0 = 0.32

base_dir = base_dir_ + '/radius_' + str(radius) + '_v0_' + str(v0) + '_dt_' + str(dt_write)
out_dir = base_dir + "/mean_flow"

if not os.path.exists(out_dir):
  os.makedirs(out_dir)

timesteps = np.arange(1000, 9100, 100)
dt = 1000

t0 = 0
t1 = 10000

filename = base_dir + '/pca_vel_' + str(t0) + '-' + str(t1) + '.csv'

if not os.path.exists(filename):
  subprocess.call('./generate_mesh_velocity.py \'' + base_dir + '/' + filename_pattern + '\' ' + str(radius) + ' ' + str(t0) + ' ' + str(t1), shell=True)


print "read from csv file",filename
data_= np.genfromtxt(filename, delimiter=',', dtype='float')





vtu_filenames = []
for nr,t in enumerate(timesteps):
  t_min = t - dt
  t_max = t + dt

  i_min = max(0, int(t_min - t0))
  i_max = int(t_max - t0)

  data = data_[:,i_min:i_max]

  M = data.shape[0]
  N = data.shape[1]

  print "calc mean flow"
  mean = [0]*M
  for i in range(0, M):
    mean[i] = np.mean(data[i,:])

  shift = M/3
  print "write mean velocity to file"

  vtu_filename = "velocity" + postfix + "_" + str(nr) + ".vtu"
  with open(out_dir + "/" + vtu_filename, "w") as outfile:
    with open("sphere/ingo/sphereIko2k_head.vtu", "r") as head:
      for line in head:
        outfile.write(line)

    outfile.write('<DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">\n')
    for i in range(0, shift):
      outfile.write(str(mean[i]) + " " + str(mean[i+shift]) + " " + str(mean[i+2*shift]) + "\n")
    outfile.write('</DataArray>')

    with open("sphere/ingo/sphereIko2k_tail.vtu", "r") as tail:
      for line in tail:
        outfile.write(line)

  vtu_filenames.append(vtu_filename)


pvd_filename = out_dir + "/velocity" + postfix + ".pvd"
with open(pvd_filename, "w") as outfile:
  outfile.write('<?xml version="1.0"?>\n')
  outfile.write('<VTKFile type="Collection" version="0.1">\n')
  outfile.write('<Collection>\n')

  for t,f in zip(timesteps, vtu_filenames):
    outfile.write('<DataSet timestep="' + str(t) + '" part="0" file="' + f + '"/>\n')

  outfile.write('</Collection>\n')
  outfile.write('</VTKFile>')






