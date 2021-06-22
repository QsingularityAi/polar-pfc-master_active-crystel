#!/usr/bin/python
import sys
import os
import glob
import fnmatch
import csv
import numpy as np
import cyflann
import matplotlib.pyplot as plt
from matplotlib2tikz import save as tikz_save

if len(sys.argv) < 3:
  print "usage: ./generate_mesh_director 'filenames_*.csv' postfix"
  exit()

print "read mesh"
areas = []
coords = []
#with open('sphere/R_100_theta_64.csv', 'r') as csvfile:
with open('sphere/ingo/sphereIko2k.csv', 'r') as csvfile:
  csvreader = csv.DictReader(csvfile, delimiter=',')

  for row in csvreader:
    coord = (float(row['Points:0']), float(row['Points:1']), float(row['Points:2']))
    area = float(row['Quality'])

    coords.append(np.array(coord))
    areas.append(area)

time_interval = [6000, 8000]
n_neighbours = 7

base_filename = sys.argv[1]
dt_write = 0.5
d = 4*np.pi/np.sqrt(3.0)

base_dir = os.path.dirname(base_filename)
base_filename = os.path.basename(base_filename)

directory = os.listdir(base_dir)
filenames = fnmatch.filter(directory, base_filename)

print len(filenames)

print "read files"
mesh_velocities = []
timesteps = []
for f_ in filenames:
  f = base_dir + '/' + f_;

  nr = f_.split('.')[0].split('_')[1]
  time = float(nr) * dt_write

  if float(nr) < time_interval[0] or float(nr) > time_interval[1]:
    continue

  timesteps.append(time)

  positions = []
  velocities = []
  print "read file",f
  with open(f, 'r') as csvfile:
    csvreader = csv.DictReader(csvfile, delimiter=',')

    for row in csvreader:
      pos = (float(row['x']), float(row['y']), float(row['z']))
      vel = (float(row['d0']), float(row['d1']), float(row['d2']))

      positions.append(np.array(pos))
      velocities.append(np.array(vel))

  kdtree = cyflann.FLANNIndex(algorithm='kdtree_single')
  kdtree.build_index(positions)

  indices, distances = kdtree.nn_index(coords, n_neighbours)

  weights = [0.0]*len(coords)
  mesh_vel = [np.array((0.,0.,0.))]*len(coords)

  for i in range(0, len(coords)):
    for j in range(0, n_neighbours):
      idx = indices[i][j]
      dist = distances[i][j]

      if dist < 2.5*d:
        w = 1.0 / max(1.e-5, dist)
        mesh_vel[i] = mesh_vel[i] + (velocities[idx] * w)
        weights[i] += w

  for i in range(0, len(mesh_vel)):
    if weights[i] > 1.e-5:
      mesh_vel[i] = mesh_vel[i] * (1.0/weights[i])

  mesh_velocities.append(mesh_vel)

print "sort by timestep"
data = zip(timesteps, mesh_velocities)
sorted_data = sorted(data, cmp=lambda x,y: cmp(x[0], y[0]))
timesteps, mesh_velocities = zip(*sorted_data)


if len(sys.argv) > 2:
  postfix = '_' + sys.argv[2]
else:
  postfix =''

print "write data to csv file"
row = [0]*(len(mesh_velocities) + 1)
with open('pca_dir_' + str(time_interval[0]) + '-' + str(time_interval[1]) + postfix + '.csv','w') as out_file:
  datawriter = csv.writer(out_file, delimiter=',')

  for r_ in range(0,3*len(areas)):
    r = r_ % len(areas)
    i = r_ / len(areas)

    row[0] = areas[r]
    for j in range(0,len(mesh_velocities)):
      row[j+1] = mesh_velocities[j][r][i]

    datawriter.writerow(row)

