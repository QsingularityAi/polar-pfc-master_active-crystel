#!/usr/bin/python
import sys
import os
import glob
import fnmatch
import csv
import numpy as np

if len(sys.argv) > 2:
  print "usage: ./generate_time_series [postfix]"
  exit()

if len(sys.argv) > 1:
  postfix = '_' + sys.argv[1]
else:
  postfix =''

filename_pattern = 'particles_*.csv'
base_dir_ = '/media/Home/projects/shtns/results/polar_pfc_new'

dt = 0.5
for radius in [10,20,40,60,80,100]:
  for V0 in [0.1, 0.2, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.425, 0.45, 0.475, 0.5,0.6,0.7,0.8,0.9,1.0]:

    base_dir = base_dir_ + '/radius_' + str(radius) + '_v0_' + str(V0) + '_dt_' + str(dt)

    timesteps = []
    director = []

    directory = os.listdir(base_dir)
    filenames = fnmatch.filter(directory, filename_pattern)

    n_files = len(filenames)

    if n_files < 500:
      continue

    print 'radius[',radius,'], v0[',V0,'], dt[',dt,']: ',n_files,'files'

    nf = 0
    for f_ in filenames:
      f = base_dir + '/' + f_;

      nr = f_.split('.')[0].split('_')[1]
      time = float(nr) * dt

      with open(f, 'r') as csvfile:
        csvreader = csv.DictReader(csvfile, delimiter=',')

        num_values = 0
        director_value = np.zeros(3)
        for row in csvreader:
          director_value = director_value + np.array([float(row['d0']), float(row['d1']), float(row['d2'])])
          num_values += 1
        # {end for csvreader}

        director_value = director_value * (1.0/float(num_values))

        if nf > 0 and nf % (n_files/10) == 0:
          print round(nf*100.0/n_files,2),'%: ',num_values

        timesteps.append(time)
        director.append(director_value)
      # {end with}

      nf += 1
    # {end for filenames}

    data = zip(timesteps, director)
    sorted_data = sorted(data, cmp=lambda x,y: cmp(x[0], y[0]))
    timesteps, director = zip(*sorted_data)

    with open(base_dir + '/director' + postfix + '.csv','w') as out_file:
      fieldnames = ['t','d0','d1','d2']
      writer = csv.DictWriter(out_file, fieldnames=fieldnames)
      writer.writeheader()
      for data in zip(timesteps, director):
        writer.writerow({'t': data[0], 'd0': data[1][0], 'd1': data[1][1], 'd2': data[1][2]})
    # {end with}
  # {end for V0}
# {end for radius}
