#!/usr/bin/python
import sys
import os
import glob
import fnmatch
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib2tikz import save as tikz_save

if len(sys.argv) < 4:
  print "usage: ./plot_order_parameter 'filenames_*.csv' timestep postfix"
  print "or..."
  print "usage: ./plot_order_parameter (this reads the file 'time_series.csv' and prints it directly)"
  exit()

if len(sys.argv) > 3:
  postfix = '_' + sys.argv[3]
else:
  postfix =''


if len(sys.argv) <= 1:
  timesteps= []
  polar = []
  nematic = []
  velocity = []

  with open('time_series.csv','r') as in_file:
    reader = csv.DictReader(in_file, delimiter=',')

    tol=1.0e-2

    for row in reader:
      t = float(row['t'])
      p = float(row['polar'])
      n = float(row['nematic'])
      v = float(row['velocity'])

      if len(timesteps) == 0:
        timesteps.append(t)
        polar.append(p)
        nematic.append(n)
        velocity.append(v)
      else:
        t_last = timesteps[len(timesteps)-1]
        p_last = polar[len(polar)-1]
        n_last = nematic[len(nematic)-1]

        if np.sqrt(((t_last - t)/5000)**2 + (p_last - p)**2 + (n_last - n)**2) > tol:
          timesteps.append(t)
          polar.append(p)
          nematic.append(n)
          velocity.append(v)


else:
  base_filename = sys.argv[1]

  dt_write = 0.1
  if len(sys.argv) > 2:
    dt_write = float(sys.argv[2])

  base_dir = os.path.dirname(base_filename)
  base_filename = os.path.basename(base_filename)

  max_files = 500

  timesteps = []
  polar = []
  nematic = []
  velocity = []

  directory = os.listdir(base_dir)
  filenames = fnmatch.filter(directory, base_filename)

  n_files = len(filenames)
  print n_files

  nf = 0
  for f_ in filenames:
    f = base_dir + '/' + f_;

    nr = f_.split('.')[0].split('_')[1]
    time = float(nr) * dt_write

    with open(f, 'r') as csvfile:
      csvreader = csv.DictReader(csvfile, delimiter=',')

      num_values = 0
      polar_value = 0.0
      nematic_value = 0.0
      velocity_value = 0.0
      for row in csvreader:
        polar_value += float(row['polar'])
        nematic_value += float(row['nematic'])
        velocity_value += np.sqrt(float(row['v0'])**2 + float(row['v1'])**2 + float(row['v2'])**2)
        num_values += 1

      polar_value /= float(num_values)
      nematic_value /= float(num_values)
      velocity_value /= float(num_values)
      print round(nf*100.0/n_files,2),'%: ',num_values

      timesteps.append(time)
      polar.append(polar_value)
      nematic.append(nematic_value)
      velocity.append(velocity_value)

    nf += 1

  print len(timesteps), len(polar), len(nematic)

  data = zip(timesteps, polar, nematic, velocity)
  sorted_data = sorted(data, cmp=lambda x,y: cmp(x[0], y[0]))
  timesteps, polar, nematic, velocity = zip(*sorted_data)

  with open('time_series' + postfix + '.csv','w') as out_file:
    fieldnames = ['t','polar','nematic','velocity']
    writer = csv.DictWriter(out_file, fieldnames=fieldnames)
    writer.writeheader()
    for data in zip(timesteps, polar, nematic, velocity):
      writer.writerow({'t': data[0], 'polar': data[1], 'nematic': data[2], 'velocity': data[3]})

# plot timeseries

plt.plot(timesteps, polar, timesteps, nematic)
plt.xlabel('time t')
plt.legend(['polar order $\\pi(t)$', 'nematic order $\\nu(t)$'], loc='lower right')

tikz_save('order_parameter' + postfix + '.tikz',
           figureheight = '\\figureheight',
           figurewidth = '\\figurewidth')

plt.show()

plt.plot(timesteps, velocity)
plt.xlabel('time t')
plt.legend(['mean particle velocity $\\mathbf{v}(t)$'], loc='upper right')

tikz_save('velocity' + postfix + '.tikz',
           figureheight = '\\figureheight',
           figurewidth = '\\figurewidth')
plt.show()
