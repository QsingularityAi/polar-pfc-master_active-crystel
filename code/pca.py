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

if len(sys.argv) < 2:
  print "usage: ./pca data.csv [postfix]"
  exit()

if len(sys.argv) > 2:
  postfix = '_' + sys.argv[2]
else:
  postfix =''


filename = sys.argv[1]

print "read from csv file"
data_= np.genfromtxt(filename, delimiter=',', dtype='float')
data = data_[:,1:]
areas = data_[:,0]

M = data.shape[0]
N = data.shape[1]

x = data.copy()

print M/3,N

print "calc mean and area"
mean = [0]*M
area = 0
for i in range(0, M):
  mean[i] = np.mean(data[i,:])
  area += areas[i]

print "shift data by time-mean"
for j in range(0, N):
  data[:,j] = data[:,j] - mean

print "rescale data by area of grid-cells"
for i in range(0, M):
  data[i,:] = data[i,:] * np.sqrt( areas[i]/area )

print "calculate covariance matrix"
cov_matrix = np.zeros((N, N))
for i in range(0, N):
  cov_matrix[i,i] = data[:,i].dot(data[:,i])
  for j in range(i+1, N):
    cov_matrix[i,j] = data[:,i].dot(data[:,j])
    cov_matrix[j,i] = cov_matrix[i,j]


print "compute eigenvalues and eigenvectors of covariance matrix"
eig_val, eig_vec = np.linalg.eig(cov_matrix)

print "Make a list of (eigenvalue, eigenvector) tuples"
eig_pairs = [(np.abs(eig_val[i]), eig_vec[:,i]) for i in range(len(eig_val))]

print "Sort the (eigenvalue, eigenvector) tuples from high to low"
eig_pairs.sort(key=lambda x: x[0], reverse=True)

eig_val, eig_vec = zip(*eig_pairs)

y0 = [0]*M
y1 = [0]*M
y2 = [0]*M
for j in range(0, N):
  for i in range(0, M):
    y0[i] += eig_vec[0][j][i] * data[i,j] * np.sqrt( areas[i]/area )
    y1[i] += eig_vec[1][j][i] * data[i,j] * np.sqrt( areas[i]/area )
    y2[i] += eig_vec[2][j][i] * data[i,j] * np.sqrt( areas[i]/area )

shift = M/3
print "write mean velocity to file"
with open("sphere/velocity_R_100_2k" + postfix + ".vtu", "w") as outfile:
  with open("sphere/ingo/sphereIko2k_head.vtu", "r") as head:
    for line in head:
      outfile.write(line)

  outfile.write('<DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">\n')
  for i in range(0, shift):
    outfile.write(str(mean[i]) + " " + str(mean[i+shift]) + " " + str(mean[i+2*shift]) + "\n")
  outfile.write('</DataArray>')

  outfile.write('<DataArray type="Float32" Name="Eig0" NumberOfComponents="3" format="ascii">\n')
  for i in range(0, shift):
    outfile.write(str(y0[i]) + " " + str(y0[i+shift]) + " " + str(y0[i+2*shift]) + "\n")
  outfile.write('</DataArray>')

  outfile.write('<DataArray type="Float32" Name="Eig1" NumberOfComponents="3" format="ascii">\n')
  for i in range(0, shift):
    outfile.write(str(y1[i]) + " " + str(y1[i+shift]) + " " + str(y1[i+2*shift]) + "\n")
  outfile.write('</DataArray>')

  outfile.write('<DataArray type="Float32" Name="Eig2" NumberOfComponents="3" format="ascii">\n')
  for i in range(0, shift):
    outfile.write(str(y2[i]) + " " + str(y2[i+shift]) + " " + str(y2[i+2*shift]) + "\n")
  outfile.write('</DataArray>')

  with open("sphere/ingo/sphereIko2k_tail.vtu", "r") as tail:
    for line in tail:
      outfile.write(line)



plt.semilogy(eig_val,'.')

tikz_save('pca_analysis' + postfix + '.tikz',
           figureheight = '\\figureheight',
           figurewidth = '\\figurewidth')
plt.show()



P0 = [0]*N
P1 = [0]*N

for j in range(0, N):
  P0[j] = 0;
  P1[j] = 0;
  for i in range(0, M):
    P0[j] += y0[i] * data[i,j] * np.sqrt( areas[i]/area )
    P1[j] += y1[i] * data[i,j] * np.sqrt( areas[i]/area )

plt.plot(P0, '-b')
plt.plot(P1, '-g')
plt.show()
