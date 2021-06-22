#!/usr/bin/python
from paraview.simple import *
import sys
import glob
import os.path

base_filename = ''
if len(sys.argv) > 1:
  base_filename = sys.argv[1]
else:
  sys.stderr.write('Error: no filenames given as first argument!\n')
  
max_files = 500

timesteps = []
if len(base_filename) > 0:
  filenames = glob.glob(base_filename)
  n_files = len(filenames)
  print n_files
  
  i = 0
  j = 0
  for f in filenames:
    i = i+1
    if '.vtu_delete' in f:
      continue
    elif '_bin.vtu' in f:
      if os.path.getsize(f) > 1000:
	idx = f.rfind('_', 0, -8)
	time = float(f[(idx+1):-8]);
	timesteps.append([time, f])
      else:
	os.rename(f, f + '_delete')
	# end(if)
      continue
      # end(if)
    
    reader =  XMLUnstructuredGridReader(FileName=[f])
    
    new_f = f.replace('.vtu', '_bin.vtu')
    #reader = OpenDataFile(f)
    w = CreateWriter(new_f, reader)
    w.DataMode = 'Appended'
    w.CompressorType = 'ZLib'
    w.UpdatePipeline()
    Delete(reader)
    del reader
    del w
    
    if os.path.isfile(new_f):
      os.remove(f)
      # end(if)
      
    if os.path.getsize(new_f) > 1000:
      idx = f.rfind('_', 0, -4)
      time = float(f[(idx+1):-4]);
      timesteps.append([time, new_f])
      j = j + 1
      print repr(i) + ': ' + f
      #end(if)
    
    if j > max_files:
      print 'maximum nr. of files reached. Please restart the script to proceed!'
      break
      # end(if)
      
    #end(for)
  
  #end(if)

if len(timesteps) > 0 and len(sys.argv) > 2:
  pvd_filename = sys.argv[2]
  new_pvd = pvd_filename.replace('.pvd', '_bin.pvd')
  
  fid = open(new_pvd, 'w')
  fid.write('<?xml version="1.0"?>\n')
  fid.write('<VTKFile type="Collection" version="0.1" ><Collection>\n')
  
  timesteps.sort()
  for t in timesteps:
    fid.write('<DataSet timestep="' + repr(t[0]) + '" part="0" file="' + t[1] + '"/>\n')
    # end(for)
  
  fid.write('</Collection></VTKFile>')
  fid.close()
  # end(if)
  
# exec(open("/home/mai/sources/paraview_scripts/tapete_slices.py").read())
 
