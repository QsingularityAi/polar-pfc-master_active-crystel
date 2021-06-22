#!/usr/bin/python
import subprocess

def cleanup(string):
  return string.replace('.', '_')
# {end cleanup}

setup0 = 'dislocations'

geometries = zip([10, 20, 40, 80], [0,1,2,3], ['02','04','08','16'])

r = -0.98
for run in range(0,2):
  for radius,res,h in geometries:
    setup = cleanup(setup0 + '_R' + str(radius))
    for V0 in [0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35]:
      postfix = cleanup('v_' + str(V0) + '_run_' + str(run))
      print 'run: ' + setup + '/' + postfix

      subprocess.call('./paraview_script.py ' + str(r) + ' ' + str(run) + ' ' + str(radius) + ' ' + str(V0), shell=True)

