#!/usr/bin/python
import sys
import os
import glob
import fnmatch
import subprocess
from timeit import default_timer as timer

basedir = '/scratch/wir/praetori/projects/polar_pfc/results'
setup = 'polar_pfc'

out_dir = '/media/Home/projects/shtns/results/polar_pfc'

use_rsync = False

if use_rsync:
  command = "rsync --prune-empty-dirs --include='*.csv' --include='*/' --exclude='*' -r --info=progress2  -e ssh "
else:
  command = "scp "

dt = 1.0
for run in range(0,5):
  for radius in [10,20,40,60,80,100]:
    for V0 in [0.1, 0.2, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.425, 0.45, 0.475, 0.5,0.6,0.7,0.8,0.9,1.0]:
      dirname = basedir + '/' + setup + '/radius_' + str(radius) + '_v0_' + str(V0) + '_run_' + str(run)
      print 'copy',dirname
      if use_rsync:
        subprocess.call(command + "praetori@taurusexport.hrsk.tu-dresden.de:" + dirname + " " + out_dir, shell=True)
      else:
        if not os.path.exists(out_dir):
          os.makedirs(out_dir)
        start = timer()
        subprocess.call(command + "praetori@taurusexport.hrsk.tu-dresden.de:" + dirname + "/data.tar.gz " + out_dir, shell=True)
        subprocess.call("tar --strip-components=7 -xzf " + out_dir + "/data.tar.gz -C " + out_dir, shell=True)
        end = timer()
        print "Time:",(end-start),"sec"
