#!/opt/software/paraview/5.1.0/bin/pvpython
from contextlib import contextmanager
import sys, os
import fnmatch
import datetime
from timeit import default_timer as timer

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

def generate_filter(in_filename, params):
  out_ = ''
  with open(in_filename,'r') as in_:
    for line in in_:
      for p,v in params.items():
        line = line.replace(p, v)
      out_ += line + '\n'
  return out_
# {end def}

def cleanup(string):
    return string.replace('.', '_')
# {end cleanup}

def reset_session():
    pxm = servermanager.ProxyManager()
    pxm.UnRegisterProxies()
    del pxm
    Disconnect()
    Connect()
# {end reset_session}

@contextmanager
def suppress_stdout():
  with open(os.devnull, "w") as devnull:
    old_stdout = sys.stdout
    sys.stdout = devnull
    try:
      yield
    finally:
      sys.stdout = old_stdout
# {end suppress_stdout}

def write_timestep(Source, out):
  with suppress_stdout():
    Input = servermanager.Fetch(Source)
  num_points = Input.GetNumberOfPoints()
  num_cells = Input.GetNumberOfCells()

  quality = [0]*num_points

  # traverse all cells in the (surface) mesh
  for cid in range(0,num_cells):
    cell = Input.GetCell(cid)
    num_edges = cell.GetNumberOfEdges()
    for eid in range(num_edges):
      edge = cell.GetEdge(eid)
      points = edge.GetPointIds()
      p0 = points.GetId(0)
      p1 = points.GetId(1)

      quality[p0] += 1
      quality[p1] += 1

  # all edges are visited twice
  for i in range(0,num_points):
    quality[i] /= 2.0

  num_defects = [0]*12
  for i in range(0,num_points):
    n = max(0, min(11, int(quality[i])))
    num_defects[n] += 1

  out.write(' '.join(map(str, num_defects)) + '\n')
# {end write_timestep}

def progress(count, total, suffix=''):
  bar_len = 60
  filled_len = int(round(bar_len * count / float(total)))

  percents = round(100.0 * count / float(total), 1)
  bar = '#' * filled_len + '-' * (bar_len - filled_len)

  sys.stdout.write('\r[%s] %s%s %s' % (bar, percents, '%', suffix))
  if count == total:
    sys.stdout.write('\n')
  sys.stdout.flush()
# {end progress}

# ------------------------------------------------------------------------------

basedir = '/media/Daten/projects/polar_pfc/results'
setup0 = 'dislocations'

r = float(sys.argv[1])
run = int(float(sys.argv[2]))
radius = int(float(sys.argv[3]))
V0 = float(sys.argv[4])

setup = cleanup(setup0 + '_R' + str(radius))
postfix = cleanup('v_' + str(V0) + '_run_' + str(run))
dirname = basedir + '/' + setup + '/' + postfix
directory = os.listdir(dirname)
csvfiles = fnmatch.filter(directory, 'p*.csv')

if os.path.exists(dirname + '/num_defects.csv'):
  os.remove(dirname + '/num_defects.csv')

csvfiles.sort()


start = timer()
with open(dirname + '/num_defects.csv', 'w') as out:
  for i,f in enumerate(csvfiles):
    reset_session()
    particles = CSVReader(FileName=dirname + '/' + f)

    # create a new 'Table To Points'
    tableToPoints1 = TableToPoints(Input=particles, XColumn='x', YColumn='y', ZColumn='z')

    # create a new 'Delaunay 3D'
    delaunay3D1 = Delaunay3D(Input=tableToPoints1)
    extractSurface1 = ExtractSurface(Input=delaunay3D1)

    write_timestep(extractSurface1, out)
    progress(i+1, len(csvfiles),'File: ' + f)

end = timer()
print "Duration:",str(datetime.timedelta(seconds=(end-start)))
