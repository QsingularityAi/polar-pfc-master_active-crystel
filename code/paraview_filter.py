input = self.GetInputDataObject(0, 0)
output = self.GetOutputDataObject(0)
num_points = input.GetNumberOfPoints()
num_cells = input.GetNumberOfCells()

quality = vtk.vtkDoubleArray()
quality.SetName("quality")
quality.SetNumberOfValues(num_points)

for i in range(0,num_points):
  quality.SetValue(i,0)

# traverse all cells in the (surface) mesh
for cid in range(0,num_cells):
  cell = input.GetCell(cid)
  num_edges = cell.GetNumberOfEdges()
  for eid in range(num_edges):
    edge = cell.GetEdge(eid)
    points = edge.GetPointIds()
    p0 = points.GetId(0)
    p1 = points.GetId(1)

    quality.SetValue(p0, quality.GetValue(p0)+1)
    quality.SetValue(p1, quality.GetValue(p1)+1)

# all edges are visited twice
for i in range(0,num_points):
  quality.SetValue(i,quality.GetValue(i)/2)

num_defects = [0]*12
for i in range(0,num_points):
  n = max(0, min(11, int(quality.GetValue(i))))
  num_defects[n] += 1

directory = '/media/Daten/projects/shtns/results/v_0_29_run_0'
with open(directory + '/num_defects.csv', 'a') as out:
  out.write(' '.join(map(str, num_defects)) + '\n')

output.GetPointData().AddArray(quality)
