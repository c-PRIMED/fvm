import sys
import os

import vtk
vtkxdmf = 0
try:
  import libvtkXdmfPython
  vtkxdmf = libvtkXdmfPython
except:
  import vtkXdmfPython
  vtkxdmf = vtkXdmfPython
import sys

if len(sys.argv) != 4:
  print "Usage: %s <data directory> <input dataset> <output file name>"
  sys.exit(1)

reader = vtk.vtkDataSetReader()
reader.SetFileName(sys.argv[2])

bp = vtk.vtkBrownianPoints()
bp.SetInput(reader.GetOutput())
bp.Update()

image = bp.GetOutput()
array = image.GetPointData().GetArray("BrownianVectors")
image.GetPointData().RemoveArray("BrownianVectors")
array.SetName("UnnamedNodeArray0")
image.GetPointData().AddArray(array)
print image

dataext = image.GetWholeExtent()

globalFactors = {
  0: ( 00, 35, 00, 35, 00, 35,  ),
  1: ( 32, 67, 00, 35, 00, 35,  ),
  2: ( 00, 35, 32, 67, 00, 35,  ),
  3: ( 32, 67, 32, 67, 00, 35,  ),
  4: ( 00, 35, 00, 35, 32, 67,  ),
  5: ( 32, 67, 00, 35, 32, 67,  ),
  6: ( 00, 35, 32, 67, 32, 67,  ),
  7: ( 32, 67, 32, 67, 32, 67,  ),
}

# Create Writers
writer_id = vtkxdmf.vtkXdmfWriter()
writer_id.SetDomainName("IronProtein")
writer_id.SetFileName("%s.ImageData.xmf" % os.path.join(sys.argv[1], sys.argv[3]))
writer_id.SetGridName("Iron0")

writer_rg = vtkxdmf.vtkXdmfWriter()
writer_rg.SetDomainName("IronRectProtein")
writer_rg.SetFileName("%s.RectilinearGrid.xmf" % os.path.join(sys.argv[1], sys.argv[3]))
writer_rg.SetGridName("Iron0")

writer_sg = vtkxdmf.vtkXdmfWriter()
writer_sg.SetDomainName("IronRectProtein")
writer_sg.SetFileName("%s.StructuredGrid.xmf" % os.path.join(sys.argv[1], sys.argv[3]))
writer_sg.SetGridName("Iron0")

origin = image.GetOrigin()
spacing = image.GetSpacing()

# Create rectilinear grid
recgrid = vtk.vtkRectilinearGrid()
recgrid.SetWholeExtent(image.GetWholeExtent())
recgrid.SetExtent(image.GetExtent())
recgrid.SetUpdateExtent(image.GetUpdateExtent())
recgrid.GetPointData().ShallowCopy(image.GetPointData())
recgrid.GetCellData().ShallowCopy(image.GetCellData())

darray = vtk.vtkDoubleArray()
for a in range(dataext[0], dataext[1]+1):
  darray.InsertNextValue( origin[0] + (a-dataext[0]) * spacing[0] )
recgrid.SetXCoordinates(darray)   
darray = vtk.vtkDoubleArray()
for a in range(dataext[2], dataext[3]+1):
  darray.InsertNextValue( origin[1] + (a-dataext[2]) * spacing[1] )
recgrid.SetYCoordinates(darray)   
darray = vtk.vtkDoubleArray()
for a in range(dataext[4], dataext[5]+1):
  darray.InsertNextValue( origin[2] + (a-dataext[4]) * spacing[2] )
recgrid.SetZCoordinates(darray)   

# Create structured grid
strgrid = vtk.vtkStructuredGrid()
strgrid.SetWholeExtent(image.GetWholeExtent())
strgrid.SetExtent(image.GetExtent())
strgrid.SetUpdateExtent(image.GetUpdateExtent())
strgrid.GetPointData().ShallowCopy(image.GetPointData())
strgrid.GetCellData().ShallowCopy(image.GetCellData())

pts = vtk.vtkPoints();
for idx in range(image.GetNumberOfPoints()):
  pts.InsertNextPoint(image.GetPoint(idx))
strgrid.SetPoints(pts)

# Split to parts
for grid in range(8):
  print "Create data: %s" % grid
  factors = globalFactors[grid]
  ext = []
  for ex in range(3):
    ext.append(factors[2*ex])
    ext.append(factors[2*ex+1])
  print dataext
  print ext

  clip_id = vtk.vtkImageClip()
  clip_id.SetInput(image)
  clip_rg = vtk.vtkRectilinearGridClip()
  clip_rg.SetInput(recgrid)
  clip_sg = vtk.vtkStructuredGridClip()
  clip_sg.SetInput(strgrid)

  clip_id.SetOutputWholeExtent(ext[0], ext[1], ext[2], ext[3], ext[4], ext[5])
  clip_id.Update()
  clip_rg.SetOutputWholeExtent(ext[0], ext[1], ext[2], ext[3], ext[4], ext[5])
  clip_rg.Update()
  clip_sg.SetOutputWholeExtent(ext[0], ext[1], ext[2], ext[3], ext[4], ext[5])
  clip_sg.Update()

  writer_id.AddInput(clip_id.GetOutput())
  writer_rg.AddInput(clip_rg.GetOutput())
  writer_sg.AddInput(clip_sg.GetOutput())

# Write files
writer_id.Write()
writer_rg.Write()
writer_sg.Write()


