#!/bin/env python

from Xdmf import *
from vtk import *
from libvtkXdmfPython import *


Reader = vtkXdmfReader()
Controller = vtkMPIController()
Reader.SetController(Controller)
ProcId = Reader.GetController().GetLocalProcessId()
NumProcs = Reader.GetController().GetNumberOfProcesses()
print 'Hello from %d of %d' % (ProcId, NumProcs)
Reader.SetFileName('Mixed.xmf')
# Reader.DebugOn()
Reader.UpdateInformation()
Reader.EnableAllGrids()
Reader.EnableAllArrays()
Reader.Update()
print 'Output = ', Reader.GetOutput()
Append = vtkAppendFilter()
print ProcId," : Ports = ", Reader.GetNumberOfOutputPorts()
RenderWindow = vtkRenderWindow()
Renderer = vtkRenderer()
RenderWindow.AddRenderer(Renderer)
RenderWindowInteractor = vtkXdmfRenderWindowInteractor()
RenderWindowInteractor.SetLightFollowCamera(0)
RenderWindowInteractor.GetInteractorStyle().SetCurrentStyleToTrackballCamera()
RenderWindowInteractor.GetInteractorStyle().SetAutoAdjustCameraClippingRange(1)
RenderWindowInteractor.SetRenderWindow(RenderWindow)
RenderWindowInteractor.Initialize()
for on in range(Reader.GetNumberOfOutputPorts()) :
        ds = Output = Reader.GetOutput(on)
        print 'Output has %d Cells', Output.GetNumberOfCells()
#        for i in range(Output.GetNumberOfCells()) :
#            cell = Output.GetCell(i)
#            print 'Cell = ', cell
        Geometry = vtkGeometryFilter()
        Geometry.SetInput(Output)
        Mapper = vtkPolyDataMapper()
        Mapper.GetLookupTable().SetHueRange( .667, 0.0 )
        Mapper.SetInput(Geometry.GetOutput())
        Mapper.SetScalarRange(Output.GetScalarRange())
        Actor = vtkActor()
        Actor.SetMapper(Mapper)
        Renderer.AddActor(Actor)
        i = 0
#        if ds :
#            print '%d : Output %d = %d Cells' % (ProcId, i, ds.GetNumberOfCells())
#            Append.AddInput(ds)
#        else :
#            print '%d : Output %d = NONE' % (ProcId, i)
#Append.Update()
#print Append.GetOutput()
RenderWindowInteractor.Start(1)

