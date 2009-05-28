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
Reader.SetFileName('Points1.xmf')
# Reader.DebugOn()
Reader.UpdateInformation()
Reader.DisableAllGrids()
Reader.EnableGrid(2)
Reader.EnableAllArrays()
Reader.Update()
print 'Output = ', Reader.GetOutput()
Append = vtkAppendFilter()
print ProcId," : Ports = ", Reader.GetNumberOfOutputPorts()
for on in range(Reader.GetNumberOfOutputPorts()) :
        Output = Reader.GetOutput(on)
        print ProcId, '  Number of Levels ',Output.GetNumberOfLevels()
        print ProcId, '  Number of Groups ',Output.GetNumberOfGroups()
        print ProcId, '  Number of DataSets in Group 0 ',Output.GetNumberOfDataSets(0)
        for i in range(Output.GetNumberOfDataSets(0)) :
            ds = Output.GetDataSet(0,i)
            if ds :
                print '%d : Output %d = %d Cells' % (ProcId, i, ds.GetNumberOfCells())
                Append.AddInput(ds)
            else :
                print '%d : Output %d = NONE' % (ProcId, i)
Append.Update()
print Append.GetOutput()

