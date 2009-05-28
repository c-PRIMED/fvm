#!/bin/env python

from Xdmf import *
from vtk import *
from libvtkXdmfPython import *


class DataParser :
    CellTypes = {
                1 : 'Vertex',
                2 : 'PolyVertex',
                3 : 'Line',
                4 : 'PolyLine',
                5 : 'Triangle',
                6 : 'TriangleStrip',
                7 : 'Polygon',
                8 : 'Pixel',
                9 : 'Quadrilateral',
                10 : 'Tetrahedron',
                11 : 'Voxel',
                12 : 'Hexahedron',
                13 : 'Wedge',
                14 : 'Pyramid',
                99 : 'Unknown'
                }
    def __init__(self) :
        self.DOM = XdmfDOM()
        self.Root = XdmfRoot()
        self.Root.SetDOM(self.DOM)
        self.Root.SetVersion(2.0)
        self.Root.Build()
        self.Root.DebugOn()
        self.Domain = XdmfDomain()
        self.Root.Insert(self.Domain)
        self.MainGrid = XdmfGrid()
        self.MainGrid.SetName('Main Grid')
        self.MainGrid.SetGridType(XDMF_GRID_TREE)
        self.Domain.Insert(self.MainGrid)
        self.Grids = []
        self.Arrays = []

    def ParseUGrid(self, UGrid) :
        if UGrid.IsHomogeneous() :
            print 'Homogeneous UGrid ', UGrid
            Fd = UGrid.GetFieldData()
            if Fd :
                NameArray = Fd.GetArray('Name')
                print 'NameArray = ', NameArray
                print 'name = ', NameArray.GetValue(0)
            CellArray = UGrid.GetCells()
            Conns = CellArray.GetData()
            print '%d Cells' % UGrid.GetNumberOfCells()
            for i in range(UGrid.GetNumberOfCells()) :
                print '(%d) %d = %s' % (i, UGrid.GetCellType(i), self.CellTypes[UGrid.GetCellType(i)])
            print '%d Points' % UGrid.GetNumberOfPoints()
            for i in range(UGrid.GetNumberOfPoints()) :
                x, y, z = UGrid.GetPoint(i)
                print '(%d) %f %f %f' % (i, x, y, z)
            print '%d Connections' % CellArray.GetNumberOfConnectivityEntries()
            i = 0
            while i < Conns.GetNumberOfTuples() :
                n = Conns.GetValue(i)
                i += 1
                for j in range(n) :
                    p = Conns.GetValue(i)
                    i += 1
                    print '%d ' % p
            Pd = UGrid.GetPointData()
            print '# of Point Attributes = %d' %  Pd.GetNumberOfArrays()
            for i in range(Pd.GetNumberOfArrays()) :
                Pa = Pd.GetArray(i)
                print Pd.GetArrayName(i), ' = ',Pa
                for j in range(Pa.GetNumberOfTuples()) :
                    print '(%d) %f' % (j, Pa.GetValue(j))
            Cd = UGrid.GetCellData()
            print '# of Cell Attributes = %d' %  Cd.GetNumberOfArrays()
            # print UGrid
        else :
            print 'Heterogeneous UGrid'

    def DataArrayToXdmfArray(self, da) :
        # print 'Converting ', da
        print 'Data Type', da.GetDataType(), ' = ', da.GetDataTypeAsString()
        Xda = XdmfArray()
        Xda.SetNumberOfElements(da.GetNumberOfTuples() * da.GetNumberOfComponents())
        Type = da.GetDataTypeAsString()
        pntr = da.GetVoidPointer(0)
        Xpntr = VoidPointerHandleToXdmfPointer(pntr)
        if Type.upper() == 'FLOAT' :
            # print 'Array is Float32'
            Xda.SetNumberType(XDMF_FLOAT32_TYPE)
        elif Type.upper() == 'INT' :
            # print 'Array is Int32'
            Xda.SetNumberType(XDMF_INT32_TYPE)
        elif Type.upper() == 'DOUBLE' :
            # print 'Array is Float64'
            Xda.SetNumberType(XDMF_FLOAT64_TYPE)
        elif Type.upper() == 'LONG' :
            # print 'Array is Int64'
            Xda.SetNumberType(XDMF_INT64_TYPE)
        elif Type.upper() == 'IDTYPE' :
            # print 'Array is Int64'
            Xda.SetNumberType(XDMF_INT64_TYPE)
        else :
            print 'Illegal NumberType : ', Type
            return None
        Xda.SetDataPointer(Xpntr)
        print 'Values ',Xda.GetValues(0, 10)
        return Xda

    def WriteIGrid(self, IGrid, Group, Index) :
        print 'Homogeneous Image UGrid ', IGrid.GetDimensions()

    def WriteUGrid(self, UGrid, Group, Index) :
        if UGrid.IsHomogeneous() :
            print 'Homogeneous UGrid '
            Fd = UGrid.GetFieldData()
            if Fd :
                NameArray = Fd.GetArray('Name')
                print 'NameArray = ', NameArray
                print 'name = ', NameArray.GetValue(0)
            Xgrid = XdmfGrid()
            Xgrid.SetName('Group %05d Index %05d' % (Group, Index))
            self.MainGrid.Insert(Xgrid)
            self.Grids.append(Xgrid)
            CellArray = UGrid.GetCells()
            Conns = CellArray.GetData()
            XConns = self.DataArrayToXdmfArray(Conns)
            print '%d Cells' % UGrid.GetNumberOfCells()
            NodesPerElement = XConns.GetValueAsInt64(0)
            print '%d NodesPerElement' % NodesPerElement
            XConns.SetShapeFromString("%d %d" % (UGrid.GetNumberOfCells(), NodesPerElement + 1))
            print 'Shape ', XConns.GetShapeAsString()
            Start = '0 1'
            Stride = '1 1'
            Count = '%d %d' % (UGrid.GetNumberOfCells(), NodesPerElement)
            XConns.SelectHyperSlabFromString(Start, Stride, Count)
            Xtop = Xgrid.GetTopology()
            Xtop.SetTopologyTypeFromString(self.CellTypes[UGrid.GetCellType(0)])
            Xtop.SetNumberOfElements(UGrid.GetNumberOfCells())
            H5 = XdmfHDF()
            H5.CopyType(XConns)
            H5.SetShapeFromString(Count)
            H5.Open('CORE:Data.h5:Group %05d/Index %05d/Conns' % (Group, Index), 'rw')
            H5.Write(XConns)
            JustConns = XdmfArray()
            JustConns.CopyType(H5)
            JustConns.CopyShape(H5)
            H5.Read(JustConns)
            Xtop.SetConnectivity(JustConns)
            H5.Close()
            self.Arrays.append(JustConns)
            
            # for i in range(UGrid.GetNumberOfCells()) :
            #     print '(%d) %d = %s' % (i, UGrid.GetCellType(i), self.CellTypes[UGrid.GetCellType(i)])
            print '%d Points' % UGrid.GetNumberOfPoints()
            pnts =  UGrid.GetPoints().GetData()
            # print "Points = ", pnts
            Xpnts = self.DataArrayToXdmfArray(pnts)
            Xgeo = Xgrid.GetGeometry()
            Xgeo.SetGeometryType(XDMF_GEOMETRY_XYZ)
            Xgeo.SetPoints(Xpnts)
            self.Arrays.append(Xpnts)
            # H5 = XdmfHDF()
            # H5.CopyType(Xpnts)
            # H5.CopyShape(Xpnts)
            # H5.Open('Data.h5:Group %05d/Index %05d/XYZ' % (Group, Index), 'rw')
            # H5.Write(Xpnts)
            # H5.Close()
            # for i in range(UGrid.GetNumberOfPoints()) :
            #     x, y, z = UGrid.GetPoint(i)
            #     print '(%d) %f %f %f' % (i, x, y, z)
            # print '%d Connections' % CellArray.GetNumberOfConnectivityEntries()
            # i = 0
            # while i < Conns.GetNumberOfTuples() :
            #     n = Conns.GetValue(i)
            #     i += 1
            #     for j in range(n) :
            #         p = Conns.GetValue(i)
            #         i += 1
            #         print '%d ' % p
            Pd = UGrid.GetPointData()
            print '# of Point Attributes = %d' %  Pd.GetNumberOfArrays()
            for i in range(Pd.GetNumberOfArrays()) :
                Pa = Pd.GetArray(i)
                Xpa = self.DataArrayToXdmfArray(Pa)
                Xattr = XdmfAttribute()
                Xattr.SetName("Point Attribute %d" %  i)
                Xattr.SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_NODE)
                Xattr.SetAttributeType(XDMF_ATTRIBUTE_TYPE_SCALAR)
                Xattr.SetValues(Xpa)
                Xgrid.Insert(Xattr)
                self.Arrays.append(Xpa)
                self.Arrays.append(Xattr)
                # self.Arrays.append(Xpa)
                # print Pd.GetArrayName(i), ' = ',Pa
                # for j in range(Pa.GetNumberOfTuples()) :
                   #  print '(%d) %f' % (j, Pa.GetValue(j))
            Cd = UGrid.GetCellData()
            print '# of Cell Attributes = %d' %  Cd.GetNumberOfArrays()
            for i in range(Cd.GetNumberOfArrays()) :
                Ca = Cd.GetArray(i)
                Xca = self.DataArrayToXdmfArray(Ca)
                Xattr = XdmfAttribute()
                Xattr.SetName("Cell Attribute %d" %  i)
                Xattr.SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_CELL)
                Xattr.SetAttributeType(XDMF_ATTRIBUTE_TYPE_SCALAR)
                Xattr.SetValues(Xca)
                Xgrid.Insert(Xattr)
                self.Arrays.append(Xca)
                self.Arrays.append(Xattr)
            # print UGrid
        else :
            print 'Heterogeneous UGrid'

    def ParseMultiGroup(self, Output) :
        print 'Parsing a vtkMultiGroupDataSet '
        NGroups = Output.GetNumberOfGroups()
        print 'Output has %d Groups' % NGroups
        for g in range(NGroups) :
            NDataSets = Output.GetNumberOfDataSets(g)
            print 'Group %d has %d DataSets' % (g + 1, NDataSets)
            for i in range(NDataSets) :
                ds = Output.GetDataSet(g,i)
                print 'Output Group %d Index %d (%s) has %d Cells' %  (g + 1, i + 1, ds.GetClassName(), ds.GetNumberOfCells())
                if ds.IsA('vtkUnstructuredGrid') :
                    self.ParseUGrid(ds)
                    self.WriteUGrid(ds, g, i)
                if ds.IsA('vtkImageData') :
                    self.WriteIGrid(ds, g, i)
                else :
                    print 'Can not handle vtk class = ', ds.GetClassName()

    def Parse(self, Output) :
        if Output.IsA('vtkMultiGroupDataSet') :
            self.ParseMultiGroup(Output)
        else :
            print 'Can not handle vtk class = ', Output.GetClassName()



if __name__ == '__main__' :
    # Reader = vtkXdmfReader()
    Reader = vtkXMLMultiGroupDataReader()
    # Reader = vtkXMLDataReader()
#    Controller = vtkMPIController()
#    Reader.SetController(Controller)
#    ProcId = Reader.GetController().GetLocalProcessId()
#    NumProcs = Reader.GetController().GetNumberOfProcesses()
#    print 'Hello from %d of %d' % (ProcId, NumProcs)
    # Reader.SetFileName('Points1.xmf')
    # Reader.SetFileName('VtkTest.pvd')
    Reader.SetFileName('VtkMulti.vtm')
    Reader.UpdateInformation()
    # Reader.DisableAllGrids()
    # Reader.EnableGrid(2)
    # Reader.DebugOn()
#    Reader.EnableAllGrids()
#    Reader.EnableAllArrays()
#    Reader.DisableGrid(0)
#    Reader.DisableGrid(1)
    Reader.Update()

    p = DataParser()
    print 'Reader has %d outputs' % Reader.GetNumberOfOutputPorts()
    for on in range(Reader.GetNumberOfOutputPorts()) :
        print 'Reading Output ', on
        Output = Reader.GetOutput(on)
        p.Parse(Output)
    p.Root.Build()
    print p.DOM.Serialize()
    p.DOM.Write('junk.xmf')

