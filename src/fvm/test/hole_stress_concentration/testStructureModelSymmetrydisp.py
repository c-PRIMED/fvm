#!/usr/bin/env python

import sys
sys.setdlopenflags(0x100|0x2)

import fvm
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
fvm.set_atype('double')

import math

if fvm.atype == 'double':
    import fvm.models_atyped_double as models
    import fvm.exporters_atyped_double as exporters
elif fvm.atype == 'tangent':
    import fvm.models_atyped_tangent_double as models
    import fvm.exporters_atyped_tangent_double as exporters
    
from FluentCase import FluentCase
from optparse import OptionParser
                        
fileBase = None
numIterations = 100
fileBase = "/home/ba01/u119/das5/memosa/src/fvm/test/shankha/structure/circularHole/hole10"
#fileBase = "/home/sm/app-memosa/src/fvm/test/cav32"
#fileBase = "/home/sm/a/data/wj"

def usage():
    print "Usage: %s filebase [outfilename]" % sys.argv[0]
    print "Where filebase.cas is a Fluent case file."
    print "Output will be in filebase-prism.dat if it is not specified."
    sys.exit(1)

def advance(smodel,niter):
    for i in range(0,niter):
        try:
            smodel.advance(1)
        except KeyboardInterrupt:
            break


def createBV7Fields(geomFields,meshes,id):
    fxx7 = fvmbaseExt.Field('bvxx7')
    fxy7 = fvmbaseExt.Field('bvxy7')
    fyy7 = fvmbaseExt.Field('bvyy7')
    fxz7 = fvmbaseExt.Field('bvxz7')
    fyz7 = fvmbaseExt.Field('bvyz7')
    
    mesh = meshes[0]
    vol = geomFields.volume[mesh.getCells()]
    
    for mesh in meshes:
        fgs = mesh.getBoundaryGroups()
        for fg in fgs:
            if fg.id==7:                                
                nFaces = fg.site.getCount()
            
                sigmaXX = vol.newSizedClone(nFaces)
                sigmaXY = vol.newSizedClone(nFaces)
                sigmaYY = vol.newSizedClone(nFaces)
                sigmaXZ = vol.newSizedClone(nFaces)
                sigmaYZ = vol.newSizedClone(nFaces)

                sigmaXXa = sigmaXX.asNumPyArray()
                sigmaXYa = sigmaXY.asNumPyArray()
                sigmaYYa = sigmaYY.asNumPyArray()
                sigmaXZa = sigmaXZ.asNumPyArray()
                sigmaYZa = sigmaYZ.asNumPyArray()

                xf =  geomFields.coordinate[fg.site].asNumPyArray()
                a = 0.5
                tractionX = 10000.

                for i in range(0,nFaces):
                    x = xf[i][0]
                    y = xf[i][1]
                    r = math.sqrt(x*x+y*y)
                    theta = math.atan(y/x)

                    term1 = 1.0
                    term2 = ((a*a)/(r*r))*((3./2.)*math.cos(2.*theta)+math.cos(4.*theta))
                    term3 = (3./2.)*(math.pow(a,4.)/math.pow(r,4.))*math.cos(4.*theta)

                    term4 = ((a*a)/(r*r))*(0.5*math.cos(2.*theta)-math.cos(4.*theta))
                    term5 = 1.5*(math.pow(a,4.)/math.pow(r,4.))*math.cos(4.*theta)

                    term6 = ((a*a)/(r*r))*(0.5*math.sin(2.*theta)+math.sin(4.*theta))
                    term7 = (3./2.)*(math.pow(a,4.)/math.pow(r,4.))*math.sin(4.*theta)

                    sigmaXXa[i] = tractionX*(term1-term2+term3)
                    sigmaYYa[i] = tractionX*(-term4-term5)
                    sigmaXYa[i] = tractionX*(-term6+term7)
                    sigmaXZa[i] = 0.0
                    sigmaYZa[i] = 0.0
                    
                fxx7[fg.site] = sigmaXX
                fxy7[fg.site] = sigmaXY
                fyy7[fg.site] = sigmaYY
                fxz7[fg.site] = sigmaXZ
                fyz7[fg.site] = sigmaYZ
                return fxx7,fxy7,fyy7,fxz7,fyz7

def createBV6Fields(geomFields,meshes,id):
    fxx6 = fvmbaseExt.Field('bvxx6')
    fxy6 = fvmbaseExt.Field('bvxy6')
    fyy6 = fvmbaseExt.Field('bvyy6')
    fxz6 = fvmbaseExt.Field('bvxz6')
    fyz6 = fvmbaseExt.Field('bvyz6')
    
    mesh = meshes[0]
    vol = geomFields.volume[mesh.getCells()]
    
    for mesh in meshes:
        fgs = mesh.getBoundaryGroups()
        for fg in fgs:
            if fg.id==6:
                nFaces = fg.site.getCount()
                
                sigmaXX = vol.newSizedClone(nFaces)
                sigmaXY = vol.newSizedClone(nFaces)
                sigmaYY = vol.newSizedClone(nFaces)
                sigmaXZ = vol.newSizedClone(nFaces)
                sigmaYZ = vol.newSizedClone(nFaces)
                
                sigmaXXa = sigmaXX.asNumPyArray()
                sigmaXYa = sigmaXY.asNumPyArray()
                sigmaYYa = sigmaYY.asNumPyArray()
                sigmaXZa = sigmaXZ.asNumPyArray()
                sigmaYZa = sigmaYZ.asNumPyArray()
                
                xf =  geomFields.coordinate[fg.site].asNumPyArray()
                a = 0.5
                tractionX = 10000.

                for i in range(0,nFaces):
                    x = xf[i][0]
                    y = xf[i][1]
                    r = math.sqrt(x*x+y*y)
                    theta = math.atan(y/x)
                    
                    term1 = 1.0
                    term2 = ((a*a)/(r*r))*((3./2.)*math.cos(2.*theta)+math.cos(4.*theta))
                    term3 = (3./2.)*(math.pow(a,4.)/math.pow(r,4.))*math.cos(4.*theta)
                    
                    term4 = ((a*a)/(r*r))*(0.5*math.cos(2.*theta)-math.cos(4.*theta))
                    term5 = 1.5*(math.pow(a,4.)/math.pow(r,4.))*math.cos(4.*theta)
                    
                    term6 = ((a*a)/(r*r))*(0.5*math.sin(2.*theta)+math.sin(4.*theta))
                    term7 = (3./2.)*(math.pow(a,4.)/math.pow(r,4.))*math.sin(4.*theta)
                    
                    sigmaXXa[i] = tractionX*(term1-term2+term3)
                    sigmaYYa[i] = tractionX*(-term4-term5)
                    sigmaXYa[i] = tractionX*(-term6+term7)
                    sigmaXZa[i] = 0.0
                    sigmaYZa[i] = 0.0

                fxx6[fg.site] = sigmaXX
                fxy6[fg.site] = sigmaXY
                fyy6[fg.site] = sigmaYY
                fxz6[fg.site] = sigmaXZ
                fyz6[fg.site] = sigmaYZ
                return fxx6,fxy6,fyy6,fxz6,fyz6
                
def createBV5Fields(geomFields,meshes,id):
    fxx5 = fvmbaseExt.Field('bvxx5')
    fxy5 = fvmbaseExt.Field('bvxy5')
    fyy5 = fvmbaseExt.Field('bvyy5')
    fxz5 = fvmbaseExt.Field('bvxz5')
    fyz5 = fvmbaseExt.Field('bvyz5')
    
    mesh = meshes[0]
    vol = geomFields.volume[mesh.getCells()]
    
    for mesh in meshes:
        fgs = mesh.getBoundaryGroups()
        for fg in fgs:
            if fg.id==5:
                nFaces = fg.site.getCount()
                
                sigmaXX = vol.newSizedClone(nFaces)
                sigmaXY = vol.newSizedClone(nFaces)
                sigmaYY = vol.newSizedClone(nFaces)
                sigmaXZ = vol.newSizedClone(nFaces)
                sigmaYZ = vol.newSizedClone(nFaces)
                
                sigmaXXa = sigmaXX.asNumPyArray()
                sigmaXYa = sigmaXY.asNumPyArray()
                sigmaYYa = sigmaYY.asNumPyArray()
                sigmaXZa = sigmaXZ.asNumPyArray()
                sigmaYZa = sigmaYZ.asNumPyArray()
                
                xf =  geomFields.coordinate[fg.site].asNumPyArray()
                a = 0.5
                tractionX = 10000.
                
                for i in range(0,nFaces):
                    x = xf[i][0]
                    y = xf[i][1]
                    r = math.sqrt(x*x+y*y)
                    theta = math.atan(y/x)
                    
                    term1 = 1.0
                    term2 = ((a*a)/(r*r))*((3./2.)*math.cos(2.*theta)+math.cos(4.*theta))
                    term3 = (3./2.)*(math.pow(a,4.)/math.pow(r,4.))*math.cos(4.*theta)
                    
                    term4 = ((a*a)/(r*r))*(0.5*math.cos(2.*theta)-math.cos(4.*theta))
                    term5 = 1.5*(math.pow(a,4.)/math.pow(r,4.))*math.cos(4.*theta)
                    term6 = ((a*a)/(r*r))*(0.5*math.sin(2.*theta)+math.sin(4.*theta))
                    term7 = (3./2.)*(math.pow(a,4.)/math.pow(r,4.))*math.sin(4.*theta)
                    
                    sigmaXXa[i] = tractionX*(term1-term2+term3)
                    sigmaYYa[i] = tractionX*(-term4-term5)
                    sigmaXYa[i] = tractionX*(-term6+term7)
                    sigmaXZa[i] = 0.0
                    sigmaYZa[i] = 0.0
                    
                fxx5[fg.site] = sigmaXX
                fxy5[fg.site] = sigmaXY
                fyy5[fg.site] = sigmaYY
                fxz5[fg.site] = sigmaXZ
                fyz5[fg.site] = sigmaYZ
                return fxx5,fxy5,fyy5,fxz5,fyz5
                    
def createBV4Fields(geomFields,meshes,id):
    fxx4 = fvmbaseExt.Field('bvxx4')
    fxy4 = fvmbaseExt.Field('bvxy4')
    fyy4 = fvmbaseExt.Field('bvyy4')
    fxz4 = fvmbaseExt.Field('bvxz4')
    fyz4 = fvmbaseExt.Field('bvyz4')
    
    mesh = meshes[0]
    vol = geomFields.volume[mesh.getCells()]
    
    for mesh in meshes:
        fgs = mesh.getBoundaryGroups()
        for fg in fgs:
            if fg.id==4:
                nFaces = fg.site.getCount()
                
                sigmaXX = vol.newSizedClone(nFaces)
                sigmaXY = vol.newSizedClone(nFaces)
                sigmaYY = vol.newSizedClone(nFaces)
                sigmaXZ = vol.newSizedClone(nFaces)
                sigmaYZ = vol.newSizedClone(nFaces)

                sigmaXXa = sigmaXX.asNumPyArray()
                sigmaXYa = sigmaXY.asNumPyArray()
                sigmaYYa = sigmaYY.asNumPyArray()
                sigmaXZa = sigmaXZ.asNumPyArray()
                sigmaYZa = sigmaYZ.asNumPyArray()
                
                xf =  geomFields.coordinate[fg.site].asNumPyArray()
                a = 0.5
                tractionX = 10000.

                for i in range(0,nFaces):
                    x = xf[i][0]
                    y = xf[i][1]
                    r = math.sqrt(x*x+y*y)
                    theta = math.atan(y/x)
                    
                    term1 = 1.0
                    term2 = ((a*a)/(r*r))*((3./2.)*math.cos(2.*theta)+math.cos(4.*theta))
                    term3 = (3./2.)*(math.pow(a,4.)/math.pow(r,4.))*math.cos(4.*theta)
                    
                    term4 = ((a*a)/(r*r))*(0.5*math.cos(2.*theta)-math.cos(4.*theta))
                    term5 = 1.5*(math.pow(a,4.)/math.pow(r,4.))*math.cos(4.*theta)
                    term6 = ((a*a)/(r*r))*(0.5*math.sin(2.*theta)+math.sin(4.*theta))
                    term7 = (3./2.)*(math.pow(a,4.)/math.pow(r,4.))*math.sin(4.*theta)
                    
                    sigmaXXa[i] = tractionX*(term1-term2+term3)
                    sigmaYYa[i] = tractionX*(-term4-term5)
                    sigmaXYa[i] = tractionX*(-term6+term7)
                    sigmaXZa[i] = 0.0
                    sigmaYZa[i] = 0.0
                    
                fxx4[fg.site] = sigmaXX
                fxy4[fg.site] = sigmaXY
                fyy4[fg.site] = sigmaYY
                fxz4[fg.site] = sigmaXZ
                fyz4[fg.site] = sigmaYZ
                return fxx4,fxy4,fyy4,fxz4,fyz4

def createBV3Fields(geomFields,meshes,id):
    fxx3 = fvmbaseExt.Field('bvxx3')
    fxy3 = fvmbaseExt.Field('bvxy3')
    fyy3 = fvmbaseExt.Field('bvyy3')
    fxz3 = fvmbaseExt.Field('bvxz3')
    fyz3 = fvmbaseExt.Field('bvyz3')
    
    mesh = meshes[0]
    vol = geomFields.volume[mesh.getCells()]
    
    for mesh in meshes:
        fgs = mesh.getBoundaryGroups()
        for fg in fgs:
            if fg.id==3:
                nFaces = fg.site.getCount()
                
                sigmaXX = vol.newSizedClone(nFaces)
                sigmaXY = vol.newSizedClone(nFaces)
                sigmaYY = vol.newSizedClone(nFaces)
                sigmaXZ = vol.newSizedClone(nFaces)
                sigmaYZ = vol.newSizedClone(nFaces)
                
                sigmaXXa = sigmaXX.asNumPyArray()
                sigmaXYa = sigmaXY.asNumPyArray()
                sigmaYYa = sigmaYY.asNumPyArray()
                sigmaXZa = sigmaXZ.asNumPyArray()
                sigmaYZa = sigmaYZ.asNumPyArray()
                
                xf =  geomFields.coordinate[fg.site].asNumPyArray()
                a = 0.5
                tractionX = 10000.
                
                for i in range(0,nFaces):
                    x = xf[i][0]
                    y = xf[i][1]
                    r = math.sqrt(x*x+y*y)
                    theta = math.atan(y/x)
                    
                    term1 = 1.0
                    term2 = ((a*a)/(r*r))*((3./2.)*math.cos(2.*theta)+math.cos(4.*theta))
                    term3 = (3./2.)*(math.pow(a,4.)/math.pow(r,4.))*math.cos(4.*theta)
                    
                    term4 = ((a*a)/(r*r))*(0.5*math.cos(2.*theta)-math.cos(4.*theta))
                    term5 = 1.5*(math.pow(a,4.)/math.pow(r,4.))*math.cos(4.*theta)
                    term6 = ((a*a)/(r*r))*(0.5*math.sin(2.*theta)+math.sin(4.*theta))
                    term7 = (3./2.)*(math.pow(a,4.)/math.pow(r,4.))*math.sin(4.*theta)
                    
                    sigmaXXa[i] = tractionX*(term1-term2+term3)
                    sigmaYYa[i] = tractionX*(-term4-term5)
                    sigmaXYa[i] = tractionX*(-term6+term7)
                    sigmaXZa[i] = 0.0
                    sigmaYZa[i] = 0.0
                    
                fxx3[fg.site] = sigmaXX
                fxy3[fg.site] = sigmaXY
                fyy3[fg.site] = sigmaYY
                fxz3[fg.site] = sigmaXZ
                fyz3[fg.site] = sigmaYZ
                return fxx3,fxy3,fyy3,fxz3,fyz3


# change as needed

# map between fvm, tecplot, and xdmf types
etype = {
    'tri' : 1,
    'quad' : 2,
    'tetra' : 3,
    'hexa' : 4
        }
tectype = {
    'tri' : 'FETRIANGLE',
    'quad' : 'FEQUADRILATERAL',
    'tetra' : 'FETETRAHEDRON',
    'hexa' : 'FEBRICK'
        }

def dumpTecplotFile(nmesh, meshes, geomFields, mtype):
    #cell sites
    cellSites = []
    for n in range(0,nmesh):
        cellSites.append( meshes[n].getCells() )
#     print "cellSites[", n, "].getCount = ", cellSites[n].getCount()
        
    #face sites
    faceSites = []
    for n in range(0,nmesh):
        faceSites.append( meshes[n].getFaces() )
        
    #node sites
    nodeSites = []
    for n in range(0,nmesh):
        nodeSites.append( meshes[n].getNodes() )
        
    #get connectivity (faceCells)
    faceCells = []
    for n in range(0,nmesh):
        faceCells.append( meshes[n].getConnectivity( faceSites[n], cellSites[n] ) )
        
    #get connectivity ( cellNodes )
    cellNodes = []
    for n in range(0,nmesh):
        cellNodes.append( meshes[n].getCellNodes() )
        
    #get Volume as array
    volumes = []
    for n in range(0,nmesh):
        volumes.append( geomFields.volume[cellSites[n]].asNumPyArray() )
    
    
    cellCentroids =[]
    for n in range(0,nmesh):
        cellCentroids.append( geomFields.coordinate[cellSites[n]].asNumPyArray() )
        
    defFields = []
    for n in range(0,nmesh):
        defFields.append( structureFields.deformation[cellSites[n]].asNumPyArray() )

    tractionXFields = []
    for n in range(0,nmesh):
        tractionXFields.append( structureFields.tractionX[cellSites[n]].asNumPyArray() )
        
    coords = []
    for n in range(0,nmesh):
        coords.append( geomFields.coordinate[nodeSites[n]].asNumPyArray() )
#     print "shape( coords[", n, "] ) = ", shape( coords[n] )

    f = open("tecplot_hole10_symmetry.dat","w")
    f.write("Title = \" tecplot file for 2D Cavity problem \" \n")
    f.write("variables = \"x\", \"y\", \"z\", \"defX\", \"defY\", \"sigmaXX\", \"sigmaXY\", \"sigmaYY\", \"cellCentroidY\" \n")
    for n in range(0,nmesh):
        title_name = "nmesh%s" % n
        ncell  = cellSites[n].getSelfCount()
        nnode  = nodeSites[n].getCount()
        f.write("Zone T = \"%s\" N = %s E = %s DATAPACKING = BLOCK, VARLOCATION = ([4-9]=CELLCENTERED), ZONETYPE=%s\n" %
                (title_name,  nodeSites[n].getCount(), ncell, tectype[mtype]))
        #write x
        for i in range(0,nnode):
            f.write(str(coords[n][i][0])+"    ")
            if ( i % 5 == 4 ):
                f.write("\n")
        f.write("\n")
                
        #write y
        for i in range(0,nnode):
            f.write(str(coords[n][i][1])+"    ")
            if ( i % 5 == 4 ):
                f.write("\n")
        f.write("\n")

        #write z
        for i in range(0,nnode):
            f.write(str(coords[n][i][2])+"    ")
            if ( i % 5 == 4 ):
                f.write("\n")
        f.write("\n")

        #write defX
        for i in range(0,ncell):
            f.write(str(defFields[n][i][0]) + "    ")
            if ( i % 5  == 4 ):
                f.write("\n")
        f.write("\n")

        #write defY
        for i in range(0,ncell):
            f.write(str(defFields[n][i][1]) + "    ")
            if ( i % 5  == 4 ):
                f.write("\n")
        f.write("\n")

        #write sigmaXX
        for i in range(0,ncell):
            f.write(str(tractionXFields[n][i][0]) + "    ")
            if ( i % 5  == 4 ):
                f.write("\n")
        f.write("\n")

        #write sigmaXY
        for i in range(0,ncell):
            f.write(str(tractionXFields[n][i][1]) + "    ")
            if ( i % 5  == 4 ):
                f.write("\n")
        f.write("\n")
                                                    
        #write sigmaYY
        for i in range(0,ncell):
            f.write(str(tractionXFields[n][i][2]) + "    ")
            if ( i % 5  == 4 ):
                f.write("\n")
        f.write("\n")
                                                                    
        #write velX
        for i in range(0,ncell):
            f.write( str(cellCentroids[n][i][1]) + "    ")
            if ( i % 5  == 4 ):
                f.write("\n")
        f.write("\n")

        #connectivity
        for i in range(0,ncell):
            nnodes_per_cell = cellNodes[n].getCount(i)
            for node in range(0,nnodes_per_cell):
                f.write( str(cellNodes[n](i,node)+1) + "     ")
            f.write("\n")
        f.write("\n")
    f.close()

parser = OptionParser()
parser.set_defaults(type='quad')
parser.add_option("--type", help="'quad'[default], 'tri', 'hexa', or 'tetra'")
parser.add_option("--xdmf", action='store_true', help="Dump data in xdmf")
parser.add_option("--time","-t",action='store_true',help="Print timing information.")
(options, args) = parser.parse_args()

outfile = None
if __name__ == '__main__' and fileBase is None:
    if len(sys.argv) < 2:
        usage()
    fileBase = sys.argv[1]
    if len(sys.argv) == 3:
        outfile = sys.argv[2]

if outfile == None:
    outfile = fileBase+"-prism.dat"
    
reader = FluentCase(fileBase+".cas")

#import debug
reader.read();

meshes = reader.getMeshList()
mesh = meshes[0]
nmesh = 1
   
import time
t0 = time.time()

geomFields =  models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

fgs = mesh.getBoundaryGroups()
for fg in fgs:
    if fg.id in [6,7]:
        fg.groupType='symmetry'

metricsCalculator.init()

cells = mesh.getCells()

rho = 7854.0
E = 1.0*math.pow(10,7)
nu = 0.3

if fvm.atype == 'tangent':
    metricsCalculator.setTangentCoords(0,7,1)

structureFields =  models.StructureFields('structure')
smodel = models.StructureModelA(geomFields,structureFields,meshes)

bcMap = smodel.getBCMap()

#left
bcID = 6
if bcID in bcMap:
    sigmaXX6,sigmaXY6,sigmaYY6,sigmaXZ6,sigmaYZ6 = createBV6Fields(geomFields,meshes,bcID)
    bc = smodel.getBCMap()[bcID]
    bc.bcType = 'Symmetry'
    bc['specifiedXXTraction']=sigmaXX6
    bc['specifiedXYTraction']=sigmaXY6
    bc['specifiedXZTraction']=sigmaXZ6
    bc['specifiedYXTraction']=sigmaXY6
    bc['specifiedYYTraction']=sigmaYY6
    bc['specifiedYZTraction']=sigmaYZ6
    bc['specifiedZXTraction']=0
    bc['specifiedZYTraction']=0
    bc['specifiedZZTraction']=0

#right 
bcID = 5
if bcID in bcMap:
    bc = smodel.getBCMap()[bcID]
    sigmaXX5,sigmaXY5,sigmaYY5,sigmaXZ5,sigmaYZ5 = createBV5Fields(geomFields,meshes,bcID)
    bc.bcType = 'SpecifiedTraction'
    bc['specifiedXXTraction']=sigmaXX5
    bc['specifiedXYTraction']=sigmaXY5
    bc['specifiedXZTraction']=sigmaXZ5
    bc['specifiedYXTraction']=sigmaXY5
    bc['specifiedYYTraction']=sigmaYY5
    bc['specifiedYZTraction']=sigmaYZ5
    bc['specifiedZXTraction']=0
    bc['specifiedZYTraction']=0
    bc['specifiedZZTraction']=0

# bottom
bcID = 7
if bcID in bcMap:
    bc = smodel.getBCMap()[bcID]
    sigmaXX7,sigmaXY7,sigmaYY7,sigmaXZ7,sigmaYZ7 = createBV7Fields(geomFields,meshes,bcID)
    bc.bcType = 'Symmetry'
    bc['specifiedXXTraction']=sigmaXX7
    bc['specifiedXYTraction']=sigmaXY7
    bc['specifiedXZTraction']=sigmaXZ7
    bc['specifiedYXTraction']=sigmaXY7
    bc['specifiedYYTraction']=sigmaYY7
    bc['specifiedYZTraction']=sigmaYZ7
    bc['specifiedZXTraction']=0
    bc['specifiedZYTraction']=0
    bc['specifiedZZTraction']=0

# top
bcID = 4
if bcID in bcMap:
    bc = smodel.getBCMap()[bcID]
    sigmaXX4,sigmaXY4,sigmaYY4,sigmaXZ4,sigmaYZ4 = createBV4Fields(geomFields,meshes,bcID)
    bc.bcType = 'SpecifiedTraction'
    bc['specifiedXXTraction']=sigmaXX4
    bc['specifiedXYTraction']=sigmaXY4
    bc['specifiedXZTraction']=sigmaXZ4
    bc['specifiedYXTraction']=sigmaXY4
    bc['specifiedYYTraction']=sigmaYY4
    bc['specifiedYZTraction']=sigmaYZ4
    bc['specifiedZXTraction']=0
    bc['specifiedZYTraction']=0
    bc['specifiedZZTraction']=0
    
#hole
bcID = 3
if bcID in bcMap:
    bc = smodel.getBCMap()[bcID]
    sigmaXX3,sigmaXY3,sigmaYY3,sigmaXZ3,sigmaYZ3 = createBV3Fields(geomFields,meshes,bcID)
    bc.bcType = 'SpecifiedTraction'
    bc['specifiedXXTraction']=sigmaXX3
    bc['specifiedXYTraction']=sigmaXY3
    bc['specifiedXZTraction']=sigmaXZ3
    bc['specifiedYXTraction']=sigmaXY3
    bc['specifiedYYTraction']=sigmaYY3
    bc['specifiedYZTraction']=sigmaYZ3
    bc['specifiedZXTraction']=0.0
    bc['specifiedZYTraction']=0.0
    bc['specifiedZZTraction']=0.0
                                        

vcMap = smodel.getVCMap()
for i,vc in vcMap.iteritems():
    vc['density'] = rho
    vc['eta'] = E/(2.*(1+nu))
    vc['eta1'] = nu*E/((1+nu)*(1-2.0*nu))

pc = fvmbaseExt.AMG()
pc.verbosity=0
defSolver = fvmbaseExt.BCGStab()
defSolver.preconditioner = pc
defSolver.relativeTolerance = 1e-9
defSolver.absoluteTolerance = 1e-30
defSolver.nMaxIterations = 5000
#defSolver.maxCoarseLevels=20
defSolver.verbosity=0

#defSolver = fvmbaseExt.AMG()
#defSolver.relativeTolerance = 1e-1
#defSolver.nMaxIterations = 1000
#defSolver.maxCoarseLevels=20
#defSolver.verbosity=1


soptions = smodel.getOptions()
soptions.deformationLinearSolver = defSolver
soptions.deformationTolerance=1e-9
soptions.setVar("deformationURF",1.0)
soptions.printNormalizedResiduals=False
soptions.transient=False

"""
if fvm.atype=='tangent':
    vcMap = fmodel.getVCMap()
    for i,vc in vcMap.iteritems():
        print vc.getVar('viscosity')
        vc.setVar('viscosity',(1.7894e-5,1))
"""
smodel.init()
smodel.printBCs
smodel.advance(numIterations)
smodel.getTractionX(mesh)

t1 = time.time()
if outfile != '/dev/stdout':
    print '\nsolution time = %f' % (t1-t0)


faceCells = mesh.getAllFaceCells()
deformation =  structureFields.deformation[mesh.getCells()].asNumPyArray()
fileName = fileBase + "displacement10_left.txt"
file = open(fileName,"w")
file.write("displacement\t\n")
for mesh in meshes:
    fgs = mesh.getBoundaryGroups()
    for fg in fgs:
        nFaces = fg.site.getCount()
        
        xf =  geomFields.coordinate[fg.site].asNumPyArray()
        a = 0.5
        k = 3.-4.*nu
#        k = (3.-nu)/(1.+nu)
        mu = E/(2.*(1.+nu))
        term=10000./(4.*mu)
        if fg.id==6:
            faceCells = mesh.getFaceCells(fg.site)            
            for i in range(0,nFaces):
                wx = deformation[faceCells(i,0)][0]
                wy = deformation[faceCells(i,0)][1]
                wmag=math.sqrt(wx*wx+wy*wy)
                x = xf[i][0]
                y = xf[i][1]
                r = math.sqrt(x*x+y*y)
                theta = math.atan(y/x)
                term1 = r*((k-1)/2.+math.cos(2.*theta))
                term2 = (a*a/r)*(1.+(1.+k)*math.cos(2.*theta))
                term3 = (math.pow(a,4.)/math.pow(r,3.))*math.cos(2.*theta)
                def0 = (term1+term2-term3)*term
                def1 = 0.
                def2 = 0.
                dmag = math.sqrt(def0*def0+def1*def1)
                file.write(" %e " % xf[i][0])
                file.write(" %e " % xf[i][1])
                file.write(" %e " % wx)
                file.write(" %e " % wy)
                file.write(" %e " % def0)
                file.write(" %e " % def1)
                file.write("\n")
file.close()

dumpTecplotFile( nmesh, meshes, geomFields, options.type)




