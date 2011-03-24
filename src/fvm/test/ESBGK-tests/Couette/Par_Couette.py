#!/usr/bin/env python
import sys

import fvm
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.fvmparallel as fvmparallel

import time
from numpy import *
import tecplotESBGK
import tecplot
import string
from mpi4py import MPI

fvm.set_atype('double')

from FluentCase import FluentCase
#fvmbaseExt.enableDebug("cdtor")
from optparse import OptionParser

#fvmbaseExt.enableDebug("cdtor")
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
xtype = {
        'tri' : 'Triangle',
        'quad' : 'Quadrilateral',
        'tetra' : 'Tetrahedron',
        'hexa' : 'Hexahedron'
        }

fileBase = None

fileBase = "/home/ba01/u116/schigull/memosa/ESBGK-tests/Feb2011/testSymmetry/PD_40by40"
def usage():
    print "Usage: %s filebase [outfilename]" % sys.argv[0]
    print "Where filebase.cas is a Fluent case file."
    print "Output will be in filebase-prism.dat if it is not specified."
    sys.exit(1)


# change as needed

outfile = None
if __name__ == '__main__' and fileBase is None:
    if len(sys.argv) < 2:
        usage()
    fileBase = sys.argv[1]
    if len(sys.argv) == 3:
        outfile = sys.argv[2]

if outfile == None:
    outfile = fileBase+"-tecplt.dat"

#import debug
parser = OptionParser()
parser.set_defaults(type='tri')
parser.add_option("--type", help="'tri'[default], 'quad', 'hexa', or 'tetra'")
parser.add_option("--xdmf", action='store_true', help="Dump data in xdmf")
parser.add_option("--time","-t",action='store_true',help="Print timing information.")
(options, args) = parser.parse_args()
if len(args) != 1:
    usage()

reader = FluentCase(args[0])
#reader = FluentCase(fileBase+".cas")

reader.read();
t0 = time.time()
meshes_fluent = reader.getMeshList()
nmesh = 1
npart = [MPI.COMM_WORLD.Get_size()]
etype = [etype[options.type]]

if not MPI.COMM_WORLD.Get_rank():
   print "parmesh is processing"

#Partitinoing of spatial mesh
#partMesh constructor and setTypes
part_mesh = fvmparallel.MeshPartitioner( meshes_fluent, npart, etype );
part_mesh.setWeightType(0);
part_mesh.setNumFlag(0);

#actions
part_mesh.partition()
part_mesh.mesh()
#part_mesh.mesh_debug()
meshes = part_mesh.meshList()
#meshes = reader.getMeshList()

geomFields =  fvm.models.GeomFields('geom')
metricsCalculator = fvm.models.MeshMetricsCalculatorA(geomFields,meshes)
metricsCalculator.init()

#user options
numTsteps = 100
numIter=1
output_interval =50
output_Coeff=100
timestep=0.5e-1;
cellno=8   #print coeffs for BGK and ESBGK &distribution function
fgam_Eq=0; #0 for BGK, ES otherwise
fgamma=2; #0 max, 1 BGK, 2ES
method =1; # 1 constant, linear otherwise
uwall=0.0297;
ktrial=500;
Pr=3.0/3.0;
SpHeatRatio=5.0/3.0;

nondimlength=1.0;
rho_init=9.28e-7;
tolx=1e-14;tolf=1e-14;

# to test quadrature
import fvm.esbgk_atyped_double as esbgk

#Quadrature
quad0=esbgk.QuadratureD(10,10,10,5.5,1.0) #cartesian
#quad1=esbgk.QuadratureD(2,2,0,12,0,5) #spherical
#quad2=esbgk.QuadratureD(16,16,1,8,1,4) #gauss-hermit quadrature and 3/8th rule

macroFields=esbgk.MacroFields('flow')
esbgk1=esbgk.KineticModelD(meshes,geomFields,macroFields,quad0)
esbgk1options = esbgk1.getOptions()

# initialize macroparameters and f
esbgk1.InitializeMacroparameters()

cellSites=meshes[0].getCells()
cellCentroid=geomFields.coordinate[cellSites].asNumPyArray()
Velocity=macroFields.velocity[cellSites].asNumPyArray()
XVelocity=Velocity[:,0]
ncells=cellSites.getSelfCount()

ncellsall=cellSites.getCount()
print ncellsall


if(method ==1):#method 1
    esbgk1.weightedMaxwellian(0.5,uwall*0.5,0.5*uwall,1.0,1.0) #initial distribution
else:
    #method 2
    for i in range(0,ncells):
	XVelocity[i]=0.0;#cellCentroid[i][1]*uwall
    esbgk1.initializeMaxwellian()
        
#initialize feq
esbgk1options.ESBGK_fgamma=fgam_Eq;  #0 for BGK, ES otherwise
esbgk1options.fgamma=fgamma;  #0 max, 1 for BGK,2  ES otherwise
esbgk1options.printCellNumber=cellno;
esbgk1options.NewtonsMethod_ktrial=ktrial;
esbgk1options.Prandtl=Pr;
esbgk1options.SpHeatRatio=SpHeatRatio;
esbgk1options.setVar('ToleranceX',tolx);
esbgk1options.setVar('ToleranceF',tolf);
esbgk1options.setVar('rho_init',rho_init);

esbgk1options.transient = False

esbgk1options.setVar('nonDimLength',nondimlength)

esbgk1.ComputeMacroparameters()
if (fgamma==0):
    esbgk1.initializeMaxwellianEq()
else:
    esbgk1.EquilibriumDistributionBGK()
if(fgamma==2):
    esbgk1.EquilibriumDistributionESBGK()

#collision frequency based on Prandlt
esbgk1.ComputeCollisionfrequency()

cellSite = meshes[0].getCells()
densityField =  macroFields.density[cellSite].asNumPyArray()
nuField = macroFields.collisionFrequency[cellSite].asNumPyArray()

bcMap = esbgk1.getBCMap()
if 4 in bcMap:
    bcTop = bcMap[4]
    bcTop.bcType = 'WallBC'
    bcTop.setVar('specifiedTemperature',1.0)
    bcTop.setVar('specifiedXVelocity',uwall)
if 6 in bcMap:
    bcBot = bcMap[6]
    bcBot.bcType = 'WallBC'
    bcBot.setVar('specifiedTemperature',1.0)
    bcBot.setVar('specifiedXVelocity',0.0)
if 3 in bcMap:
    bcLeft = bcMap[3]
    bcLeft.bcType = 'CopyBC'
    #bcLeft.setVar('specifiedTemperature',2.0)
if 5 in bcMap:
    bcRight = bcMap[5]
    bcRight.bcType = 'CopyBC'
    #bcRight.setVar('specifiedTemperature',2.0)

vcMap=esbgk1.getVCMap()
        
esbgk1options.setVar('timeStep',timestep)
esbgk1options.relativeTolerance=1e-7;
esbgk1options.absoluteTolerance=1e-10;
   
esbgk1.callBoundaryConditions()

#esbgk1.updateTime()
def advance(ntstep):
    for i in range(0,ntstep):
        print 'timestep = ',i+1
      
        esbgk1.advance(numIter)
        print 'converged';
        esbgk1.updateTime()
        
        if ((i+1)%output_Coeff == 0) :
            if(fgamma>0):
                coeff=macroFields.coeff[cellSites].asNumPyArray()
                print 'BGK:',coeff[cellno,0],'cx^2',coeff[cellno,1],'cx',coeff[cellno,2]            
                
            if(fgamma==2):
                coeffg=macroFields.coeffg[cellSites].asNumPyArray()
                print 'ESBGK:',coeffg[cellno,0],'cx^2',coeffg[cellno,1],'cx',coeffg[cellno,2]
                print '     :','cy^2',coeffg[cellno,3],'cy',coeffg[cellno,4],'cz^2',coeffg[cellno,5],'cz',coeffg[cellno,6]
                print 'cxcy',coeffg[cellno,7],'cxcz',coeffg[cellno,8],'cycz',coeffg[cellno,9]
                
        
        if ((i+1)%output_interval == 0) :
            """
            dens=macroFields.density[cellSites].asNumPyArray()
            print 'density',dens[105],dens[115],dens[125],dens[135]
            press=macroFields.pressure[cellSites].asNumPyArray()
            print 'pressure',press[105],press[115],press[125],press[135]
            """
            dsfname = "output_"+string.zfill(str(i+1),5)+"_"+str(MPI.COMM_WORLD.Get_rank())+".dat"
           #esbgk1.OutputDsfBLOCK(dsfname)
            filename = "macro_"+string.zfill(str(i+1),5)+"_"+str(MPI.COMM_WORLD.Get_rank())+".dat"
            tecplotESBGK.esbgkTecplotFile(meshes,macroFields,filename)


esbgk1.OutputDsfBLOCK("output_00000.dat")
tecplotESBGK.esbgkTecplotFile(meshes,macroFields,"macro_00000"+"_"+str(MPI.COMM_WORLD.Get_rank())+".dat")       
advance(numTsteps)        


t1=time.time()
print 'time taken for',numTsteps,'tsteps',numIter,'iters =', t1-t0

"""

cellSite = meshes[0].getCells()
velField = macroFields.velocity[ cellSite ].asNumPyArray() 
print  macroFields.velocity[ cellSite ].asNumPyArray()
rhoField = macroFields.density[ cellSite ].asNumPyArray() 
print  macroFields.density[ cellSite ].asNumPyArray()
#import debug
#fmodel.advance(100)
filename="macro1.plt"
tecplotESBGK.esbgkTecplotFile(meshes,macroFields,filename)
#tecplot.dumpTecplotFile(1,meshes,macroFields)
"""
