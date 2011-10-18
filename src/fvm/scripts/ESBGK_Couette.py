#!/usr/bin/env python
import sys

import fvm
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers
import fvm.exporters_atyped_double as exporters
import fvm.fvmparallel as fvmparallel

import time
from numpy import *
import tecplotESBGKEntireDomain

import string
from mpi4py import MPI

fvm.set_atype('double')

from FluentCase import FluentCase
#fvmbaseExt.enableDebug("cdtor")
from optparse import OptionParser
from Persistence import Persistence
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
if MPI.COMM_WORLD.Get_rank()==0:
    vFile = open("uerr.xy","w")

#user options
numTsteps=200
numIter=1
output_interval =200
output_Coeff=1000
timestep=0.5e-1;
cellno=10   #print coeffs for BGK and ESBGK &distribution function

fgamma=2; #0 max, 1 BGK, 2ES
method =1; # 1 constant, function otherwise
#uwall=0.0297;
ktrial=50;

nondimlength=0.1;
Kn=0.1;
accomCoeff=1.0;
rho_init=9.28e-6;
T_init=273.15;
relTol=1e-6;
absTol=1e-10;
tolx=1e-14;tolf=1e-14;
resint=0;
restartFileName = "test"+str(resint)+"_"+str(MPI.COMM_WORLD.Get_rank())+".hdf5" #restart from file
restartFileName = ""                                             #new run

gas='Argon'

if (gas == 'Air'):
    Pr=3.0/4.0;
    SpHeatRatio=7.0/5.0;
    muref= 1.7116e-5;
    mu_w=0.74;
    molecularWeight=28.9;
elif (gas == 'Argon'):
    Pr=2.0/3.0;
    SpHeatRatio=5.0/3.0;
    muref= 2.117e-5;
    mu_w=0.81;
    molecularWeight=39.9;
elif (gas == 'Nitrogen'):
    Pr=3.0/4.0;
    SpHeatRatio=7.0/5.0;
    muref= 1.781e-5;
    mu_w=0.74;
    molecularWeight=28.0;
        
    
print ' Pr= ',Pr,' muref= ',muref    
u_init=(2.0*8314.0*T_init/molecularWeight)**0.5
uwall=10.0/u_init;
if restartFileName != "":
    restartFile = Persistence(restartFileName,'r')
else:
    restartFile = None
# to test quadrature
import fvm.esbgk_atyped_double as esbgk

#facegroups
fgs = meshes[0].getBoundaryFaceGroups()
for fg in fgs:
    #print fg.groupType
    if  fg.id == 4 or fg.id == 6:
        fg.groupType = "realwall"
    print fg.id,"  ",fg.groupType

#Quadrature
quad0=esbgk.QuadratureD(10,10,10,5.5,1.0) #cartesian
#quad0=esbgk.QuadratureD(14,14,14,5.5,1.0) #cartesian
#quad0=esbgk.QuadratureD(8,8,0,32,0,16) #spherical
#quad2=esbgk.QuadratureD(16,16,1,8,1,4) #gauss-hermit quadrature and 3/8th rule

macroFields=esbgk.MacroFields('flow')
esbgk1=esbgk.KineticModelD(meshes,geomFields,macroFields,quad0)
esbgk1options = esbgk1.getOptions()

#initialize feq
esbgk1options.fgamma=fgamma;  #0 max, 1 for BGK,2  ES otherwise
esbgk1options.printCellNumber=cellno;
esbgk1options.NewtonsMethod_ktrial=ktrial;
esbgk1options.Prandtl=Pr;
esbgk1options.SpHeatRatio=SpHeatRatio;
esbgk1options.setVar('ToleranceX',tolx);
esbgk1options.setVar('ToleranceF',tolf);
esbgk1options.setVar('rho_init',rho_init);
esbgk1options.setVar('T_init',T_init);
esbgk1options.Prandtl=Pr;
esbgk1options.SpHeatRatio=SpHeatRatio;
esbgk1options.molecularWeight=molecularWeight;
esbgk1options.muref=muref;
esbgk1options.mu_w=mu_w;

esbgk1options.transient = False
esbgk1options.CentralDifference = False

esbgk1options.setVar('nonDimLt',nondimlength)
#Method 1:
#tSolver = fvmbaseExt.JacobiSolver()

#Method 2:
pc = fvmbaseExt.AMG()
pc.verbosity=0
tSolver = fvmbaseExt.BCGStab()
tSolver.preconditioner=pc
tSolver.verbosity=0

#Method 3: comment out the line below
esbgk1options.KineticLinearSolver=tSolver


# initialize macroparameters and f
esbgk1.InitializeMacroparameters()

cellSites=meshes[0].getCells()
cellCentroid=geomFields.coordinate[cellSites].asNumPyArray()
Velocity=macroFields.velocity[cellSites].asNumPyArray()
XVelocity=Velocity[:,0]
density =  macroFields.density[cellSites].asNumPyArray()
temperature =  macroFields.temperature[cellSites].asNumPyArray()
nuField = macroFields.collisionFrequency[cellSites].asNumPyArray()
coeffg = macroFields.coeffg[cellSites].asNumPyArray()
coeff= macroFields.coeff[cellSites].asNumPyArray()


ncells=cellSites.getSelfCount()
ncellsall=cellSites.getCount()
print ncellsall

def savevtk(n):
    writer = exporters.VTKWriterA(geomFields,meshes,"Couette"+string.zfill(str(n+1),5)+"_"+string.zfill(MPI.COMM_WORLD.Get_rank(),3)+".vtk","fix-fix beam",False,0)
    writer.init()
    writer.writeVectorField(geomFields.coordinate,"coordinates")
    writer.writeScalarField(macroFields.density,"density")
    writer.writeVectorField(macroFields.velocity,"velocity")
    writer.writeScalarField(macroFields.pressure,"pressure")
    writer.writeScalarField(macroFields.temperature,"temperature")
    writer.finish()


uana=zeros((ncells))
mACp2bAC=(2.0/accomCoeff-1.0)
for i in range(0,ncells):
    uana[i]=uwall*(cellCentroid[i][1]+mACp2bAC*Kn)/(1+2.0*mACp2bAC*Kn)
if restartFile is None:
    if(method ==1):#method 1
        esbgk1.weightedMaxwellian(0.5,0.5*uwall,0.5*uwall,1.0,1.0) #initial distribution
    else:  #method 2
        for i in range(0,ncellsall):
            XVelocity[i]=uwall*cellCentroid[i][1]
        esbgk1.initializeMaxwellian()
   
dsfList = []
ndir = quad0.getDirCount()
print 'no. of direntions ', ndir

DistFunc=esbgk1.getdsf()
for i in range(0,ndir):
    dsfList.append( DistFunc.getField(i) )
if esbgk1options.transient:
    DistFunc1=esbgk1.getdsf1()
    for i in range(0,ndir):
        dsfList.append( DistFunc1.getField(i) )
    if  esbgk1options.timeDiscretizationOrder > 1:
        DistFunc2=esbgk1.getdsf2()
        for i in range(0,ndir):
            dsfList.append( DistFunc2.getField(i) )

if restartFile is not None:
    restartFile.readKineticModel(macroFields,esbgk1,meshes,dsfList,ndir)
    restartFile.close()
    print "read f,fgamma,macropr from restart file"
    #esbgk1.initializeMaxwellian()


tecplotESBGKEntireDomain.esbgkTecplotEntireDomain(1,meshes,meshes_fluent,options.type,macroFields,"dump.dat")
esbgk1.ComputeMacroparameters()
if (fgamma==0):
    esbgk1.initializeMaxwellianEq()
else:
    esbgk1.EquilibriumDistributionBGK()
if(fgamma==2):
    esbgk1.EquilibriumDistributionESBGK()

#collision frequency based on Prandlt
esbgk1.ComputeCollisionfrequency()


bcMap = esbgk1.getBCMap()
if 4 in bcMap:
    bcTop = bcMap[4]
    bcTop.bcType = 'RealWallBC'
    bcTop.setVar('specifiedTemperature',1.0)
    bcTop.setVar('specifiedXVelocity',uwall)
    bcTop.setVar('accommodationCoefficient',accomCoeff)
if 6 in bcMap:
    bcBot = bcMap[6]
    bcBot.bcType = 'RealWallBC'
    bcBot.setVar('accommodationCoefficient',accomCoeff)
    #bcBot.setVar('specifiedTemperature',1.0)
    #bcBot.setVar('specifiedXVelocity',0.0)
if 3 in bcMap:
    bcLeft = bcMap[3]
    bcLeft.bcType = 'ZeroGradBC'
    #bcLeft.setVar('specifiedTemperature',2.0)
if 5 in bcMap:
    bcRight = bcMap[5]
    bcRight.bcType = 'ZeroGradBC'
    #bcRight.setVar('specifiedTemperature',2.0)

vcMap=esbgk1.getVCMap()
        
esbgk1options.setVar('timeStep',timestep)
esbgk1options.relativeTolerance=relTol;
esbgk1options.absoluteTolerance=absTol;
   
esbgk1.callBoundaryConditions()
#savevtk(-1)
ntc=array([0.0])
ntc1=array([1.0])
ntc=1.0*ncells+0.0*ntc1
MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[ntc, MPI.DOUBLE], op=MPI.SUM)

#initial profile
for cn in range(cellno,cellno*2):
    print cn,XVelocity[cn],uana[cn]
def calculate_rms(globaltimestep):
    fsum=array([0.0])
    for cn in range(0,ncells):
        fsum=fsum+(1-XVelocity[cn]/uana[cn])**2.0 
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[fsum, MPI.DOUBLE], op=MPI.SUM)
    #ntc=array([0.0])
    #ntc=1.0*ncells+0.0*fsum
    #MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[ntc, MPI.DOUBLE], op=MPI.SUM)
    if MPI.COMM_WORLD.Get_rank() == 0:
        fsum_global = 100.0*(fsum/ntc)**0.5
        print 'fsum',fsum_global, ' ncells',ntc
        vFile.write("%e %e\n" % (globaltimestep+1,fsum_global))
        vFile.flush()
if esbgk1options.transient:
    esbgk1.updateTime()
def advance(ntstep):
    for i in range(0,ntstep):
        print 'step = ',i+1
        
        esbgk1.advance(numIter)
        if esbgk1options.transient:
            print 'converged';
            esbgk1.updateTime()
     
        calculate_rms(i)
        if ((i+1)%output_interval == 0) :
           
            #dsfname = "output_"+string.zfill(str(i+1),5)+"_"+str(MPI.COMM_WORLD.Get_rank())+".dat"
            #esbgk1.OutputDsfBLOCK(dsfname)
            esbgk1.EntropyGeneration()
            filename = "macro_"+string.zfill(str(i+1),5)+".dat"
            tecplotESBGKEntireDomain.esbgkTecplotEntireDomain(1,meshes,meshes_fluent,options.type,macroFields,filename)
            #savevtk(i)
            
           
          
#esbgk1.OutputDsfBLOCK("output_00000.dat")
esbgk1.EntropyGeneration()
tecplotESBGKEntireDomain.esbgkTecplotEntireDomain(1,meshes,meshes_fluent,options.type,macroFields,"initial.dat")       
advance(numTsteps)        
esbgk1.EntropyGeneration()
tecplotESBGKEntireDomain.esbgkTecplotEntireDomain(1,meshes,meshes_fluent,options.type,macroFields,"final.dat")       
for cn in range(cellno,cellno*2):
    print cn,XVelocity[cn],uana[cn]
t1=time.time()
print 'time taken for',numTsteps,'tsteps',numIter,'iters =', t1-t0

"""
f = Persistence('test'+str(resint+1)+'_'+str(MPI.COMM_WORLD.Get_rank())+'.hdf5','w')
f.saveMeshes(meshes)
f.saveKineticModel(macroFields,esbgk1,meshes,dsfList,ndir)
f.close()
"""
if MPI.COMM_WORLD.Get_rank()==0:
    vFile.close()


