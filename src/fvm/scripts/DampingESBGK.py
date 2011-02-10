#!/usr/bin/env python
import sys

import fvm
import fvm.fvmbaseExt as fvmbaseExt
import fvm.importers as importers

import time
from numpy import *
import tecplotESBGK
import tecplot
import string

fvm.set_atype('double')

from FluentCase import FluentCase
#fvmbaseExt.enableDebug("cdtor")

fileBase = None

fileBase = "/home/ba01/u116/schigull/memosa/src/fvm/test/DampingESBGK/Damping100x100"
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
reader = FluentCase(fileBase+".cas")

reader.read();
t0 = time.time()

meshes = reader.getMeshList()

geomFields =  fvm.models.GeomFields('geom')
metricsCalculator = fvm.models.MeshMetricsCalculatorA(geomFields,meshes)
metricsCalculator.init()

#user options
numTsteps = 200
numIter=20
output_interval =10
output_Coeff=100
timestep=0.5e-2;
cellno=8   #print coeffs for BGK and ESBGK &distribution function

fgamma=2; #0 max, 1 BGK, 2ES
method =1; # 1 constant, linear otherwise
uwall=0.0297*0.2;  #5m/s
ktrial=500;
Pr=2.0/3.0;
SpHeatRatio=5.0/3.0;
nondimlength=1.4e-6;
rho_init=1.78e-2; #0.01atmosphere pressure
tolx=1e-14;tolf=1e-14;

# to test quadrature
import fvm.esbgk_atyped_double as esbgk

#Quadrature- velocity mesh
quad0=esbgk.QuadratureD(10,10,10,5.5,1.0) #cartesian 10by10by10
#quad0=esbgk.QuadratureD(15,15,15,5.5,1.0) #cartesian 15by15by15
#quad1=esbgk.QuadratureD(2,2,0,12,0,5) #spherical
#quad2=esbgk.QuadratureD(16,16,1,8,1,4) #gauss-hermit quadrature and 3/8th rule

#import ddd
macroFields=esbgk.MacroFields('flow')

cellSites=meshes[0].getCells()
ncells=cellSites.getSelfCount()
ncellsall=cellSites.getCount()
print ncellsall


esbgk1=esbgk.KineticModelD(meshes,geomFields,macroFields,quad0)
esbgk1options = esbgk1.getOptions()

## initialize macroparameters and f
esbgk1.InitializeMacroparameters()

esbgk1.weightedMaxwellian(0.5,0*0.5,0.5*0,1.0,1.0) #initial distribution

        
#initialize feq

esbgk1options.fgamma=fgamma;  #0 max, 1 for BGK,2  ES otherwise
esbgk1options.printCellNumber=cellno;
esbgk1options.NewtonsMethod_ktrial=ktrial;
esbgk1options.Prandtl=Pr;
esbgk1options.SpHeatRatio=SpHeatRatio;
esbgk1options.setVar('ToleranceX',tolx);
esbgk1options.setVar('ToleranceF',tolf);
esbgk1options.setVar('rho_init',rho_init);

esbgk1options.setVar('nonDimLength',nondimlength);
print "initializing macropr and feq"
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

print "setting bcs"

bcMap = esbgk1.getBCMap()
bcTop = bcMap[4]
bcTop.bcType = 'PressureInletBC'
bcTop.setVar('specifiedTemperature',1.0)
bcTop.setVar('specifiedPressure',1.0)

bcBot = bcMap[6]
bcBot.bcType = 'WallBC'
bcBot.setVar('specifiedTemperature',1.0)

bcLeft = bcMap[3]
bcLeft.bcType = 'SymmetryBC'

bcRight = bcMap[5]
bcRight.bcType = 'PressureInletBC'
bcRight.setVar('specifiedTemperature',1.0)
bcRight.setVar('specifiedPressure',1.0)

bcBeam=bcMap[7]
bcBeam.bcType = 'WallBC'
bcBeam.setVar('specifiedTemperature',1.0)
bcBeam.setVar('specifiedYVelocity',uwall)


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
            dsfname = "output_"+string.zfill(str(i+1),5)+".dat"
           #esbgk1.OutputDsfBLOCK(dsfname)
            filename = "macro_"+string.zfill(str(i+1),5)+".dat"
            tecplotESBGK.esbgkTecplotFile(meshes,macroFields,filename)


esbgk1.OutputDsfBLOCK("output_00000.dat")
tecplotESBGK.esbgkTecplotFile(meshes,macroFields,"macro_00000.dat")       
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
