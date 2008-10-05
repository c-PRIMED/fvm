import sys
sys.setdlopenflags(0x100|0x2)

import fvmbaseExt
import importers
import models_atyped_double as models
import exporters_atyped_double as exporters
from FluentCase import FluentCase
#fvmbaseExt.enableDebug("cdtor")

reader = FluentCase("/home/sm/olda/data/cond.cas")

#import debug
reader.read();


meshes = reader.getMeshList()

geomFields =  models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

thermalFields =  models.ThermalFields('therm')

tmodel = models.ThermalModelA(geomFields,thermalFields,meshes)
bcmap = tmodel.getBCMap()
bc2 = bcmap[2]

reader.importThermalBCs(tmodel)
tmodel.printBCs()
tmodel.init()
tmodel.advance(1)

writer = exporters.FluentDataExporterA(reader,"test.dat",False,0)
writer.init()
writer.writeScalarField(thermalFields.temperature,3)
writer.finish()
