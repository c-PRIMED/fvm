import fvmbaseExt
import importers
import models_atyped_double as models

#fvmbaseExt.enableDebug("cdtor")

reader = importers.FluentReader("/home/sm/uvp.cas")

#import debug
reader.readMesh();


meshes = reader.getMeshList()

geomFields =  models.GeomFields('geom')
metricsCalculator = models.MeshMetricsCalculatorA(geomFields,meshes)

metricsCalculator.init()

