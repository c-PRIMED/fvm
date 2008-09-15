import fvmbaseExt
import importers

fvmbaseExt.enableDebug("cdtor")

reader = importers.FluentReader("/home/sm/uvp.cas")

import debug
reader.readMesh();


meshes = reader.getMeshList()

metricsCalculator = fvmbaseExt.MeshMetricsCalculator_double(meshes)

metricsCalculator.init()
