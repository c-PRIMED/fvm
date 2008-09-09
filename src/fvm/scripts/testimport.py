import fvmbaseExt
import importers

fvmbaseExt.enableDebug("cdtor")

reader = importers.FluentReader("/home/sm/uvp.cas")

reader.readMesh();
