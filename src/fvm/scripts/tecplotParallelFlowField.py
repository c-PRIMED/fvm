import string
from mpi4py import MPI
from numpy import *
tectype = {
        'tri' : 'FETRIANGLE',
        'quad' : 'FEQUADRILATERAL',
        'tetra' : 'FETETRAHEDRON',
        'hexa' : 'FEBRICK'
        }

def dumpTecplotFile(nmesh, meshesLocal, meshesGlobal, mtype, flowFields, fname ):
  #cell sites
  cellSites = []
  for n in range(0,nmesh):
     cellSites.append( meshesGlobal[n].getCells() )

  #face sites
  faceSites = []
  for n in range(0,nmesh):
     faceSites.append( meshesGlobal[n].getFaces() )

  #node sites
  nodeSites = []
  for n in range(0,nmesh):
     nodeSites.append( meshesGlobal[n].getNodes() )

  #get connectivity (faceCells)
  faceCells = []
  for n in range(0,nmesh):
     faceCells.append( meshesGlobal[n].getConnectivity( faceSites[n], cellSites[n] ) )
 
  #get connectivity ( cellNodes )
  cellNodes = []
  for n in range(0,nmesh):
     cellNodes.append( meshesGlobal[n].getCellNodes() )

  cellSitesLocal = []
  for n in range(0,nmesh):
     cellSitesLocal.append( meshesLocal[n].getCells() )


  velFields = []
  for n in range(0,nmesh):
     velFields.append( flowFields.velocity[cellSitesLocal[n]].asNumPyArray() )
  densityFields = []
  
  for n in range(0,nmesh):
    densityFields.append( flowFields.density[cellSitesLocal[n]].asNumPyArray() )
  
  pressureFields = []
  for n in range(0,nmesh):
    pressureFields.append( flowFields.pressure[cellSitesLocal[n]].asNumPyArray() )
    
  coords = []
  for n in range(0,nmesh):
     coords.append( meshesGlobal[n].getNodeCoordinates().asNumPyArray() )


  #openning  global Array
  velFieldGlobal      = []
  densityFieldGlobal  = []
  pressureFieldGlobal = []
  for n in  range(0,nmesh):  
     selfCount = cellSites[n].getSelfCount()
     velFieldGlobal.append( zeros((selfCount,3), float) )
     densityFieldGlobal.append( zeros((selfCount,1), float) )
     pressureFieldGlobal.append( zeros((selfCount,1), float) )     
     
     meshLocal = meshesLocal[n]
     localToGlobal = meshLocal.getLocalToGlobalPtr().asNumPyArray()
     
     vFieldGlobal = velFieldGlobal[n]
     vFieldLocal  = velFields[n]

     dFieldGlobal = densityFieldGlobal[n]     
     dFieldLocal  = densityFields[n]
     
     pFieldGlobal = pressureFieldGlobal[n]
     pFieldLocal  = pressureFields[n]
     
     selfCountLocal = cellSitesLocal[n].getSelfCount()
     for i in range(0,selfCountLocal):
        globalID = localToGlobal[i]
	vFieldGlobal[globalID,:] = vFieldLocal[i,:]
	dFieldGlobal[globalID] = dFieldLocal[i]
	pFieldGlobal[globalID] = pFieldLocal[i]
     MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[vFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
     MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[dFieldGlobal, MPI.DOUBLE], op=MPI.SUM)    
     MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[pFieldGlobal, MPI.DOUBLE], op=MPI.SUM)   
      
  if MPI.COMM_WORLD.Get_rank() == 0:   
     f = open(fname, 'w')   
     f.write("Title = \" TECPLOT FILE FOR FLOWFIELDS \" \n")
     f.write("variables = \"x\", \"y\", \"z\", \"velX\", \"velY\", \"velZ\" , \"density\", \"pressure \"\n")
     for n in range(0,nmesh):
    	title_name = "nmesh" + str(n)
    	ncell  = cellSites[n].getSelfCount()
    	nnode  = nodeSites[n].getCount()
	f.write("Zone T = \"%s\" N = %s E = %s DATAPACKING = BLOCK, VARLOCATION = ([4-8]=CELLCENTERED), ZONETYPE=%s\n" %
            (title_name,  nodeSites[n].getCount(), ncell, tectype[mtype]))     
    	#write x	
	f.write("#x coordinates\n")
    	for i in range(0,nnode):
    	     f.write(str(coords[n][i][0])+"    ")
    	     if ( i % 5 == 4 ):
    		f.write("\n")
    	f.write("\n")	     

	f.write("#y coordinates\n")
    	#write y
    	for i in range(0,nnode):
    	     f.write(str(coords[n][i][1])+"    ")
    	     if ( i % 5 == 4 ):
    		f.write("\n")
    	f.write("\n")	     
	
	f.write("#z coordinates\n")
    	#write z
    	for i in range(0,nnode):
    	     f.write(str(coords[n][i][2])+"    ")
    	     if ( i % 5 == 4 ):
    		f.write("\n")
    	f.write("\n")	     
      
      	#write velX  	
	f.write("#velX \n")    	
    	for i in range(0,ncell):
    	   f.write( str(velFieldGlobal[n][i][0]) + "	")
    	   if ( i % 5  == 4 ):
    	       f.write("\n")
    	f.write("\n")
    		       
    	#write velY	
	f.write("#velY \n") 
    	for i in range(0,ncell):
    	   f.write( str(velFieldGlobal[n][i][1]) + "	")
    	   if ( i % 5  == 4 ):
    	       f.write("\n")
    	f.write("\n")	 
	    
    	#write velZ
	f.write("#velZ \n") 
    	for i in range(0,ncell):
    	   f.write( str(velFieldGlobal[n][i][2]) + "	")
    	   if ( i % 5  == 4 ):
    	       f.write("\n")
    	f.write("\n")
		       
        #write density
	f.write("#density \n") 	       
    	for i in range(0,ncell):
    	   f.write( str(densityFieldGlobal[n][i][0]) + "	    ")
    	   if ( i % 5  == 4 ):
    	       f.write("\n")
    	f.write("\n")	  
	
	#write pressure
	f.write("#pressure \n") 
    	for i in range(0,ncell):
    	   f.write( str(pressureFieldGlobal[n][i][0]) + "       ")
    	   if ( i % 5  == 4 ):
    	       f.write("\n")
    	f.write("\n")	   
    	       
    	#connectivity
	f.write("#connectivity \n")
    	for i in range(0,ncell):
    	   nnodes_per_cell = cellNodes[n].getCount(i)
    	   for node in range(0,nnodes_per_cell):
	       f.write( str(cellNodes[n](i,node)+1) + "	")
	   f.write("\n")	   
        f.close()

