import string
#from mpi4py import MPI

def dumpTecplotFile(nmesh, meshes, flowFields ):
  #cell sites
  cellSites = []
  for n in range(0,nmesh):
     cellSites.append( meshes[n].getCells() )
     print "cellSites[", n, "].getCount = ", cellSites[n].getCount()

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

  velFields = []
  for n in range(0,nmesh):
     velFields.append( flowFields.velocity[cellSites[n]].asNumPyArray() )
  densityFields = []
  
  for n in range(0,nmesh):
    densityFields.append( flowFields.density[cellSites[n]].asNumPyArray() )
  
  pressureFields = []
  for n in range(0,nmesh):
    pressureFields.append( flowFields.pressure[cellSites[n]].asNumPyArray() )
  coords = []
  for n in range(0,nmesh):
     coords.append( meshes[n].getNodeCoordinates().asNumPyArray() )
     
  file_name = "quadrature_nmesh" + ".plt"
  f = open(file_name, 'w')   
  #file_name = "cavity" +  str( MPI.COMM_WORLD.Get_rank() ) + ".dat"
  #f = open(file_name, 'w')

  f.write("Title = \" tecplot file for 2D Cavity problem \" \n")
  f.write("variables = \"x\", \"y\", \"z\", \"velX\", \"velY\", \"velZ\" ,\"density\",\"pressure \",\n")
  for n in range(0,nmesh):
     title_name = "nmesh" + str(n)
     ncell  = cellSites[n].getSelfCount()
     nnode  = nodeSites[n].getCount()
     zone_name = "Zone T = " + "\"" + title_name +  "\"" +      \
                 " N = " + str( nodeSites[n].getCount() ) +     \
		 " E = " + str( ncell ) +  "\n"
     f.write(zone_name)
     zone_name= " DATAPACKING = BLOCK, VARLOCATION = ([4-8]=CELLCENTERED), " + \
		 " ZONETYPE=FEQUADRILATERAL \n"
     f.write( zone_name )
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
     
     
     #write velX
     for i in range(0,ncell):
        f.write( str(velFields[n][i][0]) + "    ")
	if ( i % 5  == 4 ):
	    f.write("\n")
     f.write("\n")
     	  	    
     #write velY
     for i in range(0,ncell):
        f.write( str(velFields[n][i][1]) + "    ")
	if ( i % 5  == 4 ):
	    f.write("\n")
     f.write("\n")	  
     #write velZ
     for i in range(0,ncell):
        f.write( str(velFields[n][i][2]) + "    ")
	if ( i % 5  == 4 ):
	    f.write("\n")
     f.write("\n")
     for i in range(0,ncell):
        f.write( str(densityFields[n][i]) + "    ")
	if ( i % 5  == 4 ):
	    f.write("\n")
     f.write("\n")	
     for i in range(0,ncell):
        f.write( str(pressureFields[n][i]) + "    ")
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

