import string
from mpi4py import MPI

tectype = {
        'tri' : 'FETRIANGLE',
        'quad' : 'FEQUADRILATERAL',
        'tetra' : 'FETETRAHEDRON',
        'hexa' : 'FEBRICK'
        }

def dumpTecplotMesh(nmesh, meshes, mtype):
  #cell sites
  cellSites = []
  for n in range(0,nmesh):
     cellSites.append( meshes[n].getCells() )

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

  #coords
  coords = []
  for n in range(0,nmesh):
     coords.append( meshes[n].getNodeCoordinates().asNumPyArray() )
     
	   
  file_name = "mesh_PROC" +  str( MPI.COMM_WORLD.Get_rank() ) + ".dat"
  f = open(file_name, 'w')

  f.write("Title = \" tecplot mesh \" \n")
  f.write("variables = \"x\", \"y\", \"z\" \n")
  for n in range(0,nmesh):
     title_name = "nmesh" + str(n)
     ncellSelf  = cellSites[n].getSelfCount()
     ncell      = cellSites[n].getCount()
     nnode  = nodeSites[n].getCount()
     f.write("Zone T = \"%s\" N = %s E = %s DATAPACKING = BLOCK, ZONETYPE=%s\n" %
             (title_name,  nnode, ncell, tectype[mtype]))
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
     
     
     #connectivity
     for i in range(0,ncellSelf):
        nnodes_per_cell = cellNodes[n].getCount(i)
        for node in range(0,nnodes_per_cell):
	    f.write( str(cellNodes[n](i,node)+1) + "     ")
        f.write("\n")
	
     #connectivity (ghost cells)
     nnodes = cellNodes[n].getCount(0);
     for i in range(ncellSelf,ncell):
        nnodes_per_cell = cellNodes[n].getCount(i)
        for node in range(0,nnodes):
	    f.write( str(cellNodes[n](i,node%nnodes_per_cell)+1) + "     ")
        f.write("\n")
	

  f.close()

