import string

def esbgkTecplotFile(meshes, macroFields,filename ):
  #cell sites
  cellSites = []
  n=0
  cellSites.append( meshes[n].getCells() )
 # print "cellSites[", n, "].getCount = ", cellSites[n].getCount()

  #face sites
  faceSites = []
  faceSites.append( meshes[n].getFaces() )
  #print "faceSites[", n, "].getCount = ", faceSites[n].getCount()
  #node sites
  nodeSites = []
  nodeSites.append( meshes[n].getNodes() )
  #print "NodeSites[", n, "].getCount = ", nodeSites[n].getCount()
  #get connectivity (faceCells)
  faceCells = []
  faceCells.append( meshes[n].getConnectivity( faceSites[n], cellSites[n] ) )
 
  #get connectivity ( cellNodes )
  cellNodes = []
  cellNodes.append( meshes[n].getCellNodes() )

  velFields = []
  velFields.append( macroFields.velocity[cellSites[n]].asNumPyArray() )

  densityFields=[]
  densityFields.append(macroFields.density[cellSites[n]].asNumPyArray())
 
  pressureFields=[]
  pressureFields.append(macroFields.pressure[cellSites[n]].asNumPyArray())
  
  viscosityFields=[]
  viscosityFields.append(macroFields.viscosity[cellSites[n]].asNumPyArray())
   
  temperatureFields=[]
  temperatureFields.append(macroFields.temperature[cellSites[n]].asNumPyArray())
  
  collisionFrequencyFields=[]
  collisionFrequencyFields.append(macroFields.collisionFrequency[cellSites[n]].asNumPyArray())

  coords = []
  coords.append( meshes[n].getNodeCoordinates().asNumPyArray() )

  TxxFields=[]
  TxxFields.append(macroFields.Txx[cellSites[n]].asNumPyArray())  
  TyyFields=[]
  TyyFields.append(macroFields.Tyy[cellSites[n]].asNumPyArray())   
  TzzFields=[]
  TzzFields.append(macroFields.Tzz[cellSites[n]].asNumPyArray())
  TxyFields=[]
  TxyFields.append(macroFields.Txy[cellSites[n]].asNumPyArray())
  TxzFields=[]
  TxzFields.append(macroFields.Txz[cellSites[n]].asNumPyArray())
  TyzFields=[]
  TyzFields.append(macroFields.Tyz[cellSites[n]].asNumPyArray())
  #filename = "quadrature" + ".plt"
  f = open(filename, 'w')

  f.write("Title = \" tecplot out file\" \n")
  f.write("variables = \"x\", \"y\", \"z\", \"velX\", \"velY\", \"velZ\",\"density\",\"pressure\",\"viscosity\",\"temperature\", \"collisionFrequency\",\"Txx\",\"Tyy\",\"Tzz\",\"Txy\",\"Txz\",\"Tyz\",\n")
  #f.write("variables = \"x\", \"y\", \"z\", \"velX\", \"velY\", \"velZ\",\"density\",\"pressure\",\"viscosity\",\n")
  title_name = "nmesh" + str(n)
  ncell  = cellSites[n].getSelfCount()
  nnode  = nodeSites[n].getCount()
  zone_name = "Zone T = " + "\"" + title_name +  "\"" +      \
                 " N = " + str( nodeSites[n].getCount() ) +     \
		 " E = " + str( ncell ) +  \
		 " DATAPACKING = BLOCK, VARLOCATION = ([4-17]=CELLCENTERED), " + \
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
 
  #write density
  for i in range(0,ncell):
      f.write( str(densityFields[n][i]) + "    ")
      if ( i % 5  == 4 ):
        f.write("\n")
  f.write("\n")	
 
  #write pressure
  for i in range(0,ncell):
      f.write( str(pressureFields[n][i]) + "    ")
      if ( i % 5  == 4 ):
          f.write("\n")
  f.write("\n")	
  
  #write viscosity
  for i in range(0,ncell):
      f.write( str(viscosityFields[n][i]) + "    ")
      if ( i % 5  == 4 ):
          f.write("\n")
  f.write("\n")	 
#write temperature
  for i in range(0,ncell):
      f.write( str(temperatureFields[n][i]) + "    ")
      if ( i % 5  == 4 ):
          f.write("\n")
  f.write("\n")	 

#write collisionFrequency
  for i in range(0,ncell):
      f.write( str(collisionFrequencyFields[n][i]) + "    ")
      if ( i % 5  == 4 ):
          f.write("\n")
  f.write("\n")	 

  for i in range(0,ncell):
      f.write( str(TxxFields[n][i]) + "    ")
      if ( i % 5  == 4 ):
          f.write("\n")
  f.write("\n")	 

  for i in range(0,ncell):
      f.write( str(TyyFields[n][i]) + "    ")
      if ( i % 5  == 4 ):
          f.write("\n")
  f.write("\n")

  for i in range(0,ncell):
      f.write( str(TzzFields[n][i]) + "    ")
      if ( i % 5  == 4 ):
          f.write("\n")
  f.write("\n")

  for i in range(0,ncell):
      f.write( str(TxyFields[n][i]) + "    ")
      if ( i % 5  == 4 ):
          f.write("\n")
  f.write("\n")
  for i in range(0,ncell):
      f.write( str(TxzFields[n][i]) + "    ")
      if ( i % 5  == 4 ):
          f.write("\n")
  f.write("\n")
  for i in range(0,ncell):
      f.write( str(TyzFields[n][i]) + "    ")
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
