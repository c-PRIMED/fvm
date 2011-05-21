import string
from mpi4py import MPI
from numpy import *

tectype = {
        'tri' : 'FETRIANGLE',
        'quad' : 'FEQUADRILATERAL',
        'tetra' : 'FETETRAHEDRON',
        'hexa' : 'FEBRICK'
        }

def esbgkTecplotParallelFile(nmesh, meshesLocal, meshesGlobal, mtype, macroFields,filename ):
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

  #coords
  coords = []
  for n in range(0,nmesh):
     coords.append( meshesGlobal[n].getNodeCoordinates().asNumPyArray() )
 
  cellSitesLocal = []
  for n in range(0,nmesh):
     cellSitesLocal.append( meshesLocal[n].getCells() )


  velFields = []
  for n in range(0,nmesh):
      velFields.append( macroFields.velocity[cellSitesLocal[n]].asNumPyArray() )

  densityFields=[]
  for n in range(0,nmesh):
      densityFields.append(macroFields.density[cellSitesLocal[n]].asNumPyArray())
 
  pressureFields=[]
  for n in range(0,nmesh):
      pressureFields.append(macroFields.pressure[cellSitesLocal[n]].asNumPyArray())
  
  viscosityFields=[]
  for n in range(0,nmesh):
      viscosityFields.append(macroFields.viscosity[cellSitesLocal[n]].asNumPyArray())
   
  temperatureFields=[]
  for n in range(0,nmesh):
      temperatureFields.append(macroFields.temperature[cellSitesLocal[n]].asNumPyArray())
  
  collisionFrequencyFields=[]
  for n in range(0,nmesh):
      collisionFrequencyFields.append(macroFields.collisionFrequency[cellSitesLocal[n]].asNumPyArray())

  
  TxxFields=[]
  for n in range(0,nmesh):
      TxxFields.append(macroFields.Txx[cellSitesLocal[n]].asNumPyArray())  
  TyyFields=[]
  for n in range(0,nmesh):
      TyyFields.append(macroFields.Tyy[cellSitesLocal[n]].asNumPyArray())   
  TzzFields=[]
  for n in range(0,nmesh):
      TzzFields.append(macroFields.Tzz[cellSitesLocal[n]].asNumPyArray())
  TxyFields=[]
  for n in range(0,nmesh):
      TxyFields.append(macroFields.Txy[cellSitesLocal[n]].asNumPyArray())
  TxzFields=[]
  for n in range(0,nmesh):
      TxzFields.append(macroFields.Txz[cellSitesLocal[n]].asNumPyArray())
  TyzFields=[]
  for n in range(0,nmesh):
      TyzFields.append(macroFields.Tyz[cellSitesLocal[n]].asNumPyArray())

  #opening global Array
  densityFieldGlobal=[]
  velFieldGlobal = []
  pressureFieldGlobal = []
  viscosityFieldGlobal = []
  temperatureFieldGlobal = []
  collisionFrequencyFieldGlobal = []
  TxxFieldGlobal=[]
  TyyFieldGlobal=[]
  TzzFieldGlobal=[]
  TxyFieldGlobal=[]
  TxzFieldGlobal=[]
  TyzFieldGlobal=[]
  
  for n in range(0,nmesh):
    selfCount = cellSites[n].getSelfCount()

    #append global arrays
    densityFieldGlobal.append( zeros((selfCount,1), float) )
    velFieldGlobal.append( zeros((selfCount,3), float) ) #if it is velocity (selfCount,3)
    pressureFieldGlobal.append( zeros((selfCount,1), float) )
    viscosityFieldGlobal.append( zeros((selfCount,1), float) )
    temperatureFieldGlobal.append( zeros((selfCount,1), float) )
    collisionFrequencyFieldGlobal.append( zeros((selfCount,1), float) ) 
    TxxFieldGlobal.append( zeros((selfCount,1), float) )
    TyyFieldGlobal.append( zeros((selfCount,1), float) )
    TzzFieldGlobal.append( zeros((selfCount,1), float) )
    TxyFieldGlobal.append( zeros((selfCount,1), float) )
    TxzFieldGlobal.append( zeros((selfCount,1), float) )
    TyzFieldGlobal.append( zeros((selfCount,1), float) )
    
    meshLocal = meshesLocal[n]
    localToGlobal = meshLocal.getLocalToGlobalPtr().asNumPyArray()

    #re-lable
    dFieldGlobal  = densityFieldGlobal[n]
    dFieldLocal   = densityFields[n]
    vFieldGlobal  = velFieldGlobal[n]
    vFieldLocal   = velFields[n]    
    tFieldGlobal  = temperatureFieldGlobal[n]
    tFieldLocal   = temperatureFields[n]
    pFieldGlobal  = pressureFieldGlobal[n]
    pFieldLocal   = pressureFields[n]
    visFieldGlobal  = viscosityFieldGlobal[n]
    visFieldLocal   = viscosityFields[n]
    cfFieldGlobal  = collisionFrequencyFieldGlobal[n]
    cfFieldLocal   = collisionFrequencyFields[n]
    txxFieldGlobal  = TxxFieldGlobal[n]
    txxFieldLocal   = TxxFields[n]
    tyyFieldGlobal  = TyyFieldGlobal[n]
    tyyFieldLocal   = TyyFields[n]
    tzzFieldGlobal  = TzzFieldGlobal[n]
    tzzFieldLocal   = TzzFields[n]
    txyFieldGlobal  = TxyFieldGlobal[n]
    txyFieldLocal   = TxyFields[n]
    txzFieldGlobal  = TxzFieldGlobal[n]
    txzFieldLocal   = TxzFields[n]
    tyzFieldGlobal  = TyzFieldGlobal[n]
    tyzFieldLocal   = TyzFields[n]

    #fill local part of cells
    selfCount = cellSitesLocal[n].getSelfCount()
    for i in range(0,selfCount):
       globalID               = localToGlobal[i]
       dFieldGlobal[globalID] = dFieldLocal[i]
       vFieldGlobal[globalID] = vFieldLocal[i]
       tFieldGlobal[globalID] = tFieldLocal[i]
       pFieldGlobal[globalID] = pFieldLocal[i]
       visFieldGlobal[globalID] = visFieldLocal[i]
       cfFieldGlobal[globalID] = cfFieldLocal[i]
       txxFieldGlobal[globalID] = txxFieldLocal[i]
       tyyFieldGlobal[globalID] = tyyFieldLocal[i]
       tzzFieldGlobal[globalID] = tzzFieldLocal[i]
       txyFieldGlobal[globalID] = txyFieldLocal[i]
       txzFieldGlobal[globalID] = txzFieldLocal[i]
       tyzFieldGlobal[globalID] = tyzFieldLocal[i]
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[dFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[vFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[tFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[pFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[visFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[cfFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[txxFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[tyyFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[tzzFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[txyFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[txzFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[tyzFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
  if MPI.COMM_WORLD.Get_rank() == 0:     
  #filename = "quadrature" + ".plt"
    f = open(filename, 'w')

    f.write("Title = \" tecplot parallel file\" \n")
    f.write("variables = \"x\", \"y\", \"z\", \"velX\", \"velY\", \"velZ\",\"density\",\"pressure\",\"viscosity\",\"temperature\", \"collisionFrequency\",\"Txx\",\"Tyy\",\"Tzz\",\"Txy\",\"Txz\",\"Tyz\",\n")
    #f.write("variables = \"x\", \"y\", \"z\", \"Density\",\"velX\", \"velY\", \"velZ\",\"Temperature\",\n")
    for n in range(0,nmesh):
      title_name = "nmesh" + str(n)
      ncell  = cellSites[n].getSelfCount()
      nnode  = nodeSites[n].getCount()
      f.write("Zone T = \"%s\" N = %s E = %s DATAPACKING = BLOCK, VARLOCATION = ([4-17]=CELLCENTERED), ZONETYPE=%s\n" %
             (title_name,  nodeSites[n].getCount(), ncell, tectype[mtype]))     

 
      #zone_name = "Zone T = " + "\"" + title_name +  "\"" +      \
      #               " N = " + str( nodeSites[n].getCount() ) +     \
      #               " E = " + str( ncell ) +  \
      #               " DATAPACKING = BLOCK, VARLOCATION = ([4-17]=CELLCENTERED), " + \
      #               " ZONETYPE=FEQUADRILATERAL \n"
      #f.write( zone_name )
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
    
      
      #write velX,Y,Z
      vFieldGlobal = velFieldGlobal[n]
      for i in range(0,ncell):
          f.write( str(vFieldGlobal[i][0]) + "    ")
          if ( i % 5  == 4 ):
              f.write("\n")
      f.write("\n")                  
      #write velY
      for i in range(0,ncell):
          f.write( str(vFieldGlobal[i][1]) + "    ")
          if ( i % 5  == 4 ):
              f.write("\n")
      f.write("\n")	  
      #write velZ
      for i in range(0,ncell):
          f.write( str(vFieldGlobal[i][2]) + "    ")
          if ( i % 5  == 4 ):
            f.write("\n")
      f.write("\n")	

      #write density
      dFieldGlobal = densityFieldGlobal[n]
      for i in range(0,ncell):
         f.write( str(dFieldGlobal[i][0]) + "    ")
         if ( i % 5  == 4 ):
            f.write("\n")
      f.write("\n")
      
           
      #write pressure
      pFieldGlobal = pressureFieldGlobal[n]
      for i in range(0,ncell):
          f.write( str(pFieldGlobal[i][0]) + "    ")
          if ( i % 5  == 4 ):
              f.write("\n")
      f.write("\n")	
      
      #write viscosity
      visFieldGlobal = viscosityFieldGlobal[n]
      for i in range(0,ncell):
          f.write( str(visFieldGlobal[i][0]) + "    ")
          if ( i % 5  == 4 ):
              f.write("\n")
      f.write("\n")
      
      #write temperature (again) 
      tFieldGlobal = temperatureFieldGlobal[n]
      for i in range(0,ncell):
         f.write( str(tFieldGlobal[i][0]) + "    ")
         if ( i % 5  == 4 ):
            f.write("\n")
      f.write("\n")
      
      #write collision frequency
      cfFieldGlobal = collisionFrequencyFieldGlobal[n]
      for i in range(0,ncell):
         f.write( str(cfFieldGlobal[i][0]) + "    ")
         if ( i % 5  == 4 ):
            f.write("\n")
      f.write("\n")

      #write Txx 
      txxFieldGlobal = TxxFieldGlobal[n]
      for i in range(0,ncell):
         f.write( str(txxFieldGlobal[i][0]) + "    ")
         if ( i % 5  == 4 ):
            f.write("\n")
      f.write("\n")
      #write Tyy 
      tyyFieldGlobal = TyyFieldGlobal[n]
      for i in range(0,ncell):
         f.write( str(tyyFieldGlobal[i][0]) + "    ")
         if ( i % 5  == 4 ):
            f.write("\n")
      f.write("\n")
      #write Tzz 
      tzzFieldGlobal = TzzFieldGlobal[n]
      for i in range(0,ncell):
         f.write( str(tzzFieldGlobal[i][0]) + "    ")
         if ( i % 5  == 4 ):
            f.write("\n")
      f.write("\n")
      
      #write Txy 
      txyFieldGlobal = TxyFieldGlobal[n]
      for i in range(0,ncell):
         f.write( str(txyFieldGlobal[i][0]) + "    ")
         if ( i % 5  == 4 ):
            f.write("\n")
      f.write("\n")
      #write Txz 
      txzFieldGlobal = TxzFieldGlobal[n]
      for i in range(0,ncell):
         f.write( str(txzFieldGlobal[i][0]) + "    ")
         if ( i % 5  == 4 ):
            f.write("\n")
      f.write("\n")
      #write Tyz 
      tyzFieldGlobal = TyzFieldGlobal[n]
      for i in range(0,ncell):
         f.write( str(tyzFieldGlobal[i][0]) + "    ")
         if ( i % 5  == 4 ):
            f.write("\n")
      f.write("\n")
      
      #connectivity
      for i in range(0,ncell):
         nnodes_per_cell = cellNodes[n].getCount(i)
         for node in range(0,nnodes_per_cell):
            f.write( str(cellNodes[n](i,node)+1) + "     ")
         f.write("\n")
      
      f.close()


"""
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
 """ 
 
