import string
from mpi4py import MPI
from numpy import *

tectype = {
        'tri' : 'FETRIANGLE',
        'quad' : 'FEQUADRILATERAL',
        'tetra' : 'FETETRAHEDRON',
        'hexa' : 'FEBRICK'
        }

def esbgkTecplotEntireDomain(nmesh, meshesLocal, meshesGlobal, mtype, macroFields,filename ):
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

  entropyFields=[]
  for n in range(0,nmesh):
    entropyFields.append(macroFields.Entropy[cellSitesLocal[n]].asNumPyArray())
    
  entgenFields=[]
  for n in range(0,nmesh):
    entgenFields.append(macroFields.EntropyGenRate_Collisional[cellSitesLocal[n]].asNumPyArray())
  
  stressFields = []
  for n in range(0,nmesh):
      stressFields.append( macroFields.Stress[cellSitesLocal[n]].asNumPyArray() )
  KnqFields=[]
  for n in range(0,nmesh):
    KnqFields.append(macroFields.Knq[cellSitesLocal[n]].asNumPyArray())
    
  heatFluxFields = []
  for n in range(0,nmesh):
    heatFluxFields.append(macroFields.heatFlux[cellSitesLocal[n]].asNumPyArray())

  #opening global Array
  densityFieldGlobal=[]
  velFieldGlobal = []
  pressureFieldGlobal = []
  viscosityFieldGlobal = []
  temperatureFieldGlobal = []
  collisionFrequencyFieldGlobal = []
  stressFieldGlobal = []
  entropyFieldGlobal = []
  entgenFieldGlobal = []
  KnqFieldGlobal = []
  heatFluxFieldGlobal = []
  for n in range(0,nmesh):
    selfCount = cellSites[n].getSelfCount()

    #append global arrays
    densityFieldGlobal.append( zeros((selfCount,1), float) )
    velFieldGlobal.append( zeros((selfCount,3), float) ) #if it is velocity (selfCount,3)
    pressureFieldGlobal.append( zeros((selfCount,1), float) )
    viscosityFieldGlobal.append( zeros((selfCount,1), float) )
    temperatureFieldGlobal.append( zeros((selfCount,1), float) )
    collisionFrequencyFieldGlobal.append( zeros((selfCount,1), float) ) 
    stressFieldGlobal.append( zeros((selfCount,6), float) ) 
    entropyFieldGlobal.append( zeros((selfCount,1), float) )
    entgenFieldGlobal.append( zeros((selfCount,1), float) )
    KnqFieldGlobal.append( zeros((selfCount,1), float) )
    heatFluxFieldGlobal.append( zeros((selfCount,3), float) )
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
    sFieldGlobal  = stressFieldGlobal[n]
    sFieldLocal   = stressFields[n]    
    eFieldGlobal = entropyFieldGlobal[n]
    eFieldLocal = entropyFields[n]
    egFieldGlobal = entgenFieldGlobal[n]
    egFieldLocal = entgenFields[n]
    kFieldGlobal = KnqFieldGlobal[n]
    kFieldLocal = KnqFields[n]
    hfFieldGlobal = heatFluxFieldGlobal[n]
    hfFieldLocal = heatFluxFields[n] 
    
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
       sFieldGlobal[globalID] = sFieldLocal[i]
       eFieldGlobal[globalID] = eFieldLocal[i]
       egFieldGlobal[globalID] = egFieldLocal[i] 
       kFieldGlobal[globalID] = kFieldLocal[i]
       hfFieldGlobal[globalID] = hfFieldLocal[i]
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[dFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[vFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[tFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[pFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[visFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[cfFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[sFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[eFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[egFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[kFieldGlobal, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,[hfFieldGlobal, MPI.DOUBLE], op=MPI.SUM)

  if MPI.COMM_WORLD.Get_rank() == 0:     
  #filename = "quadrature" + ".plt"
    f = open(filename, 'w')

    f.write("Title = \" tecplot parallel file\" \n")
    f.write("variables = \"x\", \"y\", \"z\",\"velX\",\"velY\",\"velZ\",\"heatFlux_x\",\"heatFlux_y\",\"heatFlux_z\",\"density\",\"pressure\",\"viscosity\",\"temperature\", \"collisionFrequency\",\"Pxx\",\"Pyy\",\"Pzz\",\"Pxy\",\"Pyz\",\"Pzx\",\"Entropy\",\"EntGenRate\",\"Knq \",\n")


    #f.write("variables = \"x\", \"y\", \"z\", \"Density\",\"velX\", \"velY\", \"velZ\",\"Temperature\",\n")
    for n in range(0,nmesh):
      title_name = "nmesh" + str(n)
      ncell  = cellSites[n].getSelfCount()
      nnode  = nodeSites[n].getCount()
      f.write("Zone T = \"%s\" N = %s E = %s DATAPACKING = BLOCK, VARLOCATION = ([4-23]=CELLCENTERED), ZONETYPE=%s\n" %
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
      for vindex in range(0,3):
        for i in range(0,ncell):
          f.write( str(vFieldGlobal[i][vindex]) + "    ")
          if ( i % 5  == 4 ):
            f.write("\n")
        f.write("\n") 
        
      #write heat flux X,Y,Z
      hfFieldGlobal = heatFluxFieldGlobal[n]
      for hfindex in range(0,3):
         for i in range(0,ncell):
            f.write( str(hfFieldGlobal[i][hfindex]) + "    ")
            if ( i % 5 == 4 ):
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
      
      #write Pxx,Pyy,Pzz,Pxy,Pyz,Pzx
      sFieldGlobal = stressFieldGlobal[n]
      for pindex in range(0,6):
        for i in range(0,ncell):
          f.write( str(sFieldGlobal[i][pindex]) + "    ")
          if ( i % 5  == 4 ):
            f.write("\n")
        f.write("\n")                  
     
      #write entropy
      eFieldGlobal = entropyFieldGlobal[n]
      for i in range(0,ncell):
         f.write( str(eFieldGlobal[i][0]) + "    ")
         if ( i % 5  == 4 ):
            f.write("\n")
      f.write("\n")

      #write entgenrate
      egFieldGlobal = entgenFieldGlobal[n]
      for i in range(0,ncell):
         f.write( str(egFieldGlobal[i][0]) + "    ")
         if ( i % 5  == 4 ):
            f.write("\n")
      f.write("\n")
      
      #write higher order moments Eq 32 in Gallis2006
      kFieldGlobal = KnqFieldGlobal[n]
      for i in range(0,ncell):
         f.write( str(kFieldGlobal[i][0]) + "    ")
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



