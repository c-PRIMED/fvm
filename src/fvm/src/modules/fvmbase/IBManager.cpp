// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "IBManager.h"
#include "AABB.h"
#include "KSearchTree.h"
#include "Mesh.h"
#include "GradientModel.h"
#include <stack>
#include <iostream>
#include <fstream>

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif



IBManager::IBManager(GeomFields& geomFields,
                     Mesh& solidBoundaryMesh,
                     const MeshList& fluidMeshes):
  fluidNeighborsPerIBFace(50),
  fluidNeighborsPerSolidFace(50),
  solidNeighborsPerIBFace(50),
  IBNeighborsPerSolidFace(50),
  _geomFields(geomFields),
  _solidBoundaryMesh(solidBoundaryMesh),
  _fluidMeshes(fluidMeshes)
{}

void IBManager::update()
{
  AABB sMeshesAABB(_solidBoundaryMesh);
 
  const StorageSite& solidMeshFaces = _solidBoundaryMesh.getFaces();
  
  const Vec3DArray& solidMeshCoords =
    dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[solidMeshFaces]);

  KSearchTree solidMeshKSearchTree(solidMeshCoords);

  const int numFluidMeshes = _fluidMeshes.size();
  
  for (int n=0; n<numFluidMeshes; n++)
  {
      Mesh& fluidMesh = *_fluidMeshes[n];
      markIntersections(fluidMesh, sMeshesAABB);
  }


  int nIter=0;
  int nFound=0;

  // repeat till we find no more fluid cells
  do
  {
      nFound = 0;
      _geomFields.ibType.syncLocal();
      
      for (int n=0; n<numFluidMeshes; n++)
      {
          Mesh& fluidMesh = *_fluidMeshes[n];
          
          nFound += markFluid(fluidMesh);
      }
#ifdef FVM_PARALLEL
     MPI::COMM_WORLD.Allreduce( MPI::IN_PLACE, &nFound, 1, MPI::INT, MPI::SUM );
     if ( MPI::COMM_WORLD.Get_rank() == 0 ) 
        cout << "iteration " << nIter << ": found " << nFound << " fluid cells " << endl;
#endif

#ifndef FVM_PARALLEL
      cout << "iteration " << nIter << ": found " << nFound << " fluid cells " << endl;
#endif

      nIter++;
  } 
  while(nFound > 0);

  nFound=0;
  for (int n=0; n<numFluidMeshes; n++)
  {
      Mesh& fluidMesh = *_fluidMeshes[n];

      nFound += markSolid(fluidMesh);
  }
  _geomFields.ibType.syncLocal();  

  for (int n=0; n<numFluidMeshes; n++)
  {
      Mesh& fluidMesh = *_fluidMeshes[n];

      markIBTypePlus(fluidMesh);
  }

  for (int n=0; n<numFluidMeshes; n++)
  {
      Mesh& fluidMesh = *_fluidMeshes[n];

      createIBFaces(fluidMesh);
  }

 vector<NearestCell> solidFacesNearestCell(solidMeshFaces.getCount());
 // vector<NearestIBFace> solidFacesNearestIBFace(solidMeshFaces.getCount());
  
  for (int n=0; n<numFluidMeshes; n++) {
      Mesh& fluidMesh = *_fluidMeshes[n];

      if (!fluidMesh.isShell()){

	const StorageSite& cells = fluidMesh.getCells();
	const int numCells = cells.getSelfCount();
	//	StorageSite& ibFaces = fluidMesh.getIBFaces();
	//	const int nIBFaces = ibFaces.getCount();
	IntArray& cellIBType = dynamic_cast<IntArray&>(_geomFields.ibType[cells]);

	const Vec3DArray& cellCoords =
	  dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[cells]);
	//      	const Vec3DArray& ibFaceCoords =
	//	  dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[ibFaces]);
	KSearchTree fluidCellsTree;
	for(int c=0; c<numCells; c++)
	  {
	    if (cellIBType[c] == Mesh::IBTYPE_FLUID)
	      fluidCellsTree.insert(cellCoords[c],c);
	  }
	//	for(int c=0; c<nIBFaces; c++)
	//	  {
	//   	    IBFacesTree.insert(ibFaceCoords[c],c);
	//	  }
	createIBInterpolationStencil(fluidMesh,fluidCellsTree,solidMeshKSearchTree);
	
	findNearestCellForSolidFaces(fluidMesh,fluidCellsTree,solidFacesNearestCell);

      }
  }
  

#ifdef FVM_PARALLEL
    vector<doubleIntStruct>  solidFacesNearestCellMPI(solidMeshFaces.getCount());
    const int procID = MPI::COMM_WORLD.Get_rank();
    const int faceCount = solidMeshFaces.getCount();
    for( int i = 0; i < faceCount; i++ ){
       const Mesh* mesh = solidFacesNearestCell[i].mesh;
       if ( mesh != NULL){
          const int meshID = mesh->getID();
          const int tag = (std::max(procID, meshID) << 16 ) | ( std::min(procID,meshID) );
          const double val = solidFacesNearestCell[i].distanceSquared;
          solidFacesNearestCellMPI[i].VALUE = val;
          solidFacesNearestCellMPI[i].TAG   = tag;
      } else {
          const double LARGE_VALUE = 1.0e+15;
          solidFacesNearestCellMPI[i].VALUE = LARGE_VALUE;
          solidFacesNearestCellMPI[i].TAG   = -1;

      }
    } 
    //mpi comminuction
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &solidFacesNearestCellMPI[0], faceCount, MPI::DOUBLE_INT, MPI::MINLOC);
    //now update solidFAcesNearestCell
    for ( int i = 0; i < faceCount; i++ ){
       const Mesh* meshThis = solidFacesNearestCell[i].mesh;
       if ( meshThis != NULL ){
          const int meshIDThis = meshThis->getID();
          const int tagThis = (std::max(procID, meshIDThis) << 16 ) | ( std::min(procID,meshIDThis) );
          if ( tagThis != solidFacesNearestCellMPI[i].TAG ){
              solidFacesNearestCell[i].mesh = NULL;
          }

       } 
    }
#endif


  for (int n=0; n<numFluidMeshes; n++)
  {	
    Mesh& fluidMesh = *_fluidMeshes[n]; 
    if (!fluidMesh.isShell()){ 
      StorageSite& ibFaces = fluidMesh.getIBFaces();
      const StorageSite& faces = fluidMesh.getFaces();
      const Vec3DArray& faceCoords =
	  dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[faces]); 
      const int nIBFaces = ibFaces.getCount();
      const Array<int>& ibFaceIndices = fluidMesh.getIBFaceList();
      KSearchTree IBFacesTree;
      for(int c=0; c<nIBFaces; c++)
	{
	  const int gf=ibFaceIndices[c];
	  IBFacesTree.insert(faceCoords[gf],c);
	}
      //      findNearestIBFaceForSolidFaces(fluidMesh,IBFacesTree,solidFacesNearestIBFace);
      createSolidInterpolationStencil(fluidMesh,IBFacesTree,solidFacesNearestCell);
     
     }   
  }

}



void
IBManager::markIntersections(Mesh& fluidMesh, AABB& sMeshesAABB)
{
  
  const StorageSite& cells = fluidMesh.getCells();
  IntArray& cellIBType = dynamic_cast<IntArray&>(_geomFields.ibType[cells]);

  cellIBType = Mesh::IBTYPE_UNKNOWN;

  if (fluidMesh.isShell())
  {
      cellIBType = Mesh::IBTYPE_FLUID;
      return;
  }

  const Array<Vector<double,3> >& meshCoords = fluidMesh.getNodeCoordinates();
  
  const StorageSite& faces = fluidMesh.getFaces();
  const CRConnectivity& faceCells = fluidMesh.getAllFaceCells();
  const CRConnectivity& cellNodes = fluidMesh.getCellNodes();
 
  
  const int nFaces = faces.getCount();

  int nIntersections = 0;

  const bool is2D = fluidMesh.getDimension() == 2;

  if (is2D)
  {
      const int nCells = cells.getSelfCount();
      for(int n=0; n<nCells; n++)
      {
          const Vec3D& a = meshCoords[cellNodes(n,0)];
          const Vec3D& b = meshCoords[cellNodes(n,1)];
          const Vec3D& c = meshCoords[cellNodes(n,2)];
          
          if (sMeshesAABB.hasIntersectionWithTriangle(a,b,c))
          {
              cellIBType[n] = Mesh::IBTYPE_BOUNDARY;
              
              nIntersections++;
          }
          else if (cellNodes.getCount(n) == 4)
          {
              
              const Vec3D& d = meshCoords[cellNodes(n,3)];
              if (sMeshesAABB.hasIntersectionWithTriangle(c,d,a))
              {
                  cellIBType[n] = Mesh::IBTYPE_BOUNDARY;
                  nIntersections++;
              }
          }
      }
  }
  else
  {
      const CRConnectivity& faceNodes = fluidMesh.getAllFaceNodes();
      // loop through all the faces to find cells that intersect the AABB
      // faces and mark them to be of type boundary
      
      for(int f=0; f<nFaces; f++)
      {
          const int c0 = faceCells(f,0);
          const int c1 = faceCells(f,1);
          
          const Vec3D& a = meshCoords[faceNodes(f,0)];
          const Vec3D& b = meshCoords[faceNodes(f,1)];
          const Vec3D& c = meshCoords[faceNodes(f,2)];
          
          // check the triangle formed by the first three vertices
          if (sMeshesAABB.hasIntersectionWithTriangle(a,b,c))
          {
              cellIBType[c0] = Mesh::IBTYPE_BOUNDARY;
              cellIBType[c1] = Mesh::IBTYPE_BOUNDARY;
              
              nIntersections++;
          }
          else if (faceNodes.getCount(f) == 4)
          {
              // check the other triangle if this is a quad face
              const Vec3D& d = meshCoords[faceNodes(f,3)];
              if (sMeshesAABB.hasIntersectionWithTriangle(c,d,a))
              {
                  cellIBType[c0] = Mesh::IBTYPE_BOUNDARY;
                  cellIBType[c1] = Mesh::IBTYPE_BOUNDARY;
                  nIntersections++;
              }
          }
       }
  }


  // mark all cells adjacent to the boundaries as fluid unless they
  // have been found to be IBTYPE_BOUNDARY in the loop above
  
  foreach(const FaceGroupPtr fgPtr, fluidMesh.getBoundaryFaceGroups())
  {
      const FaceGroup& fg = *fgPtr;
      const StorageSite& faces = fg.site;
      
      const CRConnectivity& faceCells = fluidMesh.getFaceCells(faces);
      const int nFaces = faces.getCount();
      for(int f=0; f<nFaces; f++)
      {
          const int c0 = faceCells(f,0);
          const int c1 = faceCells(f,1);
          
          if ((cellIBType[c0] == Mesh::IBTYPE_UNKNOWN) &&
              (cellIBType[c1] == Mesh::IBTYPE_UNKNOWN))
          {
              cellIBType[c1] = Mesh::IBTYPE_FLUID;
              cellIBType[c0] = Mesh::IBTYPE_FLUID;
          }
          else if (cellIBType[c0] == Mesh::IBTYPE_BOUNDARY)
          {
              cellIBType[c1] = Mesh::IBTYPE_BOUNDARY;
          }
      }
  }
  
}


int
IBManager::markFluid(Mesh& fluidMesh)
{

  const StorageSite& cells = fluidMesh.getCells();

  IntArray& cellIBType = dynamic_cast<IntArray&>(_geomFields.ibType[cells]);


  const int nCellsTotal = cells.getCount();

  int nFound=0;

 
  const CRConnectivity& cellCells = fluidMesh.getCellCells();


  // loop over all cells including externals and if they are of type
  // fluid mark any other cells that are unknown but can be reached
  // from there
  
  for(int c=0; c<nCellsTotal; c++)
  {
      if (cellIBType[c] == Mesh::IBTYPE_FLUID)
      {
          
          stack<int> cellsToCheck;
          cellsToCheck.push(c);
          
           
           while(!cellsToCheck.empty())
          {
              int c_nb = cellsToCheck.top();
              cellsToCheck.pop();
              const int nNeighbors = cellCells.getCount(c_nb);
              for(int nn=0; nn<nNeighbors; nn++)
              {
                  const int neighbor = cellCells(c_nb,nn);
                  if (cellIBType[neighbor] == Mesh::IBTYPE_UNKNOWN)
                  {
                      cellIBType[neighbor] = Mesh::IBTYPE_FLUID;
                      nFound++;
                      cellsToCheck.push(neighbor);
                  }
              }
              
          }
      }
  }

  return nFound;
}

int
IBManager::markSolid(Mesh& fluidMesh)
{

  const StorageSite& cells = fluidMesh.getCells();

  IntArray& cellIBType = dynamic_cast<IntArray&>(_geomFields.ibType[cells]);


  const int nCellsTotal = cells.getSelfCount();

  int nFound=0;

  //everything that is not marked is solid  
  for(int c=0; c<nCellsTotal; c++)
  {
      if (cellIBType[c] == Mesh::IBTYPE_UNKNOWN)
      {
          cellIBType[c] = Mesh::IBTYPE_SOLID;
          nFound++;
      }
  }

  return nFound;
}


void
IBManager::markIBTypePlus(Mesh& fluidMesh)
{

#if 0
  int nFluid=0;
  int nSolid=0;
  int nBoundary=0;
  
  StorageSite& cells = fluidMesh.getCells();
  const int nCellsTotal = cells.getCountLevel1();
  const IntArray& cellIBType = dynamic_cast<IntArray&>(_geomFields.ibType[cells]);


  cout << " found " << nFluid << " fluid, "
       << nSolid << " solid and "
       << nBoundary << " boundary cells " << endl;
       
#endif        

#ifdef FVM_PARALLEL

#if 0
  int nCellsInner = cells.getSelfCount();
  int nFluidGlobal=0;
  int nSolidGlobal=0;
  int nBoundaryGlobal=0;
  vector<int> ibTypeCells;  
  const Array<int>& localToGlobal = fluidMesh.getLocalToGlobal();
  
 for(int c=0; c<nCellsInner; c++)
  {
      if (cellIBType[c] == Mesh::IBTYPE_FLUID)
        nFluid++;
      else if (cellIBType[c] == Mesh::IBTYPE_SOLID)
        nSolid++;
      else if (cellIBType[c] == Mesh::IBTYPE_BOUNDARY){
        nBoundary++;
	ibTypeCells.push_back(localToGlobal[c]);
      }
      else{
	throw CException("invalid ib type");
	}	
  }
  
     MPI::COMM_WORLD.Allreduce( &nFluid, &nFluidGlobal   , 1, MPI::INT, MPI::SUM);
     MPI::COMM_WORLD.Allreduce( &nSolid, &nSolidGlobal   , 1, MPI::INT, MPI::SUM);
     MPI::COMM_WORLD.Allreduce( &nBoundary, &nBoundaryGlobal, 1, MPI::INT, MPI::SUM);
     vector<int> ibTypeCellsGlobal(nBoundaryGlobal);
     const int nsize = MPI::COMM_WORLD.Get_size();
     const int rank  = MPI::COMM_WORLD.Get_rank();
     int counts [nsize];
     int offsets[nsize];
     MPI::COMM_WORLD.Allgather(&nBoundary, 1, MPI::INT, &counts[0], 1, MPI::INT);
     //form offsets
     offsets[0] = 0;
     for ( int p=1; p < nsize; p++){
        offsets[p] = counts[p-1] + offsets[p-1];
     }
     
     //gathering/
     MPI::COMM_WORLD.Gatherv(&ibTypeCells[0],counts[rank], MPI::INT,
        &ibTypeCellsGlobal[0], counts, offsets, MPI::INT, 0);


  if ( MPI::COMM_WORLD.Get_rank() == 0 ){
      set<int> ibTypeCellSet;
      for( int i = 0; i < int(ibTypeCellsGlobal.size()); i++ ){
         ibTypeCellSet.insert( ibTypeCellsGlobal[i] );
      }
      
      
      ofstream debugFile;
      debugFile.open("IBManagerDEBUG.dat");
      debugFile << " found (global) " << nFluidGlobal << " fluid, "   << nSolidGlobal << " solid and "
                << nBoundaryGlobal << " boundary cells " << endl;
      debugFile << "IBTYPE_BOUNDARY CELLS (GLOBAL NUMBERING)" << endl;
      foreach ( const set<int>::value_type cellID, ibTypeCellSet){
         debugFile << cellID << endl;
      }		

      debugFile.close(); 
       
  }      
#endif

#endif


}

void
IBManager::createIBFaces(Mesh& fluidMesh)
{

  if (fluidMesh.isShell())
    return;
  const StorageSite& faces = fluidMesh.getFaces();
  const StorageSite& cells = fluidMesh.getCells();
  
  const CRConnectivity& faceCells = fluidMesh.getAllFaceCells();

  IntArray& cellIBType = dynamic_cast<IntArray&>(_geomFields.ibType[cells]);
  IntArray& ibFaceIndex = dynamic_cast<IntArray&>(_geomFields.ibFaceIndex[faces]);

  const int nFaces = faces.getCount();

  ibFaceIndex = -1;
  
  // find number of IBFaces
  
  int nIBFaces=0;
  for(int f=0; f<nFaces; f++)
  {
      const int c0 = faceCells(f,0);
      const int c1 = faceCells(f,1);

      const int ibType0 = cellIBType[c0];
      const int ibType1 = cellIBType[c1];

      if ((ibType0 == Mesh::IBTYPE_FLUID && ibType1 == Mesh::IBTYPE_BOUNDARY) ||
          (ibType1 == Mesh::IBTYPE_FLUID && ibType0 == Mesh::IBTYPE_BOUNDARY))
      {
          ibFaceIndex[f]=nIBFaces++;
      }

      if ((ibType0 == Mesh::IBTYPE_FLUID && ibType1 == Mesh::IBTYPE_SOLID) ||
          (ibType1 == Mesh::IBTYPE_FLUID && ibType0 == Mesh::IBTYPE_SOLID))
      {
          throw CException("found face between solid and fluid cells");
      }

  }

  StorageSite& ibFaces = fluidMesh.getIBFaces();
  ibFaces.setCount(nIBFaces);

  shared_ptr<IntArray > ibFaceListPtr(new IntArray(nIBFaces));
  
  IntArray& ibFaceList = *ibFaceListPtr;

  nIBFaces = 0;
  for(int f=0; f<nFaces; f++)
  {
      const int c0 = faceCells(f,0);
      const int c1 = faceCells(f,1);

      const int ibType0 = cellIBType[c0];
      const int ibType1 = cellIBType[c1];

      if ((ibType0 == Mesh::IBTYPE_FLUID && ibType1 == Mesh::IBTYPE_BOUNDARY) ||
          (ibType1 == Mesh::IBTYPE_FLUID && ibType0 == Mesh::IBTYPE_BOUNDARY))
      {
          ibFaceList[nIBFaces++] = f;
      }
  }
  
  fluidMesh.setIBFaces(ibFaceListPtr);
  cout << " found " << nIBFaces << " ib Faces " << endl;

  GradientModelBase::clearGradientMatrix(fluidMesh);
}

/**
 * given a set of cells, add all their fluid ibtype neighbors to the
 * list if they aren't already in it.
 * 
 */

void addFluidNeighbors(set<int>& neighbors,
                       const CRConnectivity& cellCells,
                       const Array<int>& ibType)
{
  set<int> newNeighbors;
  foreach(int c, neighbors)
  {
      const int neighborCount = cellCells.getCount(c);
      for(int nnb=0; nnb<neighborCount; nnb++)
      {
          const int c_nb = cellCells(c,nnb);
          if (ibType[c_nb] == Mesh::IBTYPE_FLUID)
            newNeighbors.insert(c_nb);
      }
  }
  neighbors.insert(newNeighbors.begin(),newNeighbors.end());
}

void
IBManager::createIBInterpolationStencil(Mesh& mesh,
                                         KSearchTree& fluidCellsTree,
                                         KSearchTree& solidFacesTree)
{
  const StorageSite& cells = mesh.getCells();
  const StorageSite& solidMeshFaces = _solidBoundaryMesh.getFaces();
  StorageSite& ibFaces = mesh.getIBFaces();
  const int nIBFaces = ibFaces.getCount();

  if (nIBFaces == 0)
    return;
  
  const Vec3DArray& faceCentroid =
    dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[mesh.getFaces()]);

  Array<int> solidNeighbors(solidNeighborsPerIBFace);

  const Array<int>& ibFaceIndices = mesh.getIBFaceList();

  shared_ptr<CRConnectivity> ibFaceToCells
    (new CRConnectivity(ibFaces,cells));
  
  shared_ptr<CRConnectivity> ibFaceToSolid
    (new CRConnectivity(ibFaces,solidMeshFaces));


  //const CRConnectivity& cellCells  = mesh.getCellCells();
  const CRConnectivity& cellCells2 = mesh.getCellCells2();
  IntArray& cellIBType = dynamic_cast<IntArray&>(_geomFields.ibType[cells]);

  vector<NearestCell> nearestCellForIBFace(nIBFaces);

  const Vec3DArray& cellCoords =
    dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[cells]);

  const Vec3DArray& solidMeshCoords =
    dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[solidMeshFaces]);

  Array<int> desiredNeighbors(fluidNeighborsPerIBFace);

  ibFaceToCells->initCount();
  ibFaceToSolid->initCount();
  
  for(int f=0; f<nIBFaces; f++)
  {
      ibFaceToSolid->addCount(f,solidNeighborsPerIBFace);
  }
  
  ibFaceToSolid->finishCount();

  // in first pass we find the number of fluid cells and also set all the solid neighbors
  for(int f=0; f<nIBFaces; f++)
  {
      // the index of this ib face in the mesh faces
      const int gf = ibFaceIndices[f];
      const Vec3D& xf = faceCentroid[gf];

      // find the closest fluid cell
      Array<int> fluidNeighbors(1);
      fluidNeighbors[0] = -9999;
      fluidCellsTree.findNeighbors(xf, 1, fluidNeighbors);
      NearestCell& nc = nearestCellForIBFace[f];
      if ( fluidNeighbors[0] != -9999 ){
          nc.neighbors.push_back(fluidNeighbors[0]);
          const int c = fluidNeighbors[0];
          const int neighborCount = cellCells2.getCount(c);
          for(int nnb=0; nnb<neighborCount; nnb++) {
             const int c_nb = cellCells2(c,nnb);
             if (cellIBType[c_nb] == Mesh::IBTYPE_FLUID)
                nc.neighbors.push_back(c_nb);
          }
      }
      
      if ((int) nc.neighbors.size() > fluidNeighborsPerIBFace)
          {
	   
              // create a search tree to find the desired number of neighbors out of this set
              
              KSearchTree dtree;
              foreach(int nb, nc.neighbors)
              {
                  dtree.insert( cellCoords[nb], nb);
              }
	      bool swap = true;
	      while (swap == true){
		swap = false;
		for (int i=1; i<=(nc.neighbors.size()-1); i++){
		  const int nbLeft = nc.neighbors[i-1];
		  const int nbRight = nc.neighbors[i];
		  const double distanceLeft = mag2(cellCoords[nbLeft]-solidMeshCoords[f]);
		  const double distanceRight = mag2(cellCoords[nbRight]-solidMeshCoords[f]);
		  if (distanceLeft > distanceRight){
		    nc.neighbors[i-1] = nbRight;
		    nc.neighbors[i] = nbLeft;
		    swap = true;
		  }
		}
	      }
		
              //dtree.findNeighbors(solidMeshCoords[f], fluidNeighborsPerSolidFace, desiredNeighbors);
	      for (int i=0; i< fluidNeighborsPerIBFace; i++)
		desiredNeighbors[i] = nc.neighbors[i];
            
              // clear the current set of neighbors and add the desired ones
              nc.neighbors.clear();
              for(int i=0; i<fluidNeighborsPerIBFace; i++)
                nc.neighbors.push_back(desiredNeighbors[i]);
          }
          
      
      
#if 0
      int nLayers=0;
      // repeat till we have the required number but also protect
      // against infinite loop by capping the max number of layers
      while( ((int)nc.neighbors.size() < fluidNeighborsPerIBFace) && (nLayers < 10)) {
            addFluidNeighbors(nc.neighbors,cellCells,cellIBType);
            nLayers++;
      }
      
      if (nLayers == 10) 
        throw CException("not enough fluid cells for IB face interpolation");
#endif

      ibFaceToCells->addCount(f,nc.neighbors.size());

      // locate the required number of solid faces
      solidFacesTree.findNeighbors(xf, solidNeighborsPerIBFace, solidNeighbors);

      for(int n=0; n<solidNeighborsPerIBFace; n++)
        ibFaceToSolid->add(f,solidNeighbors[n]);
  }
 
  ibFaceToCells->finishCount();
  ibFaceToSolid->finishAdd();
 
  for(int f=0; f<nIBFaces; f++)
  {
      NearestCell& nc = nearestCellForIBFace[f];
      foreach(int nb, nc.neighbors)
      {
          ibFaceToCells->add(f,nb);
// 	  if (f==42){
// 	    cout << f << " " <<  nb << endl;
// 	    const int gf = ibFaceIndices[f];
// 	    cout << "face coordinate " << faceCentroid[gf] << endl;
// 	    cout << "cell coordinate " << cellCoords[nb] << endl;
// 	  }
      }
  }
  
  ibFaceToCells->finishAdd();
#ifdef FVM_PARALLEL
#if 0
  CRConnectivityPrintFile( *ibFaceToCells, "ibFaceToCells", MPI::COMM_WORLD.Get_rank() );
  CRConnectivityPrintFile( *ibFaceToSolid, "ibFaceToSolid", MPI::COMM_WORLD.Get_rank() );
#endif
#endif
  
  mesh.setConnectivity(ibFaces,cells,ibFaceToCells);
  mesh.setConnectivity(ibFaces,solidMeshFaces,ibFaceToSolid);
  
}

void
IBManager::findNearestCellForSolidFaces(Mesh& mesh,
                                        KSearchTree& fluidCellsTree,
                                        vector<NearestCell>& nearest)
                                           
{
  const StorageSite& cells = mesh.getCells();
  const StorageSite& solidMeshFaces = _solidBoundaryMesh.getFaces();

  const int nSolidFaces = solidMeshFaces.getCount();
  
  const Vec3DArray& cellCentroid =
    dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[cells]);

  const Vec3DArray& solidFaceCentroid =
    dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[solidMeshFaces]);


  for(int f=0; f<nSolidFaces; f++)
  {
      Array<int> fluidNeighbors(1);
      fluidNeighbors[0] = -9999;
      const Vec3D& xf = solidFaceCentroid[f];
      fluidCellsTree.findNeighbors(xf, 1, fluidNeighbors);
      if ( fluidNeighbors[0] != -9999 ){
         const int c = fluidNeighbors[0];
         const Vec3D& xc = cellCentroid[c];
         const double distanceSquared = mag2(xf-xc);
         NearestCell& nc = nearest[f];
         if ((nc.mesh == 0) || (nc.distanceSquared > distanceSquared))
         {
            nc.mesh = &mesh;
            nc.cell = c;
            nc.distanceSquared = distanceSquared;
         }
      }

   }

}



// void
// IBManager::findNearestIBFaceForSolidFaces(Mesh& mesh,
//                                         KSearchTree& IBFacesTree,
//                                         vector<NearestIBFace>& nearestIB)
                                           
// {
//   const StorageSite& cells = mesh.getCells();
//   const StorageSite& solidMeshFaces = _solidBoundaryMesh.getFaces();
//   StorageSite& ibFaces = mesh.getIBFaces();
//   const int nSolidFaces = solidMeshFaces.getCount();
//   const Vec3DArray& cellCentroid =
//     dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[cells]);

//   const Vec3DArray& solidFaceCentroid =
//     dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[solidMeshFaces]);
//   const Vec3DArray& IBFaceCentroid =
//     dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[ibFaces]);

//   for(int f=0; f<nSolidFaces; f++)
//   {
//       Array<int> IBNeighbors(1);
//       IBNeighbors[0] = -9999;
//       const Vec3D& xf = solidFaceCentroid[f];
//       IBFacesTree.findNeighbors(xf, 1, IBNeighbors);
//       if ( IBNeighbors[0] != -9999 ){
//          const int c = IBNeighbors[0];
//          const Vec3D& xc = IBFaceCentroid[c];
//          const double distanceSquared = mag2(xf-xc);
//          NearestIBFace& nc = nearestIB[f];
//          if ((nc.mesh == 0) || (nc.distanceSquared > distanceSquared))
//          {
//             nc.mesh = &mesh;
//             nc.IBFace = c;
//             nc.distanceSquared = distanceSquared;
//          }
//       }

//    }

// }


void
IBManager::createSolidInterpolationStencil(Mesh& mesh,
					   KSearchTree& IBFacesTree,
                                           vector<NearestCell>& nearest)
					   // vector<NearestIBFace>& nearestIB)
                                           
{
  const StorageSite& cells = mesh.getCells();
  const StorageSite& faces = mesh.getFaces();
  const int nFaces = faces.getCount();
  const StorageSite& solidMeshFaces = _solidBoundaryMesh.getFaces();
  const Vec3DArray& faceCoords =
    dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[faces]);
  const Vec3DArray& cellCoords =
    dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[cells]);
  StorageSite& ibFaces = mesh.getIBFaces();
  const Vec3DArray& solidMeshCoords =
    dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[solidMeshFaces]);

  const int nSolidFaces = solidMeshFaces.getCount();
  IntArray& cellIBType = dynamic_cast<IntArray&>(_geomFields.ibType[cells]);
  //const CRConnectivity& cellCells  = mesh.getCellCells();
  const CRConnectivity& cellCells2 = mesh.getCellCells2();
  

  shared_ptr<CRConnectivity> solidFacesToCells
    (new CRConnectivity(solidMeshFaces,cells));
  shared_ptr<CRConnectivity> solidToIBFaces
    (new CRConnectivity(solidMeshFaces,ibFaces));
  solidFacesToCells->initCount();
  solidToIBFaces->initCount();
  Array<int> desiredNeighbors(fluidNeighborsPerSolidFace);
  Array<int> desiredNeighborsforIBFaces(IBNeighborsPerSolidFace);
  for(int f=0; f<nSolidFaces; f++)
  {
      solidToIBFaces->addCount(f,IBNeighborsPerSolidFace);
  }

  solidToIBFaces->finishCount();
  
  for(int f=0; f<nSolidFaces; f++)
  {
      NearestCell& nc = nearest[f];
      if (nc.mesh == &mesh)
      {
          const int c = nc.cell;
          nc.neighbors.push_back(c);

          const int neighborCount = cellCells2.getCount(c);
          for(int nnb=0; nnb<neighborCount; nnb++) {
              const int c_nb = cellCells2(c,nnb);
              if (cellIBType[c_nb] == Mesh::IBTYPE_FLUID)
                 nc.neighbors.push_back(c_nb);
          }


          if ((int) nc.neighbors.size() > fluidNeighborsPerSolidFace)
          {
              // create a search tree to find the desired number of neighbors out of this set
              
              KSearchTree dtree;
              foreach(int nb, nc.neighbors)
              {
                  dtree.insert( cellCoords[nb], nb);
              }
	      // sort out the neighborlist by distance to ibface
	      bool swap = true;
	      while (swap == true){
		swap = false;
		for (int i=1; i<=(nc.neighbors.size()-1); i++){
		  const int nbLeft = nc.neighbors[i-1];
		  const int nbRight = nc.neighbors[i];
		  const double distanceLeft = mag2(cellCoords[nbLeft]-solidMeshCoords[f]);
		  const double distanceRight = mag2(cellCoords[nbRight]-solidMeshCoords[f]);
		  if (distanceLeft > distanceRight){
		    nc.neighbors[i-1] = nbRight;
		    nc.neighbors[i] = nbLeft;
		    swap = true;
		  }
		}
	      }
		

              //dtree.findNeighbors(solidMeshCoords[f], fluidNeighborsPerSolidFace, desiredNeighbors);
	      for (int i=0; i< fluidNeighborsPerSolidFace; i++)
		desiredNeighbors[i] = nc.neighbors[i];

              // clear the current set of neighbors and add the desired ones
              nc.neighbors.clear();
              for(int i=0; i<fluidNeighborsPerSolidFace; i++)
                nc.neighbors.push_back(desiredNeighbors[i]);
          }
          

#if 0
          int nLayers=0;
          // repeat till we have the required number but also protect 
          // against infinite loop by capping the max number of layers
          while( ((int)nc.neighbors.size() < fluidNeighborsPerIBFace) && (nLayers < 10)) {
             addFluidNeighbors(nc.neighbors,cellCells,cellIBType);
             nLayers++;
	  }

          //if (nLayers == 10)
          //  throw CException("not enough fluid cells for solid face interpolation");
#endif	  
          solidFacesToCells->addCount(f,nc.neighbors.size());
	
      }
          
  }
  
  solidFacesToCells->finishCount();


  for(int f=0; f<nSolidFaces; f++)
  {
      NearestCell& nc = nearest[f];
      if (nc.mesh == &mesh)
      {
          foreach(int nb, nc.neighbors)
          {
              solidFacesToCells->add(f,nb);
	      //if (f==100){
	      //   cout << f << " " <<  nb << endl;
		// cout << "face coordinate " << solidMeshCoords[f] << endl;
		// cout << "cell coordinate " << cellCoords[nb] << endl;
	      //}
          }
      }
  }
  
  solidFacesToCells->finishAdd();
  mesh.setConnectivity(solidMeshFaces,cells,solidFacesToCells);
 
 for(int f=0; f<nSolidFaces; f++)
  {
      // the index of this ib face in the mesh faces
    const Vec3D& xf = solidMeshCoords[f];
    const Array<int>& ibFaceIndices = mesh.getIBFaceList();
    IBFacesTree.findNeighbors(xf, IBNeighborsPerSolidFace, desiredNeighborsforIBFaces);
    for(int n=0; n< IBNeighborsPerSolidFace; n++){
      solidToIBFaces->add(f, desiredNeighborsforIBFaces[n]);
 //    if (f==20){
//        cout << f << " " <<  desiredNeighborsforIBFaces[n]<< endl;
//        cout << "face coordinate " << solidMeshCoords[f] << endl;
//        cout << "cell coordinate " << faceCoords[ibFaceIndices[desiredNeighborsforIBFaces[n]]] << endl;
      //   }
    }
  }

 solidToIBFaces->finishAdd();
 mesh.setConnectivity(solidMeshFaces,ibFaces,solidToIBFaces);

#ifdef FVM_PARALLEL
#if 0
  CRConnectivityPrintFile( *solidFacesToCells, "solidFacesToCells", MPI::COMM_WORLD.Get_rank() );
#endif
#endif 
}


void 
IBManager::CRConnectivityPrintFile(const CRConnectivity& conn, const string& name, const int procID) const
{
#ifdef FVM_PARALLEL
    if ( MPI::COMM_WORLD.Get_rank() == procID ){
       ofstream   debugFile;
       stringstream ss(stringstream::in | stringstream::out);
       ss <<  procID;
       string  fname = name +  ss.str() + ".dat";
       debugFile.open( fname.c_str() );
       ss.str("");

      debugFile <<  name << " :" << endl;
      debugFile << endl;
      const Array<int>& row = conn.getRow();
      const Array<int>& col = conn.getCol();
      for ( int i = 0; i < row.getLength()-1; i++ ){
         debugFile << " i = " << i << ",    ";
         for ( int j = row[i]; j < row[i+1]; j++ )
            debugFile << col[j] << "  ";
         debugFile << endl;
      }
      debugFile << endl;
   }
#endif
}
