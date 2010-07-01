#include "IBManager.h"
#include "AABB.h"
#include "KSearchTree.h"
#include "Mesh.h"
#include <stack>

IBManager::IBManager(GeomFields& geomFields,
                     Mesh& solidBoundaryMesh,
                     const MeshList& fluidMeshes):
  fluidNeighborsPerIBFace(2),
  solidNeighborsPerIBFace(4),
  fluidNeighborsPerSolidFace(4),
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

      markIB(fluidMesh, sMeshesAABB);

      const StorageSite& cells = fluidMesh.getCells();
      const int numCells = cells.getSelfCount();
      
      const Array<int>& cellIBType = fluidMesh.getOrCreateIBType();
      const Vec3DArray& cellCoords =
        dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[cells]);
      
      KSearchTree fluidCellsTree;

      for(int c=0; c<numCells; c++)
      {
          if (cellIBType[c] == Mesh::IBTYPE_FLUID)
            fluidCellsTree.insert(cellCoords[c],c);
      }

      createIBInterpolationStencil(fluidMesh,fluidCellsTree,solidMeshKSearchTree);
      createSolidInterpolationStencil(fluidMesh,fluidCellsTree);
                                              
  }
}



void
IBManager::markIB(Mesh& fluidMesh, AABB& sMeshesAABB)
{
  const Array<Vector<double,3> >& meshCoords = fluidMesh.getNodeCoordinates();
  const StorageSite& faces = fluidMesh.getFaces();
  const StorageSite& cells = fluidMesh.getCells();
  
  const CRConnectivity& faceCells = fluidMesh.getAllFaceCells();
  const CRConnectivity& cellNodes = fluidMesh.getCellNodes();

  Array<int>& cellIBType = fluidMesh.getOrCreateIBType();

  cellIBType = Mesh::IBTYPE_UNKNOWN;
  
  const int nFaces = faces.getCount();
  const int nCells = cells.getCount();
  
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

  // find one cell that is definitely inside or outside by brute force search

  int c=0;
  int iBType=Mesh::IBTYPE_UNKNOWN;
  
  for(c=0; c<nCells; c++)
  {
      if (cellIBType[c] == Mesh::IBTYPE_UNKNOWN)
      {
          const int nNodes = cellNodes.getCount(c);

          // the side the first node is on. if all the other nodes of
          // the cell are on the same side then we know the cell is on
          // the same side. if any node is on the AABB faces we skip
          // checking it
          
          int s0 = sMeshesAABB.findOrientedSide(meshCoords[cellNodes(c,0)]);

          bool located = true;
          for (int nn=1; nn<nNodes; nn++)
          {
              int s1 = sMeshesAABB.findOrientedSide(meshCoords[cellNodes(c,nn)]);
              if (s1 !=0)
              {
                  // handle the possibility that the first node was on the AABB faces
                  if (s0==0)
                    s0 = s1;
                  else if (s0 != s1)
                  {
                      located = false;
                      break;
                  }
              }
          }
          if (located)
          {
              if (s0 == 1)
                iBType = Mesh::IBTYPE_FLUID;
              else if (s0 == -1)
                iBType = Mesh::IBTYPE_SOLID;
              else
                throw CException("found cell with all nodes on solid boundary");
          }
          else
            throw CException("found cell with intersections with solid boundary");

          // now that we know the ib type of cell c we can set the type to be
          // the same for all the unmarked cells we can recursively reach
          // through the neighbours list.
          
          cellIBType[c] = iBType;
          
          stack<int> cellsToCheck;
          cellsToCheck.push(c);
          
          const CRConnectivity& cellCells = fluidMesh.getCellCells();
          
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
                      cellIBType[neighbor] = iBType;
                      cellsToCheck.push(neighbor);
                  }
              }
              
          }
      }
  }
  
  int nFluid=0;
  int nSolid=0;
  int nBoundary=0;

  for(int c=0; c<nCells; c++)
  {
      if (cellIBType[c] == Mesh::IBTYPE_FLUID)
        nFluid++;
      else if (cellIBType[c] == Mesh::IBTYPE_SOLID)
        nSolid++;
      else if (cellIBType[c] == Mesh::IBTYPE_BOUNDARY)
        nBoundary++;
      else
        throw CException("invalid ib type");
  }

  cout << " found " << nFluid << " fluid, "
       << nSolid << " solid and "
       << nBoundary << " boundary cells " << endl;

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
          nIBFaces++;
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

  const Vec3DArray& faceCentroid =
    dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[mesh.getFaces()]);

  Array<int> fluidNeighbors(fluidNeighborsPerIBFace);
  Array<int> solidNeighbors(solidNeighborsPerIBFace);

  const Array<int>& ibFaceIndices = mesh.getIBFaceList();

  shared_ptr<CRConnectivity> ibFaceToCells
    (new CRConnectivity(ibFaces,cells));
  
  shared_ptr<CRConnectivity> ibFaceToSolid
    (new CRConnectivity(ibFaces,solidMeshFaces));
  
  ibFaceToCells->initCount();
  ibFaceToSolid->initCount();

  for(int f=0; f<nIBFaces; f++)
  {
      ibFaceToCells->addCount(f,fluidNeighborsPerIBFace);
      ibFaceToSolid->addCount(f,solidNeighborsPerIBFace);
  }

  ibFaceToCells->finishCount();
  ibFaceToSolid->finishCount();

  for(int f=0; f<nIBFaces; f++)
  {
      // the index of this ib face in the mesh faces
      const int gf = ibFaceIndices[f];
      const Vec3D& xf = faceCentroid[gf];
      fluidCellsTree.findNeighbors(xf, fluidNeighborsPerIBFace, fluidNeighbors);
      solidFacesTree.findNeighbors(xf, solidNeighborsPerIBFace, solidNeighbors);

      for(int n=0; n<fluidNeighborsPerIBFace; n++)
        ibFaceToCells->add(f,fluidNeighbors[n]);
      for(int n=0; n<solidNeighborsPerIBFace; n++)
        ibFaceToSolid->add(f,solidNeighbors[n]);
  }

  ibFaceToCells->finishAdd();
  ibFaceToSolid->finishAdd();

  mesh.setConnectivity(ibFaces,cells,ibFaceToCells);
  mesh.setConnectivity(ibFaces,solidMeshFaces,ibFaceToSolid);
  
}


void
IBManager::createSolidInterpolationStencil(Mesh& mesh,
                                           KSearchTree& fluidCellsTree)
                                           
{
  const StorageSite& cells = mesh.getCells();
  const StorageSite& solidMeshFaces = _solidBoundaryMesh.getFaces();

  const int nSolidFaces = solidMeshFaces.getCount();
  
  const Vec3DArray& cellCentroid =
    dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[cells]);

  const Vec3DArray& solidFaceCentroid =
    dynamic_cast<const Vec3DArray&>(_geomFields.coordinate[solidMeshFaces]);


  Array<int> fluidNeighbors(fluidNeighborsPerSolidFace);

  AABB meshAABB(mesh);

  shared_ptr<CRConnectivity> solidFacesToCells
    (new CRConnectivity(solidMeshFaces,cells));
  
  solidFacesToCells->initCount();

  Array<bool> inThisMesh(nSolidFaces);
  
  for(int f=0; f<nSolidFaces; f++)
  {
      const Vec3D& xf = solidFaceCentroid[f];
      fluidCellsTree.findNeighbors(xf, fluidNeighborsPerSolidFace,
                                   fluidNeighbors);

      // check if the line joining face to nearest cell intersects the
      // mesh boundaries
      
      const Vec3D& xc = cellCentroid[fluidNeighbors[0]];
      if (!meshAABB.hasIntersectionWithSegment(xf,xc))
      {
          solidFacesToCells->addCount(f,fluidNeighborsPerSolidFace);
          inThisMesh[f] = true;
      }
      else
        inThisMesh[f] = false;
  }

  solidFacesToCells->finishCount();

  for(int f=0; f<nSolidFaces; f++)
  {
      if (inThisMesh[f])
      {
          const Vec3D& xf = solidFaceCentroid[f];
          fluidCellsTree.findNeighbors(xf, fluidNeighborsPerSolidFace,
                                       fluidNeighbors);
          for(int n=0; n<fluidNeighborsPerSolidFace; n++)
            solidFacesToCells->add(f,fluidNeighbors[n]);
      }
  }
  
  solidFacesToCells->finishAdd();
  mesh.setConnectivity(solidMeshFaces,cells,solidFacesToCells);
}
