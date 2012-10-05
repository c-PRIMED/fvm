// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "Grid.h"

typedef Vector<double,3>  VecD3;

typedef Array<VecD3> VecD3Array;

shared_ptr<Array<VecD3> >
readVectors(const char *file)
{
  FILE *fp;
  int nGrid;
  double vx=0, vy=0, vz=0;
   
  fp=fopen(file,"r");
  fscanf(fp,"%i\n",&nGrid);
  shared_ptr<Array<VecD3> > Grid_Points  (new Array<VecD3> (nGrid));
   
  //read in velocity
  for(int i=0; i<nGrid; i++){
      fscanf(fp,"%lf\t%lf\t%lf\n", &vx, &vy, &vz);
      (*Grid_Points)[i][0]=vx;
      (*Grid_Points)[i][1]=vy;
      (*Grid_Points)[i][2]=vz;
  }
  fclose(fp);   
  return (Grid_Points);
}

Grid::Grid(GeomFields& geomFields, FlowFields& flowFields, string coordFileName,
           string velocityFileName):
  _geomFields(geomFields),
  _flowFields(flowFields),
  _nodes(0),
  _cells(0),
  _coordinates(),
  _velocities()
{

    //set up grid data//  
  //grid.setandwriteGrids(fileBase);

  //read grid coordinates

  _coordinates = readVectors(coordFileName.c_str());

  //read grid velocity
  _velocities = readVectors(velocityFileName.c_str());

  _nodes.setCount(_velocities->getLength());

  //store coordinates in geomField
  geomFields.coordinate.addArray(_nodes, _coordinates);

  //store velocity in flowField
  flowFields.velocity.addArray(_nodes, _velocities);

  // create cell to node connectivity
  createCellToNodeConnectivity();

#if 0
  //test findCells using single point
  VecD3 point;
  point[0]=0.;
  point[1]=0.6;
  point[2]=0.0;
  const StorageSite& gridCells = grid.getCells();
  
  const int nGridCells = gridCells.getCount();

  const CRConnectivity& cellToGrids = mesh.getConnectivity(gridCells, grids);

  vector<int> nb = findNeighborsByCells(point, (*Grid_Coordinates), nGridCells, cellToGrids);
  for(int n=0; n<nb.size();n++){
      cout<<"point is in cell "<<nb[n]<<endl;
  }
#endif
}

Grid::~Grid() { }

void Grid::setandwriteGrids(const string fileBase)
{
 
  int nX=7, nY=3, nZ=1;
  double gapX=1.4/nX, gapY=1.2/nY, gapZ=0.0/nZ;
    
  VecD3 center;
  center[0]=0.501;
  center[1]=0.501;
  center[2]=0.0;

  int count=0;
  VecD3 gridCoord[nX*nY*nZ];
  VecD3 gridVelocity[nX*nY*nZ];

  //set up grid coordinate
  for(int i=0; i<nX; i++){
      for(int j=0; j<nY; j++){
          for(int k=0; k<nZ; k++){
              gridCoord[count][0]=i*gapX-0.023;//+(center[0]-(nX*gapX/2.));
              gridCoord[count][1]=-j*gapY+0.92;//+(center[1]-(nY*gapY/2.));
              gridCoord[count][2]=k*gapZ;
              count+=1;
          }
      }
  }

  //set up grid velocity
  const double vMax=1.0;
  const double vMin=0.0;
  const double vGap=(vMax-vMin)/nX;
  count=0;
  for(int i=0; i<nX; i++){
      for(int j=0; j<nY; j++){
          for(int k=0; k<nZ; k++){
              gridVelocity[count][0]=vMin+vGap*i;
              gridVelocity[count][1]=0.0;
              gridVelocity[count][2]=0.0;
              count++;
          }
      }
  }
  //write out coordinate into file
  string fileName1=fileBase+"Grid_Coord.dat";
  char* file1;
  file1=&fileName1[0];
  FILE *fp1=fopen(file1,"w");
  fprintf(fp1,"%i\n",count);
  for(int p=0; p<count; p++){
      fprintf(fp1, "%lf\t%lf\t%lf\n", gridCoord[p][0], gridCoord[p][1], gridCoord[p][2]);
  } 
  fclose(fp1);

  //write out velocity into file
  string fileName2=fileBase+"Grid_Velocity.dat";
  char* file2;
  file2=&fileName2[0];
  FILE *fp2=fopen(file2,"w");
  fprintf(fp2,"%i\n",count);
  for(int p=0; p<count; p++){
      fprintf(fp2, "%lf\t%lf\t%lf\n", gridVelocity[p][0],gridVelocity[p][1],gridVelocity[p][2]);
  }     
  fclose(fp2);
}


void Grid::createCellToNodeConnectivity()
{
  const int nXGrid = 7;
  const int nYGrid = 3;

  //init storagesite for cells
  const int nGridCells = (nXGrid-1)*(nYGrid-1)*2;
  _cells.setCount(nGridCells);

  //create cell to grids connectivity
  //since it is tri cell, each cell connects to 3 grid points
  _cellNodes = shared_ptr<CRConnectivity>(new CRConnectivity (_cells, _nodes));
  (*_cellNodes).initCount();
  for (int n=0; n<nGridCells; n++){
      (*_cellNodes).addCount(n, 3);
  }
  (*_cellNodes).finishCount();
  int node1, node2, node3; 
  int id = 0;
  for (int n=0; n<nGridCells/2; n+=2){
      node1 = id;
      node2 = id+1;
      node3 = id+4;
      (*_cellNodes).add(n, node1);
      (*_cellNodes).add(n, node2);
      (*_cellNodes).add(n, node3);
      node1 = id;
      node2 = id+4;
      node3 = id+3;
      (*_cellNodes).add(n+1, node1);
      (*_cellNodes).add(n+1, node2);
      (*_cellNodes).add(n+1, node3);
      id+=3;
  }
  id = 1;
  for (int n=nGridCells/2; n<nGridCells;  n+=2){
      node1 = id;
      node2 = id+1;
      node3 = id+4;
      (*_cellNodes).add(n, node1);
      (*_cellNodes).add(n, node2);
      (*_cellNodes).add(n, node3);
      node1 = id;
      node2 = id+4;
      node3 = id+3;
      (*_cellNodes).add(n+1, node1);
      (*_cellNodes).add(n+1, node2);
      (*_cellNodes).add(n+1, node3);
      id+=3;
  }
  
  (*_cellNodes).finishAdd();

}

//giving a point, find out which cell contains this point
//since the number of cells is not large, we search by loop over all cells
vector<int>
Grid::findNeighborsByCells(const VecD3& point)
{
  vector<int> neighborList;
  const int nCells = _cells.getCount();
  const CRConnectivity& cellNodes = *_cellNodes;
  const Array<VecD3>& nodeCoords = *_coordinates;
  
  for ( int c=0; c<nCells; c++){
    const int nsize = cellNodes.getCount(c);
    int node[nsize];
    for(int n=0; n<nsize; n++){
      node[n] = cellNodes(c,n);
    }
    Array<VecD3> faceArea(nsize);
    Array<VecD3> faceCentroid(nsize);
    for(int n=0; n<nsize-1; n++){
      faceArea[n]=nodeCoords[node[n+1]]-nodeCoords[node[n]];
      faceCentroid[n]=(nodeCoords[node[n+1]]+nodeCoords[node[n]])/2.;
    }
    faceArea[nsize-1]=nodeCoords[node[0]]-nodeCoords[node[nsize-1]];
    faceCentroid[nsize-1]=(nodeCoords[node[0]]+nodeCoords[node[nsize-1]])/2.;
    
    int find = 1;
    for(int n=0; n<nsize; n++){
      VecD3 faceNorm;
      faceNorm[0]=faceArea[n][1];
      faceNorm[1]=-faceArea[n][0];
      faceNorm[2]=faceArea[n][2];

      VecD3 dr = point - faceCentroid[n];
      dr[2] = 0.0;  //difference in z direction is neglected. just check x-y plane
      if (dot(faceNorm,dr) > 0.0) {
	find = 0;
	break;
      }
    }
    //point falls into a cell
    if (find){
      for(int n=0; n<nsize; n++){
	node[n]=cellNodes(c,n);
	neighborList.push_back(node[n]);
      }
      return(neighborList);
    }
  }

  //point falls outside of grid boundary
  //find the cloest cell to this point
  //use that cell to interpolate

  double distMin = 1.0e10;
  int closestCell = 0; 
  for ( int c=0; c<nCells; c++){
    const int nsize = cellNodes.getCount(c);
    VecD3 cellCentroid;
    cellCentroid[0]=0.0;
    cellCentroid[1]=0.0;
    cellCentroid[2]=0.0;
    for(int n=0; n<nsize; n++){
      cellCentroid += nodeCoords[cellNodes(c,n)];
    }
    cellCentroid /= 3.;
    VecD3 dr = point - cellCentroid;
    double distance = mag(dr);
    if(distance < distMin){
      distMin = distance; 
      closestCell = c;
    }
  }
  const int nsize = cellNodes.getCount(closestCell);
  for(int n=0; n<nsize; n++){
    neighborList.push_back(cellNodes(closestCell,n));
  }
  return(neighborList);
}
    



//giving a point, find out its neighboring nodes
//!!! only applicable for rectangular mesh
vector<int>
Grid::findNeighbors(const VecD3& point)

{
  const Array<VecD3>& nodeCoords = *_coordinates;
  const int nNodes = nodeCoords.getLength();
  
  vector<int> neighborList;

  //first, judge if this point falls into nodeCoords
  //if not, return NULL vector
  //if yes, continue to find the neighbors
  Vector<int,3> mFlag, pFlag;
  mFlag=pFlag=0;
  
  for (int i=0; i<nNodes; i++){
      for(int k=0; k<3; k++){
          if(point[k]<nodeCoords[i][k])	{
              mFlag[k]++;
          }
          if (point[k]>nodeCoords[i][k]) {
              pFlag[k]++;
          }
      }
  }
 

  if(mFlag[0]==nNodes||mFlag[1]==nNodes||mFlag[2]==nNodes||pFlag[0]==nNodes||pFlag[1]==nNodes||pFlag[2]==nNodes)
    return (neighborList);

  else{	 
   
      VecD3 pLimit, mLimit, allow, dR;
      for(int k=0; k<3; k++){
          pLimit[k]=1e10;
          mLimit[k]=-1e10;
          allow[k]=1.e-8;
          dR[k]=0.0;
      }
      //find the bounds for point
      for(int i=0; i<nNodes; i++){
          dR = point - nodeCoords[i];
          for(int k=0; k<3; k++){
              if( dR[k]<=0.0 && dR[k]>mLimit[k] ){
                  mLimit[k]=dR[k];
              }
              if (dR[k]>=0.0 && dR[k]<pLimit[k]){
                  pLimit[k]=dR[k];
              }
          }
      }
    
      //find who bounds the point
  
      for(int i=0; i<nNodes; i++){
          if(nodeCoords[i][0]>(point[0]-pLimit[0]-allow[0])&& nodeCoords[i][0]<(point[0]-mLimit[0]+allow[0])){
              if(nodeCoords[i][1]>(point[1]-pLimit[1]-allow[1])&& nodeCoords[i][1]<(point[1]-mLimit[1]+allow[1])){
                  if(nodeCoords[i][2]>(point[2]-pLimit[2]-allow[2])&& nodeCoords[i][2]<(point[2]-mLimit[2]+allow[2])){
                      neighborList.push_back(i);
                  }
              }
          }
      }

      return(neighborList);
  }

}


//set up pointToNodes connectivity

const shared_ptr<CRConnectivity>
Grid::createConnectivity(const StorageSite& pointSite, 
                         const VecD3Array& points)
{
  shared_ptr<CRConnectivity> pointToNodes (new CRConnectivity (pointSite, _nodes));

  (*pointToNodes).initCount();

  const int nPoints = pointSite.getCount();
  
  for(int i=0; i<nPoints; i++){
   
      vector<int> neighborList = findNeighborsByCells(points[i]);
      int count = neighborList.size();
      if(count!=0)     
        (*pointToNodes).addCount(i, count);
  }
    
  (*pointToNodes).finishCount();

  for(int i=0; i<nPoints; i++){
      vector<int> neighborList = findNeighborsByCells(points[i]);
      int count = neighborList.size();
      for(int j=0; j<count; j++){
          (*pointToNodes).add(i, neighborList[j]);
      }
  }
  (*pointToNodes).finishAdd();

  return(pointToNodes);

}



void Grid::setConnFaceToGrid(Mesh& mesh,
			     const StorageSite& faces)
{
  const VecD3Array& faceCentroid = 
    dynamic_cast<const VecD3Array& > (_geomFields.coordinate[faces]);
  
  shared_ptr<CRConnectivity> faceGridsCR = 
    createConnectivity(faces, faceCentroid);
 
  mesh.setConnectivity(faces, _nodes, faceGridsCR);
}

  
shared_ptr<ArrayBase>
Grid::computeInterpolatedVelocity(const StorageSite& faces)
{
  typedef CRMatrixTranspose<double,double,double> IMatrix;
  typedef CRMatrixTranspose<double,VecD3,VecD3> IMatrixV3;


  const int nFaces = faces.getCount();
    
    
  GeomFields::SSPair key(&faces,&_nodes);
  const IMatrix& mIC =
    dynamic_cast<const IMatrix&> (*(_geomFields._interpolationMatrices[key]));

  IMatrixV3 mICV(mIC);
	  
  shared_ptr<VecD3Array> faceV(new VecD3Array(nFaces));
	  
  faceV->zero();

  mICV.multiplyAndAdd(*faceV,*_velocities);

  return (faceV);
      
}



 
