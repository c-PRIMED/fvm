#include "Grid.h"



typedef Vector<double,3>  VecD3;

typedef Array<VecD3> VecD3Array;


Grid::Grid():
  _grids(0),
  _gridCells(0),
  _coordinates(),
  _velocities()
{}

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

const shared_ptr<Array<VecD3> > Grid::readCoordinates(const char *file)

{
    FILE *fp;
    int nGrid;
    double x=0, y=0, z=0;
  
    fp=fopen(file,"r");
    fscanf(fp,"%i\n",&nGrid);
    
    shared_ptr<Array<VecD3> > Grid_Points ( new Array<VecD3> (nGrid));
    //read in coordinate
    for(int i=0; i<nGrid; i++){
      fscanf(fp,"%lf\t%lf\t%lf\n", &x, &y, &z);
      (*Grid_Points)[i][0]=x;
      (*Grid_Points)[i][1]=y;
      (*Grid_Points)[i][2]=z;
    }
    fclose(fp);   
    return (Grid_Points);
}

const shared_ptr<Array<VecD3> >  Grid::readVelocities(const char *file)

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

void Grid::setGridMesh(Mesh& mesh)
{
  const int nXGrid = 7;
  const int nYGrid = 3;

  //init storagesite for cells
  const int nGridCells = (nXGrid-1)*(nYGrid-1)*2;
  _gridCells.setCount(nGridCells);

  //create cell to grids connectivity
  //since it is tri cell, each cell connects to 3 grid points
  shared_ptr<CRConnectivity> cellToGrids (new CRConnectivity (_gridCells, _grids));
  (*cellToGrids).initCount();
  for (int n=0; n<nGridCells; n++){
    (*cellToGrids).addCount(n, 3);
  }
  (*cellToGrids).finishCount();
  int node1, node2, node3; 
  int id = 0;
  for (int n=0; n<nGridCells/2; n+=2){
    node1 = id;
    node2 = id+1;
    node3 = id+4;
    (*cellToGrids).add(n, node1);
    (*cellToGrids).add(n, node2);
    (*cellToGrids).add(n, node3);
    node1 = id;
    node2 = id+4;
    node3 = id+3;
    (*cellToGrids).add(n+1, node1);
    (*cellToGrids).add(n+1, node2);
    (*cellToGrids).add(n+1, node3);
    id+=3;
  }
  id = 1;
  for (int n=nGridCells/2; n<nGridCells;  n+=2){
    node1 = id;
    node2 = id+1;
    node3 = id+4;
    (*cellToGrids).add(n, node1);
    (*cellToGrids).add(n, node2);
    (*cellToGrids).add(n, node3);
    node1 = id;
    node2 = id+4;
    node3 = id+3;
    (*cellToGrids).add(n+1, node1);
    (*cellToGrids).add(n+1, node2);
    (*cellToGrids).add(n+1, node3);
    id+=3;
  }
  
  (*cellToGrids).finishAdd();

  //save cellToGrids connectivity;
  mesh.setConnectivity(_gridCells, _grids, cellToGrids);

  //test this connectivity
  const CRConnectivity& gridCellToGrids = mesh.getConnectivity(_gridCells, _grids);
  /*
  cout<<"cellToGrids connectivity is"<<endl;
  for (int n=0; n<nGridCells; n++){
    for (int c=0; c<3; c++){
      cout<<"cell "<<n<<" has nodes  "<<gridCellToGrids(n,c)<<endl;
    }
    }*/

}





void Grid::Init(const shared_ptr<Array<VecD3> > coordinates,
	       const shared_ptr<Array<VecD3> > velocities )
{

  const int n = (*coordinates).getLength();  //number of particles
  _grids.setCount(n);
  
  Grid::setCoordinates(coordinates);
  Grid::setVelocities(velocities);

}



//giving a point, find out which cell contains this point
//since the number of cells is not large, we search by loop over all cells
vector<int> Grid::findNeighborsByCells(const VecD3& point,
		    const Array<VecD3>& grids,
		    const int nCells, 
		    const CRConnectivity& cellToGrids)
{
  vector<int> neighborList;
  
  int find = 0;
  for ( int c=0; c<nCells; c++){
    const int nsize = cellToGrids.getCount(c);
    int node[nsize];
    for(int n=0; n<nsize; n++){
      node[n] = cellToGrids(c,n);
    }
    Array<VecD3> faceArea(nsize);
    Array<VecD3> faceCentroid(nsize);
    Array<VecD3> faceNorm(nsize);
    for(int n=0; n<nsize-1; n++){
      faceArea[n]=grids[node[n+1]]-grids[node[n]];
      faceCentroid[n]=(grids[node[n+1]]+grids[node[n]])/2.;
    }
    faceArea[nsize-1]=grids[node[0]]-grids[node[nsize-1]];
    faceCentroid[nsize-1]=(grids[node[0]]+grids[node[nsize-1]])/2.;
    for (int n=0; n<nsize; n++){
      faceNorm[n][0]=faceArea[n][1];
      faceNorm[n][1]=-faceArea[n][0];
      faceNorm[n][2]=faceArea[n][2];
    }

    int sum = 0;
    int flag[nsize];
    for(int n=0; n<nsize; n++){
      VecD3 dr = point - faceCentroid[n];
      double product = dot(faceNorm[n],dr);
      if(product > 0.0) flag[n] = 1;
      else flag[n] = -1;
      sum+=flag[n];
    }
    //point falls into a cell
    if (sum == -nsize){
      find=1;
      for(int n=0; n<nsize; n++){
	node[n]=cellToGrids(c,n);
	neighborList.push_back(node[n]);
      }
      return(neighborList);
    }
  }
  //point falls outside of grid boundary
  //find the cloest cell to this point
  //use that cell to interpolate
 
  
 if (find==0) {
    double distMin = 1.0e10;
    int cloestCell = 0; 
    for ( int c=0; c<nCells; c++){
      const int nsize = cellToGrids.getCount(c);
      VecD3 cellCentroid;
      cellCentroid[0]=0.0;
      cellCentroid[1]=0.0;
      cellCentroid[2]=0.0;
      for(int n=0; n<nsize; n++){
	cellCentroid += grids[cellToGrids(c,n)];
      }
      cellCentroid /= 3.;
      VecD3 dr = point - cellCentroid;
      double distance = mag(dr);
      if(distance < distMin){
	distMin = distance; 
	cloestCell = c;
      }
    }
    const int nsize = cellToGrids.getCount(cloestCell);
    for(int n=0; n<nsize; n++){
      neighborList.push_back(cellToGrids(cloestCell,n));
    }
    return(neighborList);
  }
      
}
    



//giving a point, find out its neighboring grids
//!!! only applicable for rectangular mesh
vector<int> Grid::findNeighbors(const VecD3& point, 
				const Array<VecD3>& grids,
				const int nGrid)
{
  vector<int> neighborList;

  //first, judge if this point falls into grids
  //if not, return NULL vector
  //if yes, continue to find the neighbors
  Vector<int,3> mFlag, pFlag;
  mFlag=pFlag=0;
  
  for (int i=0; i<nGrid; i++){
    for(int k=0; k<3; k++){
      if(point[k]<grids[i][k])	{
	mFlag[k]++;
      }
      if (point[k]>grids[i][k]) {
	pFlag[k]++;
      }
    }
  }
 

  if(mFlag[0]==nGrid||mFlag[1]==nGrid||mFlag[2]==nGrid||pFlag[0]==nGrid||pFlag[1]==nGrid||pFlag[2]==nGrid)
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
    for(int i=0; i<nGrid; i++){
      dR = point - grids[i];
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
  
    for(int i=0; i<nGrid; i++){
      if(grids[i][0]>(point[0]-pLimit[0]-allow[0])&& grids[i][0]<(point[0]-mLimit[0]+allow[0])){
	if(grids[i][1]>(point[1]-pLimit[1]-allow[1])&& grids[i][1]<(point[1]-mLimit[1]+allow[1])){
	  if(grids[i][2]>(point[2]-pLimit[2]-allow[2])&& grids[i][2]<(point[2]-mLimit[2]+allow[2])){
	    neighborList.push_back(i);
	  }
	}
      }
    }

   return(neighborList);
  }

}





//set up pointToGrid connectivity

const shared_ptr<CRConnectivity>  Grid::setConnectivity
                             (const StorageSite& pointSite, 
			      const StorageSite& gridSite,
			      const VecD3Array& points, 
			      const VecD3Array& grids,
			      Grid& grid,
			      const int nGridCells, 
			      const CRConnectivity& cellToGrids )
{

  shared_ptr<CRConnectivity> pointToGrids (new CRConnectivity (pointSite, gridSite));

  (*pointToGrids).initCount();

  const int nPoint = pointSite.getCount();
  const int nGrid = gridSite.getCount();
  
  for(int i=0; i<nPoint; i++){
   
    vector<int> neighborList = grid.findNeighborsByCells(points[i], grids,nGridCells, cellToGrids);
    int count = neighborList.size();
    if(count!=0)     
    (*pointToGrids).addCount(i, count);
  }
    
  (*pointToGrids).finishCount();

  for(int i=0; i<nPoint; i++){
    vector<int> neighborList = grid.findNeighborsByCells(points[i], grids,nGridCells, cellToGrids);
    int count = neighborList.size();
    for(int j=0; j<count; j++){
      (*pointToGrids).add(i, neighborList[j]);
    }
  }
  (*pointToGrids).finishAdd();

  return(pointToGrids);

}



void Grid::Impl(Mesh& mesh, GeomFields& geomFields, 
		FlowFields& flowFields, Grid& grid, 
		const string fileBase)
{
   const StorageSite& cells = mesh.getCells();

   const int nCells = cells.getSelfCount();         

   const VecD3Array& cellCentroid = 
     dynamic_cast<const VecD3Array& > (geomFields.coordinate[cells]);   

   const StorageSite& faces = mesh.getFaces();

   const int nFaces = faces.getSelfCount();        

   const VecD3Array& faceCentroid = 
     dynamic_cast<const VecD3Array& > (geomFields.coordinate[faces]);

  //set up grid data//  
  grid.setandwriteGrids(fileBase);

  //read grid coordinates
  string fileName1=fileBase+"Grid_Coord.dat";
  char *file1;
  file1=&fileName1[0];
  shared_ptr<VecD3Array> coordinates = grid.readCoordinates(file1);

  //read grid velocity
  string fileName2=fileBase+"Grid_Velocity.dat";
  char *file2;
  file2=&fileName2[0];
  shared_ptr<VecD3Array>  velocities = grid.readVelocities(file2);

  //store them in grid class
  grid.Init(coordinates,  velocities);

  const StorageSite& grids = grid.getGrids();

  const int nGrid = grids.getCount();

  cout<<"nGrid is "<<nGrid<<endl;

  const shared_ptr<VecD3Array> Grid_Coordinates = grid.getCoordinates();
  
  const shared_ptr<VecD3Array> Grid_Velocities = grid.getVelocities();

  //store coordinates in geomField
  geomFields.coordinate.addArray(grids, Grid_Coordinates);

  //store velocity in flowField
  flowFields.velocity.addArray(grids, Grid_Velocities);

  grid.setGridMesh(mesh);

  #if 0
  //test findCells using single point
  VecD3 point;
  point[0]=0.;
  point[1]=0.6;
  point[2]=0.0;
  const StorageSite& gridCells = grid.getCells();
  
  const int nGridCells = gridCells.getCount();

  const CRConnectivity& cellToGrids = mesh.getConnectivity(gridCells, grids);

  vector<int> nb = grid.findNeighborsByCells(point, (*Grid_Coordinates), nGridCells, cellToGrids);
  for(int n=0; n<nb.size();n++){
    cout<<"point is in cell "<<nb[n]<<endl;
  }
  #endif
#if 1

 //set up connectivity FaceToGrid
  const StorageSite& gridCells = grid.getCells();
  
  const int nGridCells = gridCells.getCount();

  const CRConnectivity& cellToGrids = mesh.getConnectivity(gridCells, grids);
 
  shared_ptr<CRConnectivity> faceGridsCR = 
    grid.setConnectivity(faces, grids, faceCentroid, *Grid_Coordinates, grid, nGridCells, cellToGrids);
 
  mesh.setConnectivity(faces, grids, faceGridsCR);

  const CRConnectivity& faceGrids = mesh.getConnectivity(faces, grids);

  #if 1
    string fileName=fileBase+"Face_Coord.dat";
    char *file;
    file=&fileName[0];
    FILE *fp;
    fp=fopen(file, "w");
    for(int n=0; n<nFaces; n++){
      if(faceGrids.getCount(n)!=0){
	fprintf(fp,"%i\t%f\t%f\t%f\n", n, faceCentroid[n][0],faceCentroid[n][1],
		faceCentroid[n][2]);
      }
    }
    fclose(fp);
    #endif
 

#endif
 

}

  
void Grid::computeInterpolatedVelocity(const StorageSite& grids, 
				   Grid& grid,
				   const Mesh& mesh,
				       const GeomFields& geomFields,
				       const string fileBase)
  {
    typedef CRMatrixTranspose<double,double,double> IMatrix;
    typedef CRMatrixTranspose<double,VecD3,VecD3> IMatrixV3;

    //const VectorT3Array& pV =
    //dynamic_cast<const VectorT3Array&>(_flowFields.velocity[grids]);
    const shared_ptr<Array<VecD3> >& gridV = grid.getVelocities();

    const StorageSite& faces = mesh.getFaces();

    const int nFaces = faces.getCount();

    const VecD3Array& faceCentroid = 
     dynamic_cast<const VecD3Array& > (geomFields.coordinate[faces]);

    GeomFields::SSPair key(&faces,&grids);
    const IMatrix& mIC =
      dynamic_cast<const IMatrix&> (*(geomFields._interpolationMatrices[key]));

    IMatrixV3 mICV(mIC);

    shared_ptr<VecD3Array> faceV(new VecD3Array(nFaces));

    faceV->zero();

    mICV.multiplyAndAdd(*faceV,*gridV);

    #if 1
    string fileName=fileBase+"Face_vel.dat";
    char *file;
    file=&fileName[0];
    FILE *fp;
    fp=fopen(file, "w");
    for(int n=0; n<nFaces; n++){
      fprintf(fp,"%i\t%f\t%f\t%f\t%f\t%f\t%f\n", n, faceCentroid[n][0],faceCentroid[n][1],
	      faceCentroid[n][2],(*faceV)[n][0],(*faceV)[n][1],(*faceV)[n][2]);
    }
    fclose(fp);
    #endif
  }
       




