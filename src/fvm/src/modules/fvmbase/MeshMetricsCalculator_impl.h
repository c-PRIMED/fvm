#include "Model.h"

#include "NumType.h"
#include "Array.h"
#include "Vector.h"
#include "Field.h"
#include "CRConnectivity.h"
#include "StorageSite.h"
#include "GlobalFields.h"
#include "CRMatrixTranspose.h"
#include "Mesh.h"
#include "MatrixOperation.h"

  
/**
 * The mesh always has coordinates in double. The coordField that
 * the rest of the code uses to get the coordinates needs to have
 * them in the current AType so this function just creates a copy in
 * that AType.
 * 
 */

template<class T>
void
MeshMetricsCalculator<T>::calculateNodeCoordinates(const Mesh& mesh)
{
  const StorageSite& nodes = mesh.getNodes();
  const Array<Vector<double,3> >& meshCoords = mesh.getNodeCoordinates();

  const int count = nodes.getCount();
  shared_ptr<VectorT3Array> ncPtr(new VectorT3Array(count));
  VectorT3Array& nc = *ncPtr;

  for(int i=0; i<count; i++)
    for(int j=0; j<3; j++)
      nc[i][j] = T(meshCoords[i][j]);

  _coordField.addArray(nodes,ncPtr);
}

/**
 * calculates the face centroids. Needs face area and magnitude for
 * non-planar corrections.
 * 
 */

template<class T>
void
MeshMetricsCalculator<T>::calculateFaceCentroids(const Mesh& mesh)
{
  const StorageSite& faces = mesh.getFaces();
  const StorageSite& nodes = mesh.getNodes();
  const CRConnectivity& faceNodes = mesh.getAllFaceNodes();

  const int count = faces.getCount();
  shared_ptr<VectorT3Array> fcPtr(new VectorT3Array(count));
  VectorT3Array& fc = *fcPtr;

  const VectorT3Array& nodeCoord =
    dynamic_cast<const VectorT3Array&>(_coordField[nodes]);
  for(int f=0; f<count; f++)
  {
      const int numNodes = faceNodes.getCount(f);
      if (numNodes > 0)
      {
          fc[f] = nodeCoord[faceNodes(f,0)];
          for (int nf=1; nf<numNodes; nf++)
            fc[f] += nodeCoord[faceNodes(f,nf)];
          fc[f] /= T(numNodes);
      }
  }

  const VectorT3Array& faceArea =
    dynamic_cast<const VectorT3Array&>(_areaField[faces]);
  const TArray& faceAreaMag =
    dynamic_cast<const TArray&>(_areaMagField[faces]);
      
  // corrections for non-planar quad and polygonal faces
  const T twoThirds(2./3.);
  const T half(0.5);
  for(int f=0; f<count; f++)
  {
      const int numNodes = faceNodes.getCount(f);
      if (numNodes > 3)
      {
          const VectorT3 en = faceArea[f]/faceAreaMag[f];
          T denom(0.0);
          VectorT3 cfc(NumTypeTraits<VectorT3>::getZero());

          for (int nn=0; nn<numNodes; nn++)
          {
              const int n0 = faceNodes(f,nn);
              const int n1 = faceNodes(f,(nn+1)%numNodes);
              const VectorT3 rc0 = nodeCoord[n0] - fc[f];
              const VectorT3 rc1 = nodeCoord[n1] - fc[f];
              const VectorT3 triArea = half*cross(rc0,rc1);
              const T triAreaP = dot(triArea,en); 
              VectorT3 xm = half*(nodeCoord[n0]+nodeCoord[n1]);

              cfc += twoThirds*(xm-fc[f])*triAreaP;
              denom += triAreaP;
          }
          cfc /= denom;
          fc[f] += cfc;
      }
  }

  _coordField.addArray(faces,fcPtr);
}

/**
 * calculate cell centroids. Needs to have the face centroid and
 * areas computed already.
 * 
 */

template<class T>
void
MeshMetricsCalculator<T>::calculateCellCentroids(const Mesh &mesh)
{
  const StorageSite& cells = mesh.getCells();
      
  const int cellCount = cells.getCount();
  const int selfCellCount = cells.getSelfCount();

  shared_ptr<VectorT3Array> ccPtr(new VectorT3Array(cellCount));
  VectorT3Array& cellCentroid =  *ccPtr;

  const StorageSite& faces = mesh.getFaces();

  const VectorT3Array& faceCentroid =
    dynamic_cast<const VectorT3Array&>(_coordField[faces]);
  const TArray& faceAreaMag = dynamic_cast<const TArray&>(_areaMagField[faces]);
  const int faceCount = faces.getCount();
  const CRConnectivity& faceCells = mesh.getAllFaceCells();

  TArray weight(cellCount);

  weight.zero();
  cellCentroid.zero();
      
  for(int f=0; f<faceCount; f++)
  {
      for(int nc=0; nc<faceCells.getCount(f); nc++)
      {
          const int c = faceCells(f,nc);
          cellCentroid[c] += faceCentroid[f]*faceAreaMag[f];
          weight[c] += faceAreaMag[f];
      }
  }

  for(int c=0; c<selfCellCount; c++)
    cellCentroid[c] /= weight[c];

  // boundary cells have the corresponding face's centroid
  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
  {
      const FaceGroup& fg = *fgPtr;
      const StorageSite& faces = fg.site;
      const VectorT3Array& faceCentroid =
        dynamic_cast<const VectorT3Array&>(_coordField[faces]);
      const CRConnectivity& faceCells = mesh.getFaceCells(faces);
      const int faceCount = faces.getCount();
          
      for(int f=0; f<faceCount; f++)
      {
          const int c1 = faceCells(f,1);
          cellCentroid[c1] = faceCentroid[f];
      }
  }
  _coordField.addArray(cells,ccPtr);
}
    
template<class T>
void
MeshMetricsCalculator<T>::calculateFaceAreas(const Mesh& mesh)
{
  const StorageSite& faces = mesh.getFaces();
  const StorageSite& nodes = mesh.getNodes();
  const CRConnectivity& faceNodes = mesh.getAllFaceNodes();
    
  const int count = faces.getCount();
  shared_ptr<VectorT3Array> faPtr(new VectorT3Array(count));
  VectorT3Array& fa = *faPtr;
    
  const T half(0.5);
  const VectorT3Array& nodeCoord =
    dynamic_cast<const VectorT3Array&>(_coordField[nodes]);
    
  for(int f=0; f<count; f++)
  {
      const int numNodes = faceNodes.getCount(f);
        
      if (numNodes == 2)
      {
          const int n0 = faceNodes(f,0);
          const int n1 = faceNodes(f,1);
          VectorT3 dr = nodeCoord[n1]-nodeCoord[n0];
          fa[f][0] = dr[1];
          fa[f][1] = -dr[0];
          fa[f][2] = 0.;
      }
      else if (numNodes == 3)
      {
          const int n0 = faceNodes(f,0);
          const int n1 = faceNodes(f,1);
          const int n2 = faceNodes(f,2);
          VectorT3 dr10 = nodeCoord[n1]-nodeCoord[n0];
          VectorT3 dr20 = nodeCoord[n2]-nodeCoord[n0];
          fa[f] = half*cross(dr10,dr20);
      }
      else if (numNodes == 4)
      {
          const int n0 = faceNodes(f,0);
          const int n1 = faceNodes(f,1);
          const int n2 = faceNodes(f,2);
          const int n3 = faceNodes(f,3);
          VectorT3 dr20 = nodeCoord[n2]-nodeCoord[n0];
          VectorT3 dr31 = nodeCoord[n3]-nodeCoord[n1];
          fa[f] = half*cross(dr20,dr31);
      }
      else if (numNodes > 0)
      {
          fa[f].zero();
          for (int nn=0; nn<numNodes; nn++)
          {
              const int n0 = faceNodes(f,nn);
              const int n1 = faceNodes(f,(nn+1)%numNodes);
              VectorT3 xm = T(0.5)*(nodeCoord[n1]+nodeCoord[n0]);
              VectorT3 dr = (nodeCoord[n1]-nodeCoord[n0]);

              fa[f][0] += xm[1]*dr[2];
              fa[f][1] += xm[2]*dr[0];
              fa[f][2] += xm[0]*dr[1];
          }
      }
  }
      
  _areaField.addArray(faces,faPtr);
}


  
template<class T>
void
MeshMetricsCalculator<T>::calculateFaceAreaMag(const Mesh& mesh)
{
  const StorageSite& faces = mesh.getFaces();
  const VectorT3Array& faceArea =
    dynamic_cast<const VectorT3Array&>(_areaField[faces]);

  const int count = faces.getCount();
  shared_ptr<TArray> famPtr(new TArray(count));
  TArray& fam = *famPtr;

  for(int f=0; f<count; f++)
  {
      fam[f] = mag(faceArea[f]);
  }

  _areaMagField.addArray(faces,famPtr);
}
      
  
template<class T>
void
MeshMetricsCalculator<T>::calculateCellVolumes(const Mesh& mesh)
{
  const StorageSite& faces = mesh.getFaces();
  const StorageSite& cells = mesh.getCells();

  const int cellCount = cells.getCount();
  shared_ptr<TArray> vPtr(new TArray(cellCount));
  TArray& cellVolume = *vPtr;
      
  const VectorT3Array& faceCentroid =
    dynamic_cast<const VectorT3Array&>(_coordField[faces]);
  const VectorT3Array& cellCentroid =
    dynamic_cast<const VectorT3Array&>(_coordField[cells]);
  const VectorT3Array& faceArea =
    dynamic_cast<const VectorT3Array&>(_areaField[faces]);
  const CRConnectivity& faceCells = mesh.getAllFaceCells();
      
  const T dim(mesh.getDimension());
  const int faceCount = faces.getCount();

  cellVolume.zero();

   
  for(int f=0; f<faceCount; f++)
  {
      const int c0 = faceCells(f,0);
      cellVolume[c0] += dot(faceCentroid[f]-cellCentroid[c0],faceArea[f])/dim;
      const int c1 = faceCells(f,1);
      cellVolume[c1] -= dot(faceCentroid[f]-cellCentroid[c1],faceArea[f])/dim;
  }
    
  T volumeSum(0);
  for(int c=0; c<cells.getSelfCount(); c++)
    volumeSum += cellVolume[c];

  cout << "volume sum for Mesh " << mesh.getID() << " = " << volumeSum << endl;


  // update boundary cells with adjacent interior cells values
  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
  {
      const FaceGroup& fg = *fgPtr;
      const StorageSite& faces = fg.site;
      const CRConnectivity& faceCells = mesh.getFaceCells(faces);
      const int faceCount = faces.getCount();
          
      for(int f=0; f<faceCount; f++)
      {
          const int c0 = faceCells(f,0);
          const int c1 = faceCells(f,1);
          cellVolume[c1] = cellVolume[c0];
      }
  }
  _volumeField.addArray(cells,vPtr);
}
  
    
template<class T>
void
MeshMetricsCalculator<T>::computeIBInterpolationMatrices
(const Mesh& mesh,
 const StorageSite& mpmParticles)
{
  typedef CRMatrixTranspose<T,T,T> IMatrix;
  
  const StorageSite& ibFaces = mesh.getIBFaces();
  const StorageSite& cells = mesh.getCells();
  const StorageSite& faces = mesh.getFaces();
  const CRConnectivity& ibFaceToCells = mesh.getConnectivity(ibFaces,cells);
  const CRConnectivity& ibFaceToParticles
    = mesh.getConnectivity(ibFaces,mpmParticles);

  const VectorT3Array& xFaces =
    dynamic_cast<const VectorT3Array&>(_coordField[faces]);

  const VectorT3Array& xCells =
    dynamic_cast<const VectorT3Array&>(_coordField[cells]);

  const VectorT3Array& xParticles =
    dynamic_cast<const VectorT3Array&>(_coordField[mpmParticles]);

  const Array<int>& ibFCRow = ibFaceToCells.getRow();
  const Array<int>& ibFCCol = ibFaceToCells.getCol();

  const Array<int>& ibFPRow = ibFaceToParticles.getRow();
  const Array<int>& ibFPCol = ibFaceToParticles.getCol();
  
  const int nIBFaces = ibFaces.getCount();

  const Array<int>& ibFaceIndices = mesh.getIBFaceList();
  
  shared_ptr<IMatrix> cellToIB(new IMatrix(ibFaceToCells));
  shared_ptr<IMatrix> particlesToIB(new IMatrix(ibFaceToParticles));

  Array<T>& cellToIBCoeff = cellToIB->getCoeff();
  Array<T>& particlesToIBCoeff = particlesToIB->getCoeff();

#if 0
  /******distance weighted interpolation******/

  for(int n=0; n<nIBFaces; n++)
  {
      const int f = ibFaceIndices[n];
      T wtSum(0);
      int nnb(0);
      
      for(int nc=ibFCRow[n]; nc<ibFCRow[n+1]; nc++)
      {
          const int c = ibFCCol[nc];
          VectorT3 dr(xCells[c]-xFaces[f]);
          T wt = 1.0/dot(dr,dr);
          cellToIBCoeff[nc] = wt;
          wtSum += wt;
          nnb++;
      }
      for(int np=ibFPRow[n]; np<ibFPRow[n+1]; np++)
      {
          const int p = ibFPCol[np];
          VectorT3 dr(xParticles[p]-xFaces[f]);
          T wt = 1.0/dot(dr,dr);
          particlesToIBCoeff[np] = wt;
          wtSum += wt;
          nnb++;
      }

      if (nnb == 0)
	throw CException("no cell or particle neighbors for ib face");
      
      for(int nc=ibFCRow[n]; nc<ibFCRow[n+1]; nc++)
      {
	cellToIBCoeff[nc] /= wtSum;
      }

      for(int np=ibFPRow[n]; np<ibFPRow[n+1]; np++)
      {
	particlesToIBCoeff[np] /= wtSum;
      }
  }
#endif

#if 1

  /**********linear least square interpolation*********/
  // X=x-xf  Y=y-yf  Z=Z-zf
  //matrix M=[1 X1 Y1 Z1]
  //         [1 X2 Y2 Z2]
  //         ...........
  //         [1 Xn Yn Zn]
  //coefficient matrix A=[a, b, c, d]T
  //velocity element v = M * A
  //linear relation v = a + b*X + c*Y + d*Z
  //to make least square
  //matrix A = (M(T)*M)^(-1)*M(T)*v
  //so, velocity at face is vface = a + b*Xf + c*Yf + d*Zf = a
  //which is the first row of matrix A
  //so, the weight function should be the first row of  (M(T)*M)^(-1)*M(T)
  //note Q = M(T)*M  and Qinv = Q^(-1)
  //the following code is to calculate it
  //insteading of doing full matrix operation, only nessesary operation on entries are performed

 

  for(int n=0; n<nIBFaces; n++)
  {
      const int f = ibFaceIndices[n];
      T wt(0);
       int nnb(0);
       //calculate Q (4x4)
      T Q[4][4];
      T Qinv[4][4];

      for(int i=0; i<4; i++){
	for(int j=0; j<4; j++){
	  Q[i][j]=Qinv[i][j]=0;
	}
      }
      

       for(int nc=ibFCRow[n]; nc<ibFCRow[n+1]; nc++)
      {
          const int c = ibFCCol[nc];
          VectorT3 dr(xCells[c]-xFaces[f]);
	  //cout<<n<<" cells "<<c<<" "<<xCells[c]<<endl;
	  //cout<<n<<" faces "<<f<<" "<<xFaces[f]<<endl;
	  Q[0][0] += 1.0;
	  Q[0][1] += dr[0];
	  Q[0][2] += dr[1];
	  Q[0][3] += dr[2];
	  Q[1][1] += dr[0]*dr[0];
	  Q[1][2] += dr[0]*dr[1];
	  Q[1][3] += dr[0]*dr[2];
	  Q[2][2] += dr[1]*dr[1];
	  Q[2][3] += dr[1]*dr[2];
	  Q[3][3] += dr[2]*dr[2];          
          nnb++;
      }
      for(int np=ibFPRow[n]; np<ibFPRow[n+1]; np++)
      {
          const int p = ibFPCol[np];
          VectorT3 dr(xParticles[p]-xFaces[f]);
	  // cout<<n<<" particles "<<p<<" "<<xParticles[p]<<endl;
	  Q[0][0] += 1.0;
	  Q[0][1] += dr[0];
	  Q[0][2] += dr[1];
	  Q[0][3] += dr[2];
	  Q[1][1] += dr[0]*dr[0];
	  Q[1][2] += dr[0]*dr[1];
	  Q[1][3] += dr[0]*dr[2];
	  Q[2][2] += dr[1]*dr[1];
	  Q[2][3] += dr[1]*dr[2];
	  Q[3][3] += dr[2]*dr[2];      
	  nnb++;
      }
      
      if (nnb < 4)
	throw CException("not enough cell or particle neighbors for ib face to interpolate!");

      //symetric matrix
      for(int i=0; i<4; i++){
	for(int j=0; j<i; j++){
	  Q[i][j]=Q[j][i];
	}
      }
     
      //calculate the inverse of Q(4x4)
     
      matrix<T> matrix;
      matrix.Invert4x4(Q, Qinv);

     


      //calculate Qinv*M(T) get the first row element, put in coeffMatrix
       for(int nc=ibFCRow[n]; nc<ibFCRow[n+1]; nc++)
      {
          const int c = ibFCCol[nc];
          VectorT3 dr(xCells[c]-xFaces[f]);
	  wt = Qinv[0][0];
	  for (int i=1; i<=3; i++){
	    wt += Qinv[0][i]*dr[i-1];
	  }
	  cellToIBCoeff[nc] = wt;
	  //cout<<n<<" cells "<<nc<<" "<<cellToIBCoeff[nc]<<endl;
      }
       for(int np=ibFPRow[n]; np<ibFPRow[n+1]; np++)
      {
          const int p = ibFPCol[np];
          VectorT3 dr(xParticles[p]-xFaces[f]);
	  wt = Qinv[0][0];
	  for (int i=1; i<=3; i++){
	    wt += Qinv[0][i]*dr[i-1];
	  }
	  particlesToIBCoeff[np] = wt;
	  // cout<<n<<" particles  "<<np<<" "<<particlesToIBCoeff[np]<<endl;
      }
  }
#endif

  GeomFields::SSPair key1(&ibFaces,&cells);
  this->_geomFields._interpolationMatrices[key1] = cellToIB;

  GeomFields::SSPair key2(&ibFaces,&mpmParticles);
  this->_geomFields._interpolationMatrices[key2] = particlesToIB;
  
}

template<class T>
MeshMetricsCalculator<T>::MeshMetricsCalculator(GeomFields& geomFields,const MeshList& meshes) :
  Model(meshes),
  _geomFields(geomFields),
  _coordField(geomFields.coordinate),
  _areaField(geomFields.area),
  _areaMagField(geomFields.areaMag),
  _volumeField(geomFields.volume)
{
  logCtor();
}
  

template<class T>
void
MeshMetricsCalculator<T>::init()
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
  {
      const Mesh& mesh = *_meshes[n];
      calculateNodeCoordinates(mesh);
      calculateFaceAreas(mesh);
      calculateFaceAreaMag(mesh);
      calculateFaceCentroids(mesh);
      calculateCellCentroids(mesh);
      calculateCellVolumes(mesh);
  }

  _volumeField.syncLocal();
  _coordField.syncLocal();
}
//***********************************************************************//

template<class T>
void
MeshMetricsCalculator<T>::computeGridInterpolationMatrices
(const Mesh& mesh,
 const StorageSite& grids, 
 const StorageSite& faces)
{

  typedef CRMatrixTranspose<T,T,T> IMatrix;
  
  const CRConnectivity& faceToGrids = mesh.getConnectivity(faces, grids );

  const VectorT3Array& xFaces =
    dynamic_cast<const VectorT3Array&>(_coordField[faces]);

  const VectorT3Array& xGrids =
    dynamic_cast<const VectorT3Array&>(_coordField[grids]);

  const Array<int>& FGRow = faceToGrids.getRow();
  const Array<int>& FGCol = faceToGrids.getCol();
  
  //  const int nGrids = grids.getCount();

  const int nFaces = faces.getCount();
  
  shared_ptr<IMatrix> gridToFaces (new IMatrix(faceToGrids));
 
  Array<T>& gridToFaceCoeff = gridToFaces->getCoeff(); 

#if 0
  /* distance weighted */
  for(int n=0; n<nFaces; n++)
  {
     
      T wtSum(0);
      int nnb(0);
      
      for(int nc = FGRow[n]; nc < FGRow[n+1]; nc ++)
      {
	  
          const int c = FGCol[nc ];	 
          VectorT3 dr(xGrids[c]-xFaces[n]);	 
          T wt = 1.0/dot(dr,dr);
          gridToFaceCoeff[nc ] = wt;
          wtSum += wt;
          nnb++;
      }
     
     
      
      for(int nc=FGRow[n]; nc<FGRow[n+1]; nc++)
      {
          gridToFaceCoeff[nc] /= wtSum;	
      }
     
  }
#endif 


#if 1
   /* linear interpolation */
   // v = a + b*x +c*y
   //solve a, b, c 
   T Q[3][3];
   T Qinv[3][3];

 
   for(int n=0; n<nFaces; n++)
     {
      T wt(0);
      int nnb(0);
         
      for(int i=0; i<3; i++){
	for(int j=0; j<3; j++){
	  Q[i][j]=Qinv[i][j]=0;
	}
      }
      //      int size = FGRow[n+1]-FGRow[n];
     
      // if(size == 3){
      for(int nc=FGRow[n]; nc<FGRow[n+1]; nc++)
	{
          const int c = FGCol[nc];
	  VectorT3 dr(xGrids[c]-xFaces[n]);
	  //VectorT3 dr(xGrids[c]);
	  Q[nnb][0]=1.0;
	  Q[nnb][1]=dr[0];
	  Q[nnb][2]=dr[1];
          nnb++;
      }
       
      matrix<T> matrix;
      matrix.Invert3x3(Q, Qinv);
      nnb = 0;
      for(int nc=FGRow[n]; nc<FGRow[n+1]; nc++)
	{	  
	  //wt = Qinv[0][nnb]+xFaces[n][0]*Qinv[1][nnb]+xFaces[n][1]*Qinv[2][nnb];
	  wt = Qinv[0][nnb];
	  gridToFaceCoeff[nc] = wt;
	  nnb++;

	}
      //}
     }
 
 
 
#endif



#if 0
   //bi-linear average 
   for(int n=0; n<nFaces; n++)
     {
       const int ncSize = FGRow[n+1] - FGRow[n];
       T dX[ncSize];
       T dY[ncSize];
       T dZ[ncSize];
       T wt[ncSize];
       int count = 0;
       for(int nc=FGRow[n]; nc<FGRow[n+1]; nc++)
	{
          const int c = FGCol[nc];
          VectorT3 dr(xGrids[c]-xFaces[n]);
	  dX[count] = fabs(dr[0]);
	  dY[count] = fabs(dr[1]);
	  dZ[count] = fabs(dr[2]);
	  count++;	  
	}
       const T u = dX[0]/(dX[0]+dX[2]);
       const T v = dY[0]/(dY[0]+dY[1]);
       wt[3] = u*v;
       wt[2] = u*(1.-v);
       wt[1] = (1.-u)*v;
       wt[0] = (1.-u)*(1.-v);
       count = 0;
       for(int nc=FGRow[n]; nc<FGRow[n+1]; nc++)
       {
	 gridToFaceCoeff[nc] = wt[count];
	 count++;
       }
       
     }


#endif
 
  GeomFields::SSPair key1(&faces,&grids);
  this->_geomFields._interpolationMatrices[key1] = gridToFaces;
  
}
//**************************************************************************//

template<class T>
void
MeshMetricsCalculator<T>::computeIBInterpolationMatrices(const StorageSite& p)
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
  {
      const Mesh& mesh = *_meshes[n];
      computeIBInterpolationMatrices(mesh,p);
  }
}
template<class T>
void
MeshMetricsCalculator<T>::computeGridInterpolationMatrices(const StorageSite& g,const StorageSite& f )
{
  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
  {
      const Mesh& mesh = *_meshes[n];
      computeGridInterpolationMatrices(mesh,g, f);
  }
}


template<class T>
MeshMetricsCalculator<T>::~MeshMetricsCalculator()
{
  logDtor();
}

#ifdef USING_ATYPE_TANGENT
template<class T>
void
MeshMetricsCalculator<T>::setTangentCoords(int meshID, int faceGroupID,
                                           int dim)
{
  const Mesh& mesh = *_meshes[meshID];

  VectorT3Array& nodeCoords =
    dynamic_cast<VectorT3Array&>(_coordField[mesh.getNodes()]);
  
  bool found = false;
  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
  {
      const FaceGroup& fg = *fgPtr;
      if (fg.id == faceGroupID)
      {
          const StorageSite& faces = fg.site;
          const CRConnectivity& faceNodes = mesh.getFaceNodes(faces);
          const int nFaces = faces.getCount();
          for(int f=0; f<nFaces; f++)
          {
              const int nFaceNodes = faceNodes.getCount(f);
              for(int nf=0; nf<nFaceNodes; nf++)
              {
                  VectorT3& xn = nodeCoords[faceNodes(f,nf)];
                  xn[dim]._dv=1.0;
                  //xn[dim]._v=1e-7;
              }
          }
          found = true;
      }
  }

  if (!found)
    throw CException("setTangentCoords: invalid faceGroupID");

  const int numMeshes = _meshes.size();
  for (int n=0; n<numMeshes; n++)
  {
      const Mesh& mesh = *_meshes[n];
      calculateFaceAreas(mesh);
      calculateFaceAreaMag(mesh);
      calculateFaceCentroids(mesh);
      calculateCellCentroids(mesh);
      calculateCellVolumes(mesh);
  }

}

#endif
