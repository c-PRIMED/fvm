#ifndef _MOVINGMESHMODEL_H_
#define _MOVINGMESHMODEL_H_

#include "Model.h"
#include "GeomFields.h"
#include "FlowFields.h"
#include "NumType.h"
#include "Mesh.h"
#include "Array.h"
#include "Field.h"
#include "Vector.h"
#include "CRConnectivity.h"
#include "StorageSite.h"
#include "MovingMeshBC.h"


template<class T>
class MovingMeshModel : public Model
{
public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  void advance() 
  {
    const int numMeshes = _meshes.size();
    for(int n=0;n<numMeshes;n++)
    {
	const Mesh& mesh = *_meshes[n];
        const StorageSite& nodes = mesh.getNodes();
	const int nNodes = nodes.getCount();
	shared_ptr<Array<int> > GlobalToLocalPtr = mesh.createAndGetBNglobalToLocal();
        Array<int>& GlobalToLocal = *GlobalToLocalPtr;
	const StorageSite& boundaryNodes = mesh.getBoundaryNodes();
	VectorT3Array& nodeCoordinate = 
          dynamic_cast<VectorT3Array&> (_geomFields.coordinate[nodes]);
	VectorT3Array& nodeDisplacement = 
          dynamic_cast<VectorT3Array&> (_geomFields.nodeDisplacement[nodes]);
	VectorT3Array& nodeNormal = 
          dynamic_cast<VectorT3Array&> (_geomFields.boundaryNodeNormal[boundaryNodes]);          
	
        const CRConnectivity& cellNodes = mesh.getCellNodes();
	shared_ptr<CRConnectivity> nodeCellsPtr = cellNodes.getTranspose();
	shared_ptr<CRConnectivity> nodeNodesPtr = nodeCellsPtr->multiply(cellNodes,false);
	CRConnectivity& nodeNodes = *nodeNodesPtr;
	nodeDisplacement.zero();
        const T one(1.0);
	const T underrelaxation = _options["underrelaxation"];
	const T small(1e-10);

        
        for(int i=0;i<_options.nNodeDisplacementSweeps;i++)
	{

            int nDirichlet =0;
            T averageDirichletDisplacement(0.);

	    shared_ptr<VectorT3Array> _previousNodeDisplacementPtr =
              dynamic_pointer_cast<VectorT3Array>(nodeDisplacement.newCopy());

            for(int j=0;j<nNodes;j++)
	    {
	        VectorT3 dr(NumTypeTraits<VectorT3>::getZero());
	        T weight(0.0);
	        for(int k=0;k<nodeNodes.getCount(j);k++)
		{
		    const int num = nodeNodes(j,k);
		    if(num!=j)
		    {
		        const VectorT3 ds = nodeCoordinate[num]-nodeCoordinate[j];
			if(mag(ds)!=T(0.0))
			{
			    dr += nodeDisplacement[num]/mag(ds);
                            weight += one/mag(ds);
			}
			else
			{
			    dr += nodeDisplacement[num]/small;
                            weight += one/small;
			}
		    }
		}
		dr = dr/weight;
		if((*_displacementOptionsPtr)[j] == 0)
		{
		    nodeDisplacement[j] = T(0.0);
		    nodeCoordinate[j] += nodeDisplacement[j] - (*(_previousNodeDisplacementPtr))[j];
		}
		else if((*_displacementOptionsPtr)[j] == 1)
		{
                    averageDirichletDisplacement += mag((*_dirichletNodeDisplacementPtr)[j]);
                    nDirichlet ++;
                    
		    nodeDisplacement[j] = (*_dirichletNodeDisplacementPtr)[j];
		    nodeCoordinate[j] += nodeDisplacement[j] - (*(_previousNodeDisplacementPtr))[j];
		}
		else if((*_displacementOptionsPtr)[j] == 2)
		{
		    T temp = dot(dr,nodeNormal[GlobalToLocal[j]]);
		    nodeDisplacement[j] = dr - temp*nodeNormal[GlobalToLocal[j]];
		    //if(dot(nodeDisplacement[j],nodeNormal[GlobalToLocal[j]])!=zero)
		    //  nodeDisplacement[j] = dr + temp*nodeNormal[GlobalToLocal[j]];
		    nodeDisplacement[j] = (*(_previousNodeDisplacementPtr))[j] +
                      underrelaxation*(nodeDisplacement[j]-(*(_previousNodeDisplacementPtr))[j]);
		    nodeCoordinate[j] += nodeDisplacement[j] - (*(_previousNodeDisplacementPtr))[j];
		}
		else if((*_displacementOptionsPtr)[j] == 3)
		{
		    nodeDisplacement[j] = dr;
		    nodeDisplacement[j] = (*(_previousNodeDisplacementPtr))[j] +
         	      underrelaxation*(nodeDisplacement[j]-(*(_previousNodeDisplacementPtr))[j]);
		    nodeCoordinate[j] += nodeDisplacement[j] - (*(_previousNodeDisplacementPtr))[j];
		}
	    }

            if (nDirichlet > 0)
              averageDirichletDisplacement /= nDirichlet;
            else
              averageDirichletDisplacement = T(1.);
            
            T maxChangeInDisplacement = T(0.0);
	    for(int j=0;j<nNodes;j++)
	    {
	        const T changeInDisplacement = (mag(nodeDisplacement[j]-
                                              (*(_previousNodeDisplacementPtr))[j]));
                if (maxChangeInDisplacement < changeInDisplacement)
                  maxChangeInDisplacement = changeInDisplacement;
	    }

            T maxChangeRelative = maxChangeInDisplacement / averageDirichletDisplacement;
            cout<<"\nsweep  "<<i<<" max change is "<< maxChangeInDisplacement
                <<" and ratio is "<< maxChangeRelative<<"\n";
	    if((maxChangeInDisplacement<=_options.absTolerance)||
               (maxChangeRelative<=_options.relativeTolerance))
	    {
		
		return;
	    }
	}
    }   	
  }
 
  MovingMeshModel(const MeshList& meshes, GeomFields& geomFields,
                  FlowFields& flowFields) :
    Model(meshes),
    _geomFields(geomFields),
    _flowFields(flowFields),
    _meshes(meshes),
    _displacementOptionsPtr(0),
    _dirichletNodeDisplacementPtr(0)
  {
    logCtor();
  }
                        
  virtual ~MovingMeshModel() {}
  
  MovingMeshModelOptions<T>& getOptions() {return _options;}

  ArrayBase& getDisplacementOptions()
  {
      // assumes number of meshes is 1
      const int numMeshes = _meshes.size();
      for(int n=0;n<numMeshes;n++)
      {
          const Mesh& mesh = *_meshes[n];
	  const StorageSite& nodes = mesh.getNodes();
	  const int nNodes = nodes.getCount();
	  if(!_displacementOptionsPtr)
	    _displacementOptionsPtr = new Array<int>(nNodes);
	  *_displacementOptionsPtr = 3;
      }
      return *_displacementOptionsPtr;
  }
  

  ArrayBase& setNodeDisplacement()
  {
      // assumes number of meshes is 1
      const int numMeshes = _meshes.size();
      for(int n=0;n<numMeshes;n++)
      {
	  const Mesh& mesh = *_meshes[n];
	  const StorageSite& nodes = mesh.getNodes();
	  const int nNodes = nodes.getCount();
	  if(!_dirichletNodeDisplacementPtr)
	    _dirichletNodeDisplacementPtr = new VectorT3Array(nNodes);
	  (*_dirichletNodeDisplacementPtr).zero();
      }
      return *_dirichletNodeDisplacementPtr;
  }

  void volChange()
  {
      const int numMeshes = _meshes.size();
      for (int m=0;m<numMeshes;m++)
      {
          const Mesh& mesh = *_meshes[m];
	  const StorageSite& nodes = mesh.getNodes();
	  const StorageSite& faces = mesh.getFaces();
	  const StorageSite& cells = mesh.getCells();
	  const int nNodes = nodes.getCount();
	  const int nFaces = faces.getCount();
	  const int nCells = cells.getCount();
	  const CRConnectivity& faceNodes = mesh.getAllFaceNodes();
	  const CRConnectivity& faceCells = mesh.getAllFaceCells();
          const TArray& rho =
            dynamic_cast<const TArray&>(_flowFields.density[cells]);
          const VectorT3Array& faceArea =
	    dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
	  VectorT3Array& faceAreaN1 =
	    dynamic_cast<VectorT3Array&>(_geomFields.areaN1[faces]);
	  const VectorT3Array& nodeCoord =
	    dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[nodes]);
	  VectorT3Array& nodeCoordN1 =
	    dynamic_cast<VectorT3Array&>(_geomFields.coordinateN1[nodes]);
	  TArray& cellVolume = 
	    dynamic_cast<TArray&>(_geomFields.volume[cells]);
	  TArray& gridFlux = 
            dynamic_cast<TArray&>(_geomFields.gridFlux[faces]);
          VectorT3Array& faceVel =
            dynamic_cast<VectorT3Array&>(_geomFields.faceVel[faces]);
	  TArray& sweptVolDot = 
	    dynamic_cast<TArray&>(_geomFields.sweptVolDot[faces]);
	  shared_ptr<VectorT3Array> gridVelPtr(new VectorT3Array(nNodes));
	  VectorT3Array& gridVel = *gridVelPtr;
	  TArray volChangeDot(cells.getCount());


	  faceVel.zero();
	  gridVel.zero();
	  volChangeDot.zero();
	  sweptVolDot.zero();
	  gridFlux.zero();
	  const T deltaT = _options["timeStep"];
	  const T half(0.5);
	  const T onepointfive(1.5);
	  for(int i=0;i<nNodes;i++)
	  {
	      VectorT3 dr = (nodeCoord[i] - nodeCoordN1[i]);
              gridVel[i] = dr/deltaT;
	  }
	  for (int f=0;f<nFaces;f++)
	  {
	      const int numNodes = faceNodes.getCount(f);
	      for (int n=0;n<numNodes;n++)
	      {
		  const int num = faceNodes(f,n);
		  faceVel[f] += gridVel[num];
	      }
	      faceVel[f] /= T(numNodes);
	      T temp = dot((half*(faceArea[f]+faceAreaN1[f])),faceVel[f]);
	      sweptVolDot[f] = temp;
	      const int c0 = faceCells(f,0);
	      volChangeDot[c0] += temp;
	      const int c1 = faceCells(f,1);
	      volChangeDot[c1] -= temp;
	      if (_geomFields.sweptVolDotN1.hasArray(faces))
	      {
	          TArray& sweptVolDotN1 = 
	            dynamic_cast<TArray&>(_geomFields.sweptVolDotN1[faces]);
	          gridFlux[f] = onepointfive*sweptVolDot[f] - 
                    half*sweptVolDotN1[f];
	      }
	      else
	        gridFlux[f] = sweptVolDot[f];
	      gridFlux[f] *= half*(rho[c0]+rho[c1]);
	  }
	  T volumeChange(0.0);
	  for(int c=0;c<cells.getSelfCount();c++)
	    volumeChange += volChangeDot[c];
	  volumeChange *= deltaT;
	  cout << "volume change for Mesh " << mesh.getID() << " = " << volumeChange << endl;


          // update boundary cells with adjacent interior cells values
          foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	      const FaceGroup& fg = *fgPtr;
              const StorageSite& bfaces = fg.site;
              const CRConnectivity& bfaceCells = mesh.getFaceCells(bfaces);
              const int faceCount = bfaces.getCount();
              for(int f=0; f<faceCount; f++)
	      {
		  const int c0 = bfaceCells(f,0);
                  const int c1 = bfaceCells(f,1);
                  volChangeDot[c1] = volChangeDot[c0];
	      }
	  }
	  for (int c=0;c<nCells;c++)
	    cellVolume[c] += volChangeDot[c]*deltaT;
      }
  }

  void updateTime()
  {
      const int numMeshes = _meshes.size();
      for (int n=0;n<numMeshes;n++)
      {
	  const Mesh& mesh = *_meshes[n];
          const StorageSite& cells = mesh.getCells();
	  const StorageSite& faces = mesh.getFaces();
	  const StorageSite& nodes = mesh.getNodes();
       	  const TArray& cellVol =
            dynamic_cast<TArray&>(_geomFields.volume[cells]);
	  TArray& cellVolN1 =
            dynamic_cast<TArray&>(_geomFields.volumeN1[cells]);
	  VectorT3Array& nodeCoord =
            dynamic_cast<VectorT3Array&>(_geomFields.coordinate[nodes]);
	  VectorT3Array& nodeCoordN1 =
            dynamic_cast<VectorT3Array&>(_geomFields.coordinateN1[nodes]);
	  const VectorT3Array& faceArea =
            dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
	  VectorT3Array& faceAreaN1 =
            dynamic_cast<VectorT3Array&>(_geomFields.areaN1[faces]);
	  if (_options.timeDiscretizationOrder > 1)
	  {
              TArray& cellVolN2 =
                dynamic_cast<TArray&> (_geomFields.volumeN2[cells]);
	      cellVolN2 = cellVolN1;
              TArray& sweptVolDot =
                dynamic_cast<TArray&>(_geomFields.sweptVolDot[faces]);
	      TArray& sweptVolDotN1 = 
	        dynamic_cast<TArray&>(_geomFields.sweptVolDotN1[faces]);
	      sweptVolDotN1 = sweptVolDot;
	  }
	  cellVolN1 = cellVol;
	  faceAreaN1 = faceArea;
	  nodeCoordN1 = nodeCoord;
      }
  }

  void init() 
  {
      const int numMeshes = _meshes.size();
      for (int n=0;n<numMeshes;n++)
      {
	  const Mesh& mesh = *_meshes[n];
	  const StorageSite& faces = mesh.getFaces();
	  const StorageSite& cells = mesh.getCells();
	  const StorageSite& nodes = mesh.getNodes();
          const TArray& cellVolume =
            dynamic_cast<TArray&>(_geomFields.volume[cells]);
	  _geomFields.volumeN1.addArray(cells,
	    dynamic_pointer_cast<ArrayBase>(cellVolume.newCopy()));
	  VectorT3Array& nodeCoord =
            dynamic_cast<VectorT3Array&>(_geomFields.coordinate[nodes]);
          _geomFields.coordinateN1.addArray(nodes,
	    dynamic_pointer_cast<ArrayBase>(nodeCoord.newCopy()));
	  shared_ptr<TArray> sweptVolume(new TArray(faces.getCount()));
          sweptVolume->zero();
          _geomFields.sweptVolDot.addArray(faces,sweptVolume);
	  shared_ptr<TArray> gridFlux(new TArray(faces.getCount()));
          gridFlux->zero();
          _geomFields.gridFlux.addArray(faces,gridFlux);
          shared_ptr<VectorT3Array> faceVel(new VectorT3Array(faces.getCount()));
          faceVel->zero();
          _geomFields.faceVel.addArray(faces,faceVel);
	  const VectorT3Array& faceArea =
            dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
	  _geomFields.areaN1.addArray(faces,
	    dynamic_pointer_cast<ArrayBase>(faceArea.newCopy()));
          if (_options.timeDiscretizationOrder > 1)
	  {
	      _geomFields.volumeN2.addArray(cells,
                dynamic_pointer_cast<ArrayBase>(cellVolume.newCopy()));
              _geomFields.sweptVolDotN1.addArray(faces,
                dynamic_pointer_cast<ArrayBase>(sweptVolume->newCopy()));
	  }
      }
  }

private:
  GeomFields& _geomFields;
  FlowFields& _flowFields;
  const MeshList _meshes;
  MovingMeshModelOptions<T> _options;
  Array<int>* _displacementOptionsPtr;
  VectorT3Array* _dirichletNodeDisplacementPtr;
};


#endif

