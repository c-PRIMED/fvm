// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

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
        VectorT3Array& dirichletNodeDisplacement =
          dynamic_cast<VectorT3Array&> (_geomFields.dirichletNodeDisplacement[nodes]);
	VectorT3Array& nodeNormal = 
          dynamic_cast<VectorT3Array&> (_geomFields.boundaryNodeNormal[boundaryNodes]);          
	const Array<int>& displacementOptions =
	  dynamic_cast<Array<int>& > (_geomFields.displacementOptions[nodes]);
	
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
		if(displacementOptions[j] == 0)
		{
		    nodeDisplacement[j] = T(0.0);
		    nodeCoordinate[j] += nodeDisplacement[j] - (*(_previousNodeDisplacementPtr))[j];
		}
		else if(displacementOptions[j] == 1)
		{
		    averageDirichletDisplacement += mag((dirichletNodeDisplacement)[j]);
                    nDirichlet ++;
                    
		    nodeDisplacement[j] = (dirichletNodeDisplacement)[j];
		    nodeCoordinate[j] += nodeDisplacement[j] - (*(_previousNodeDisplacementPtr))[j];
		}
		else if(displacementOptions[j] == 2)
		{
		    T temp = dot(dr,nodeNormal[GlobalToLocal[j]]);
		    nodeDisplacement[j] = dr - temp*nodeNormal[GlobalToLocal[j]];
		    //if(dot(nodeDisplacement[j],nodeNormal[GlobalToLocal[j]])!=zero)
		    //  nodeDisplacement[j] = dr + temp*nodeNormal[GlobalToLocal[j]];
		    nodeDisplacement[j] = (*(_previousNodeDisplacementPtr))[j] +
                      underrelaxation*(nodeDisplacement[j]-(*(_previousNodeDisplacementPtr))[j]);
		    nodeCoordinate[j] += nodeDisplacement[j] - (*(_previousNodeDisplacementPtr))[j];
		}
		else if(displacementOptions[j] == 3)
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
            //cout<<"\nsweep  "<<i<<" max change is "<< maxChangeInDisplacement
            //    <<" and ratio is "<< maxChangeRelative<<"\n";
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
    _meshes(meshes)
  {
    logCtor();
  }
                        
  virtual ~MovingMeshModel() {}
  
  MovingMeshModelOptions<T>& getOptions() {return _options;}

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

	  const T deltaT = _options["timeStep"];

	  if (mesh.getDimension() == 2)
	  {

	      faceVel.zero();
	      gridVel.zero();
	      volChangeDot.zero();
	      sweptVolDot.zero();
	      gridFlux.zero();
	      //	      const T deltaT = _options["timeStep"];
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
	  }
	  else if (mesh.getDimension() == 3)
	  {

	      faceVel.zero();
	      gridVel.zero();
	      volChangeDot.zero();
	      sweptVolDot.zero();
	      gridFlux.zero();
	      //	      const T deltaT = _options["timeStep"];
	      const T half(0.5);
	      const T onepointfive(1.5);
	      const T onesixth(1./6.);
	      const T onethird(1./3.);
	      T temp(0.);
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
		  if (numNodes == 3)
		  {
		      const int n0 = faceNodes(f,0);
		      const int n1 = faceNodes(f,1);
		      const int n2 = faceNodes(f,2);
		      VectorT3 dr10 = nodeCoord[n1]-nodeCoord[n0];
		      VectorT3 dr20 = nodeCoord[n2]-nodeCoord[n0];
                      VectorT3 dr10N1 = nodeCoordN1[n1]-nodeCoordN1[n0];
                      VectorT3 dr20N1 = nodeCoordN1[n2]-nodeCoordN1[n0];
		      VectorT3 eta = onesixth*(cross(dr10,dr20)+cross(dr10N1,dr20N1)
					+half*(cross(dr10,dr20N1)+cross(dr10N1,dr20)));
                      VectorT3 del0 = nodeCoord[n0]-nodeCoordN1[n0];
                      VectorT3 del1 = nodeCoord[n1]-nodeCoordN1[n1];
                      VectorT3 del2 = nodeCoord[n2]-nodeCoordN1[n2];
		      VectorT3 avg = onethird*(del0 + del1 + del2);
		      temp = dot(avg,eta)/deltaT;
		  }
                  else if (numNodes == 4)
		  {
                      const int n0 = faceNodes(f,0);
                      const int n1 = faceNodes(f,1);
                      const int n2 = faceNodes(f,2);
                      const int n3 = faceNodes(f,3);

                      VectorT3 del0 = nodeCoord[n0]-nodeCoordN1[n0];
                      VectorT3 del1 = nodeCoord[n1]-nodeCoordN1[n1];
                      VectorT3 del2 = nodeCoord[n2]-nodeCoordN1[n2];
                      VectorT3 del3 = nodeCoord[n3]-nodeCoordN1[n3];

                      VectorT3 dr10 = nodeCoord[n1]-nodeCoord[n0];
                      VectorT3 dr20 = nodeCoord[n2]-nodeCoord[n0];
                      VectorT3 dr10N1 = nodeCoordN1[n1]-nodeCoordN1[n0];
                      VectorT3 dr20N1 = nodeCoordN1[n2]-nodeCoordN1[n0];
                      VectorT3 eta1 = onesixth*(cross(dr10,dr20)+cross(dr10N1,dr20N1)
					       +half*(cross(dr10,dr20N1)+cross(dr10N1,dr20)));
                      VectorT3 avg1 = onethird*(del0+del1+del2);
                      T temp1 = dot(avg1,eta1)/deltaT;

		      //                      VectorT3 dr20 = nodeCoord[n2]-nodeCoord[n0];
                      VectorT3 dr30 = nodeCoord[n3]-nodeCoord[n0];
		      //                      VectorT3 dr20N1 = nodeCoordN1[n2]-nodeCoordN1[n0];
                      VectorT3 dr30N1 = nodeCoordN1[n3]-nodeCoordN1[n0];
                      VectorT3 eta2 = onesixth*(cross(dr20,dr30)+cross(dr20N1,dr30N1)
						+half*(cross(dr20,dr30N1)+cross(dr20N1,dr30)));
                      VectorT3 avg2 = onethird*(del0+del2+del3);
                      T temp2 = dot(avg2,eta2)/deltaT;

		      temp = temp1+temp2;
		  }
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

	  }

          // update boundary cells with adjacent interior cells values
          foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
	  {
	      const FaceGroup& fg = *fgPtr;
	      if(fg.groupType != "interior")
	      {
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
          shared_ptr<Array<int> > displacementOptions(new Array<int>(nodes.getCount()));
          *displacementOptions = 3;
          _geomFields.displacementOptions.addArray(nodes,displacementOptions);
          shared_ptr<VectorT3Array> dirichletNodeDisplacement
	    (new VectorT3Array(nodes.getCount()));
          dirichletNodeDisplacement->zero();
          _geomFields.dirichletNodeDisplacement.addArray(nodes,dirichletNodeDisplacement);
          shared_ptr<VectorT3Array> nodeDisplacement
            (new VectorT3Array(nodes.getCount()));
          nodeDisplacement->zero();
          _geomFields.nodeDisplacement.addArray(nodes,nodeDisplacement);
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
};


#endif

