// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _PLATEDEFORMATIONMODEL_H_
#define _PLATEDEFORMATIONMODEL_H_

#include "Model.h"
#include "GeomFields.h"
#include "PlateFields.h"
#include "NumType.h"
#include "Mesh.h"
#include "Array.h"
#include "Field.h"
#include "Vector.h"
#include "CRConnectivity.h"
#include "StorageSite.h"


template<class T>
class PlateDeformationModel : public Model
{
public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Vector<T,3> VecD3;
  typedef Array<VectorT3> VectorT3Array;  
 
  PlateDeformationModel(GeomFields& geomFields,
			PlateFields& plateFields, 
			const MeshList& meshes) :
    Model(meshes),
    _geomFields(geomFields),
    _plateFields(plateFields),
    _meshes(meshes)
  {
    logCtor();
  }
                        
  virtual ~PlateDeformationModel() {}
  
  void calculateNodeDisplacement()
  {
      const int numMeshes = _meshes.size();
      for (int m=0;m<numMeshes;m++)
      {
          const Mesh& mesh = *_meshes[m];
	  const StorageSite& nodes = mesh.getNodes();
	  const StorageSite& cells = mesh.getCells();
	  const int nNodes = nodes.getCount();
	  const CRConnectivity& cellNodes = mesh.getCellNodes();
	  shared_ptr<CRConnectivity> nodeCellsPtr = cellNodes.getTranspose();
	  CRConnectivity& nodeCells = *nodeCellsPtr;
          const VectorT3Array& nodeCoordinate0 =
            dynamic_cast<const VectorT3Array&>(_geomFields.coordinate0[nodes]);
	  const VectorT3Array& nodeCoordinate =
	    dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[nodes]);
          const VectorT3Array& cellCentroid =
            dynamic_cast<const VectorT3Array&> (_geomFields.coordinate[cells]);
          const VectorT3Array& cellDisplacement =
            dynamic_cast<const VectorT3Array&>(_plateFields.deformation[cells]);
          VectorT3Array& nodeDisplacement =
            dynamic_cast<VectorT3Array&>(_geomFields.nodeDisplacement[nodes]);
	  const T one(1.0);
          const T zero(0.0);

	  for(int j=0;j<nNodes;j++)
          {
	      VectorT3 dr(NumTypeTraits<VectorT3>::getZero());
	      T weight(0.0);
	      for(int k=0;k<nodeCells.getCount(j);k++)
              {
		  const int num = nodeCells(j,k);
		  const VectorT3 ds = cellCentroid[num]-nodeCoordinate[j];
		  dr += cellDisplacement[num]/mag(ds);
		  weight += one/mag(ds);
	      }
	      dr = dr/weight;
	      //nodeDisplacement[j][0] = nodeCoordinate0[j][2]*nodeDisplacement[j][0];
	      //nodeDisplacement[j][1] = nodeCoordinate0[j][2]*nodeDisplacement[j][1];
	      nodeDisplacement[j][0] = zero;
              nodeDisplacement[j][1] = zero;
	      nodeDisplacement[j][2] = dr[2];
	  }
      }
  }

  void deformPlate()
  {
      const int numMeshes = _meshes.size();
      for (int n=0;n<numMeshes;n++)
      {
	  const Mesh& mesh = *_meshes[n];
	  const StorageSite& nodes = mesh.getNodes();
          const int nNodes = nodes.getCount();
       	  VectorT3Array& nodeCoord =
            dynamic_cast<VectorT3Array&>(_geomFields.coordinate[nodes]);
          //Array<VecD3>& nodeCoordMesh = mesh.getNodeCoordinates();
	  const VectorT3Array& nodeCoord0 =
            dynamic_cast<const VectorT3Array&>(_geomFields.coordinate0[nodes]);
          VectorT3Array& nodeDisplacement =
            dynamic_cast<VectorT3Array&>(_geomFields.nodeDisplacement[nodes]);
          VectorT3Array& nodeCoordK1 =
            dynamic_cast<VectorT3Array&>(_geomFields.coordinateK1[nodes]);
	  nodeCoordK1 = nodeCoord;
	  for (int i=0;i<nNodes;i++)
	  {
	      nodeCoord[i] = nodeCoord0[i] + nodeDisplacement[i];
	      //nodeCoordMesh[i] = nodeCoord[i]; 
	  }
	  
      }
  }

  void deformMeshStructure()
  {
      const int numMeshes = _meshes.size();
      for (int n=0;n<numMeshes;n++)
      {
	  const Mesh& mesh = *_meshes[n];
	  const StorageSite& nodes = mesh.getNodes();
	  const int nNodes = nodes.getCount();
	  const Array<int>& displacementOptions =
	    dynamic_cast<Array<int>& > (_geomFields.displacementOptions[nodes]);
	  VectorT3Array& nodeCoord =
	    dynamic_cast<VectorT3Array&>(_geomFields.coordinate[nodes]);
	  const VectorT3Array& nodeCoord0 =
	    dynamic_cast<const VectorT3Array&>(_geomFields.coordinate0[nodes]);
	  VectorT3Array& nodeDisplacement =
	    dynamic_cast<VectorT3Array&>(_geomFields.nodeDisplacement[nodes]);
          VectorT3Array& nodeCoordK1 =
            dynamic_cast<VectorT3Array&>(_geomFields.coordinateK1[nodes]);
          nodeCoordK1 = nodeCoord;	  
	  for (int i=0;i<nNodes;i++)
	  {
	      if(displacementOptions[i] == 3)
	      {
		  nodeCoord[i] = nodeCoord0[i] + nodeDisplacement[i];
	      }
	  }
      }
  }


  void updateBoundaryMesh(const Mesh& mesh, Mesh& bMesh, const double thickness)
  {
    // update node coords
    const StorageSite& bMeshNodes = bMesh.getNodes();
    const int nBNodes = bMeshNodes.getCount();
    
    const StorageSite& nodes = mesh.getNodes();

    const VectorT3Array& nodeCoord =
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[nodes]);

    VectorT3Array& bNodeCoord =
      dynamic_cast<VectorT3Array&>(_geomFields.coordinate[bMeshNodes]);

    VectorT3Array& bMeshCoord = bMesh.getNodeCoordinates();
    // the node x,y coordinates keep the same; 
    // z coordinate gets updated by beam node coodinate;
    for(int n=0; n<nBNodes/2; n++)
    {
        bNodeCoord[n][2] = -thickness/2 + nodeCoord[n][2];
        bMeshCoord[n][2] = -thickness/2 + nodeCoord[n][2];
    }
   
    for(int n=nBNodes/2; n<nBNodes; n++)
    {
        bNodeCoord[n][2] = thickness/2 + nodeCoord[n-nBNodes/2][2];
        bMeshCoord[n][2] = thickness/2 + nodeCoord[n-nBNodes/2][2];
    }
  }

  
  void updateBoundaryMesh(const Mesh& mesh, Mesh& bMesh,
			  const double thickness, 
			  const double timeStep, 
			  Field& velocityField)
  {
    updateBoundaryMesh(mesh, bMesh, thickness);

    // update face velocity
    const StorageSite& nodes = mesh.getNodes();
    const StorageSite& bMeshFaces = bMesh.getFaces();

    shared_ptr<VectorT3Array>
      bVelocity(new VectorT3Array(bMeshFaces.getCount()));

    velocityField.addArray(bMeshFaces,bVelocity);

    VectorT3Array& w =
      dynamic_cast<VectorT3Array&>(_geomFields.nodeDisplacement[nodes]);
    VectorT3Array& wN1 =
      dynamic_cast<VectorT3Array&>(_geomFields.nodeDisplacementN1[nodes]);
    
    const StorageSite& cells = mesh.getCells();
    const int nLocalCells = cells.getSelfCount();
    const int nCells = cells.getCount();
    
    for (int c=0; c<nLocalCells; c++){
        const CRConnectivity& cellNodes = mesh.getCellNodes();
	const int nCellNodes = cellNodes.getCount(c);
	VectorT3 vf(VectorT3::getZero());
	for (int n=0; n<nCellNodes; n++){
	  const int node = cellNodes(c, n);
	  VectorT3 vn = (w[node]-wN1[node])/timeStep;
	  vf += vn;
	}
	vf /= nCellNodes;
	
	(*bVelocity)[c] = vf;
	(*bVelocity)[c+nLocalCells] = vf;
    }

    for (int c=nLocalCells; c<nCells; c++){
        const CRConnectivity& cellNodes = mesh.getCellNodes();
	const int nCellNodes = cellNodes.getCount(c);
	VectorT3 vf(VectorT3::getZero());
	for (int n=0; n<nCellNodes; n++){
	  const int node = cellNodes(c, n);
	  VectorT3 vn = (w[node]-wN1[node])/timeStep;
	  vf += vn;
	}
	vf /= nCellNodes;
	
	(*bVelocity)[c+nLocalCells] = vf;
    }  

   
  }
  
  
  const ArrayBase& getCommon(const StorageSite& site, const StorageSite& osite)
  {
    
    const StorageSite::CommonMap& commonMap = site.getCommonMap();    
    cout<<"\n The size of commonMap is "<<commonMap.size()<<"\n";

    const Array<int>& myNodeIndices =
      *(site.getCommonMap().find(&osite)->second);
    return myNodeIndices;
  }
    
  void init() 
  {
      const int numMeshes = _meshes.size();
      for (int n=0;n<numMeshes;n++)
      {
	  const Mesh& mesh = *_meshes[n];
	  const StorageSite& nodes = mesh.getNodes();
	  VectorT3Array& nodeCoord =
            dynamic_cast<VectorT3Array&>(_geomFields.coordinate[nodes]);
          _geomFields.coordinateK1.addArray(nodes,
	    dynamic_pointer_cast<ArrayBase>(nodeCoord.newCopy()));

          _geomFields.coordinate0.addArray(nodes,
	    dynamic_pointer_cast<ArrayBase>(nodeCoord.newCopy()));

          shared_ptr<VectorT3Array>
            nodeDisplacement(new VectorT3Array(nodes.getCount()));
          nodeDisplacement->zero();
          _geomFields.nodeDisplacement.addArray(nodes,nodeDisplacement);

          shared_ptr<VectorT3Array>
            nodeDisplacementN1(new VectorT3Array(nodes.getCount()));
          nodeDisplacementN1->zero();

          _geomFields.nodeDisplacementN1.addArray(nodes,
                                                  nodeDisplacementN1);
      }
  }

  void updateTime()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& nodes = mesh.getNodes();
        VectorT3Array& w =
          dynamic_cast<VectorT3Array&>(_geomFields.nodeDisplacement[nodes]);
        VectorT3Array& wN1 =
          dynamic_cast<VectorT3Array&>(_geomFields.nodeDisplacementN1[nodes]);
	wN1 = w;
    }
  }

private:
  GeomFields& _geomFields;
  PlateFields& _plateFields;
  const MeshList _meshes;
};


#endif

