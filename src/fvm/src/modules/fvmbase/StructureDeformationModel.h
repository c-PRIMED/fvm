// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _STRUCTUREDEFORMATIONMODEL_H_
#define _STRUCTUREDEFORMATIONMODEL_H_

#include "Model.h"
#include "GeomFields.h"
#include "StructureFields.h"
#include "NumType.h"
#include "Mesh.h"
#include "Array.h"
#include "Field.h"
#include "Vector.h"
#include "CRConnectivity.h"
#include "StorageSite.h"


template<class T>
class StructureDeformationModel : public Model
{
public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;  
 
  StructureDeformationModel(GeomFields& geomFields,
                            StructureFields& structureFields, 
			    const MeshList& meshes) :
    Model(meshes),
    _geomFields(geomFields),
    _structureFields(structureFields),
    _meshes(meshes)
  {
    logCtor();
  }
                        
  virtual ~StructureDeformationModel() {}
  
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
	  const VectorT3Array& nodeCoordinate =
	    dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[nodes]);
          const VectorT3Array& cellCentroid =
            dynamic_cast<const VectorT3Array&> (_geomFields.coordinate[cells]);
          const VectorT3Array& cellDisplacement =
            dynamic_cast<const VectorT3Array&>(_structureFields.deformation[cells]);
          VectorT3Array& nodeDisplacement =
            dynamic_cast<VectorT3Array&>(_geomFields.nodeDisplacement[nodes]);
	  const T one(1.0);

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
	      nodeDisplacement[j] = dr;
	  }
      }
  }

  void deformStructure()
  {
      const int numMeshes = _meshes.size();
      for (int n=0;n<numMeshes;n++)
      {
	  const Mesh& mesh = *_meshes[n];
	  const StorageSite& nodes = mesh.getNodes();
          const int nNodes = nodes.getCount();
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
	      nodeCoord[i] = nodeCoord0[i] + nodeDisplacement[i];
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


  // update node coordinates and face velocity for the boundary mesh
  // corresponding to the given mesh
  void updateBoundaryMesh(const Mesh& mesh, Mesh& bMesh,
                          Field& velocityField,
                          const double timeStep)
  {
    // update node coords
    const StorageSite& bMeshNodes = bMesh.getNodes();
    
    
    const StorageSite& nodes = mesh.getNodes();

    const VectorT3Array& nodeCoord =
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[nodes]);

    VectorT3Array& bNodeCoord =
      dynamic_cast<VectorT3Array&>(_geomFields.coordinate[bMeshNodes]);

    VectorT3Array& bMeshCoord = bMesh.getNodeCoordinates();

    const Array<int>& myNodeIndices =
      *(bMeshNodes.getCommonMap().find(&nodes)->second);
    const Array<int>& bNodeIndices =
      *(nodes.getCommonMap().find(&bMeshNodes)->second);

    const int nCommon = myNodeIndices.getLength();
    for(int n=0; n<nCommon; n++)
    {
        bNodeCoord[bNodeIndices[n]] = nodeCoord[myNodeIndices[n]];
        bMeshCoord[bNodeIndices[n]] = nodeCoord[myNodeIndices[n]];
    }
   


    // update face velocity
    
    const StorageSite& bMeshFaces = bMesh.getFaces();

    shared_ptr<VectorT3Array>
      bVelocity(new VectorT3Array(bMeshFaces.getCount()));

    velocityField.addArray(bMeshFaces,bVelocity);

    VectorT3Array& w =
      dynamic_cast<VectorT3Array&>(_geomFields.nodeDisplacement[nodes]);
    VectorT3Array& wN1 =
      dynamic_cast<VectorT3Array&>(_geomFields.nodeDisplacementN1[nodes]);

    int bMeshFaceCount=0;

    foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
    {
        const FaceGroup& fg = *fgPtr;
        const StorageSite& faces = fg.site;
        if (fg.groupType!="interior")
        {
            const int nFaces = faces.getCount();
            const CRConnectivity& faceNodes = mesh.getFaceNodes(faces);
            for(int f=0; f<nFaces; f++)
            {
                const int nFaceNodes = faceNodes.getCount(f);
                VectorT3 vf(VectorT3::getZero());
                
                for(int nn=0; nn<nFaceNodes; nn++)
                {
                    const int node = faceNodes(f,nn);
                    VectorT3 vn = (w[node] - wN1[node])/timeStep;
                    vf += vn;
                }
                vf /= nFaceNodes;
                (*bVelocity)[bMeshFaceCount] = vf;
                 bMeshFaceCount++;
            }
        }
    }
  }
  

  // update node coordinates and face velocity for the boundary mesh
  // corresponding to the given mesh
  void updateBoundaryMesh(const Mesh& mesh, Mesh& bMesh,
                          Field& velocityField, const map<int,int>& commonFacesMap,
                          const double timeStep)
  {
    // update face velocity
    const StorageSite& bMeshFaces = bMesh.getFaces();
    shared_ptr<VectorT3Array>
      bVelocity(new VectorT3Array(bMeshFaces.getCount()));
      bVelocity->zero();

    velocityField.addArray(bMeshFaces,bVelocity);
  
    const StorageSite& nodes = mesh.getNodes();
    VectorT3Array& w =
      dynamic_cast<VectorT3Array&>(_geomFields.nodeDisplacement[nodes]);
    VectorT3Array& wN1 =
      dynamic_cast<VectorT3Array&>(_geomFields.nodeDisplacementN1[nodes]);

    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups()){
       const FaceGroup& fg = *fgPtr;
       const StorageSite& faces = fg.site;
       const int nFaces = faces.getCount();
       const CRConnectivity& faceNodes = mesh.getFaceNodes(faces);
       for(int f=0; f<nFaces; f++){
          const int faceID = f + faces.getOffset();
          const int nFaceNodes = faceNodes.getCount(f);
          VectorT3 vf(VectorT3::getZero());
              
          for(int nn=0; nn<nFaceNodes; nn++){
             const int node = faceNodes(f,nn);
             VectorT3 vn = (w[node] - wN1[node])/timeStep;
             vf += vn;
          }
          vf /= nFaceNodes;
          (*bVelocity)[ commonFacesMap.find(faceID)->second ] = vf;
       }
    }
  } 

  
  const ArrayBase& getCommon(const StorageSite& site, const StorageSite& osite)
  {
    
    const StorageSite::CommonMap& commonMap = site.getCommonMap();    
    //cout<<"\n The size of commonMap is "<<commonMap.size()<<"\n";

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
  StructureFields& _structureFields;
  const MeshList _meshes;
};


#endif

