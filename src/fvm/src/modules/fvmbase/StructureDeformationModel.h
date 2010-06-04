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
  typedef map<const StorageSite*, shared_ptr<ArrayBase> > ArrayMap;
  typedef pair<const StorageSite*, const StorageSite*> EntryIndex;
  
 
  StructureDeformationModel(GeomFields& geomFields,StructureFields& structureFields, 
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
            dynamic_cast<VectorT3Array&>(_structureFields.nodeDisplacement[nodes]);
	  const T one(1.0);
	  const T small(1e-10);


	  for(int j=0;j<nNodes;j++)
          {
	      VectorT3 dr(NumTypeTraits<VectorT3>::getZero());
	      T weight(0.0);
	      for(int k=0;k<nodeCells.getCount(j);k++)
              {
		  const int num = nodeCells(j,k);
		  const VectorT3 ds = cellCentroid[num]-nodeCoordinate[j];
		  if(mag(ds)!=T(0.0))
		  {
		      dr += cellDisplacement[num]/mag(ds);
			  weight += one/mag(ds);
		  }
		  else
		  {
		      dr += nodeDisplacement[num]/small;
		      weight += one/small;
		  }
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
	  VectorT3Array& nodeCoordN1 =
            dynamic_cast<VectorT3Array&>(_geomFields.coordinateN1[nodes]);
          VectorT3Array& nodeDisplacement =
            dynamic_cast<VectorT3Array&>(_structureFields.nodeDisplacement[nodes]);
	  nodeCoordN1 = nodeCoord;
	  for (int i=0;i<nNodes;i++)
	      nodeCoord[i] = nodeCoord[i] + nodeDisplacement[i];
      }
  }
  
  const ArrayBase& getCommon(const StorageSite& site)
  {
    
    const StorageSite::CommonMap& commonMap = site.getCommonMap();    
    cout<<"\n The size of commonMap is "<<commonMap.size()<<"\n";

    foreach(const StorageSite::CommonMap::value_type& mpos, commonMap)
    {
	const StorageSite& oSite = *mpos.first;
	const Array<int>& fromIndices = *(mpos.second);
	cout<<"\n Common Nodes length is  "<<fromIndices.getLength()<<"\n";
	return fromIndices;
    }
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
          _geomFields.coordinateN1.addArray(nodes,
	    dynamic_pointer_cast<ArrayBase>(nodeCoord.newCopy()));
          shared_ptr<VectorT3Array> nodeDisplacement(new VectorT3Array(nodes.getCount()));
          nodeDisplacement->zero();
          _structureFields.nodeDisplacement.addArray(nodes,nodeDisplacement);
      }
  }

private:
  GeomFields& _geomFields;
  StructureFields& _structureFields;
  const MeshList _meshes;
};


#endif

