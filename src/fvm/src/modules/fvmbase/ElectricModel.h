#ifndef _ELECTRICMODEL_H_
#define _ELECTRICMODEL_H_

#include "Model.h"
#include "GeomFields.h"
#include "ElectricFields.h"
#include "Mesh.h"
#include "ElectricBC.h"
#include "Octree.h"
#include "Array.h"
#include "Vector.h"
#include "PhysicsConstant.h"

template<class T>
class ElectricModel : public Model
{
public:

  typedef std::map<int,ElectricBC<T>*> ElectricBCMap;
  typedef std::map<int,ElectricVC<T>*> ElectricVCMap;

  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Vector<double, 3> VectorD3;
  typedef Array<VectorT3> VectorT3Array;

  class Impl;
  
  
  ElectricModel(const GeomFields& geomFields,
               ElectricFields& electricFields, const MeshList& meshes);
  
  virtual ~ElectricModel();

  virtual void init();
  
  ElectricBCMap& getBCMap();
  ElectricBC<T>& getBC(const int id);

  ElectricVCMap& getVCMap();
  ElectricVC<T>& getVC(const int id);

  ElectricModelOptions<T>& getOptions();

  ElectricModelConstants<T>& getConstants();
  
  //computeIBFacePotential(const StorageSite& particles);
    
  void printBCs();

  bool advance(const int niter);

  void updateTime();

  void calculateEquilibriumParameters();

  void  dielectricOneDimModelPrep(const int nXCol, const int nYCol, const int nGrid,
				  const VectorD3 corner1_1, 
				  const VectorD3 corner1_2, 
				  const VectorD3 corner1_3,
				  const VectorD3 corner1_4,
				  const VectorD3 corner2_1, 
				  const VectorD3 corner2_2, 
				  const VectorD3 corner2_3,
				  const VectorD3 corner2_4);
  

 

/*
  const int findClosestPoint(const VectorT3 point, Octree& O);

  boost::shared_ptr<ArrayBase> createPathAndDiscretize (
      const VectorT3 startPoint, const VectorT3 endPoint, const int N);

  boost::shared_ptr<ArrayBase > pathPointsInterpolation
    (boost::shared_ptr<ArrayBase> pathPoints, 
     Octree& O,
     const double radius);
*/


private:
  shared_ptr<Impl> _impl;
};

#endif

