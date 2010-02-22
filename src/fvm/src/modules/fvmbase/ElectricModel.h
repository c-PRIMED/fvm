#ifndef _ElectricMODEL_H_
#define _ElectricMODEL_H_

#include "Model.h"
#include "GeomFields.h"
#include "ElectricFields.h"
#include "Mesh.h"
#include "ElectricBC.h"
#include "Octree.h"
#include "Array.h"
#include "Vector.h"

template<class T>
class ElectricModel : public Model
{
public:

  typedef std::map<int,ElectricBC<T>*> ElectricBCMap;
  typedef std::map<int,ElectricVC<T>*> ElectricVCMap;

  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  class Impl;
  
  
  ElectricModel(const GeomFields& geomFields,
               ElectricFields& ElectricFields, const MeshList& meshes);
  
  virtual ~ElectricModel();

  virtual void init();
  
  ElectricBCMap& getBCMap();
  ElectricBC<T>& getBC(const int id);

  ElectricVCMap& getVCMap();
  ElectricVC<T>& getVC(const int id);

  ElectricModelOptions<T>& getOptions();
  
  //computeIBFacePotential(const StorageSite& particles);
    
  void printBCs();

  bool advance(const int niter);

  void updateTime();

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

