#ifndef _ELECTRONICMODEL_H_
#define _ELECTRONICMODEL_H_

#include "Model.h"
#include "GeomFields.h"
#include "ElectronicFields.h"
#include "Mesh.h"
#include "ElectronicBC.h"
#include "Octree.h"
#include "Array.h"
#include "Vector.h"

template<class T>
class ElectronicModel : public Model
{
public:

  typedef std::map<int,ElectronicBC<T>*> ElectronicBCMap;
  
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  class Impl;
  
  
  ElectronicModel(const GeomFields& geomFields,
               ElectronicFields& electronicFields, const MeshList& meshes);
  
  virtual ~ElectronicModel();

  virtual void init();
  
  ElectronicBCMap& getBCMap();
  ElectronicBC<T>& getBC(const int id);

  ElectronicModelOptions<T>& getOptions();
  
  //computeIBFacePotential(const StorageSite& particles);
    
  void printBCs();

  bool advance(const int niter);

  void updateTime();

  const int findClosestPoint(const VectorT3 point, Octree& O);

  boost::shared_ptr<ArrayBase> createPathAndDiscretize (
      const VectorT3 startPoint, const VectorT3 endPoint, const int N);

  boost::shared_ptr<ArrayBase > pathPointsInterpolation
    (boost::shared_ptr<ArrayBase> pathPoints, 
     Octree& O,
     const double radius);



private:
  shared_ptr<Impl> _impl;
};

#endif

