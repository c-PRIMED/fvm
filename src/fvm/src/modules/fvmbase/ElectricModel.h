// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _ELECTRICMODEL_H_
#define _ELECTRICMODEL_H_

#include "Model.h"
#include "GeomFields.h"
#include "ElectricFields.h"
#include "Mesh.h"
#include "ElectricBC.h"

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

  virtual map<string,shared_ptr<ArrayBase> >& getPersistenceData();
  virtual void restart(); 

  ElectricBCMap& getBCMap();
  ElectricBC<T>& getBC(const int id);

  ElectricVCMap& getVCMap();
  ElectricVC<T>& getVC(const int id);

  ElectricModelOptions<T>& getOptions();

  ElectricModelConstants<T>& getConstants();
  
  void computeIBFacePotential(const StorageSite& solid);

  void computeSolidSurfaceForce(const StorageSite& particles);

  void computeSolidSurfaceForcePerUnitArea(const StorageSite& particles);  

  void printBCs();

  bool advance(const int niter);

  void updateTime();

  void calculateEquilibriumParameters();
  
  vector<T> getTunnelCurrent();

  T getPotentialFluxIntegral(const Mesh& mesh, const int faceGroupId);
 
private:
  shared_ptr<Impl> _impl;
};

#endif

