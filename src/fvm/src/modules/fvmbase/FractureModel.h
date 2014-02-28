// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _FRACTUREMODEL_H_
#define _FRACTUREMODEL_H_

#include "Model.h"

#include "GeomFields.h"
#include "FractureFields.h"

#include "Mesh.h"


#include "FractureBC.h"

template<class T>
class FractureModel : public Model
{
public:

  typedef std::map<int,FractureBC<T>*> FractureBCMap;
  typedef std::map<int,FractureVC<T>*> FractureVCMap;
  class Impl;
  
  
  FractureModel(const GeomFields& geomFields,
               FractureFields& fractureFields, const MeshList& meshes);
  
  virtual ~FractureModel();

  virtual void init();
  
  FractureBCMap& getBCMap();
  FractureVCMap& getVCMap();
  
  FractureBC<T>& getBC(const int id);

  FractureModelOptions<T>& getOptions();
  
  //void computeIBFaceTemperature(const StorageSite& particles);

  //T getHeatFluxIntegral(const Mesh& mesh, const int faceGroupId);
  
  void printBCs();
//#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC))
//  void dumpMatrix(const string fileBase);
//#endif
  void advance(const int niter);

  void updateTime();
private:
  shared_ptr<Impl> _impl;
};

#endif

