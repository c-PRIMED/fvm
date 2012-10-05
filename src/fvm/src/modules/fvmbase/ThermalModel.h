// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _THERMALMODEL_H_
#define _THERMALMODEL_H_

#include "Model.h"

#include "GeomFields.h"
#include "ThermalFields.h"

#include "Mesh.h"


#include "ThermalBC.h"

template<class T>
class ThermalModel : public Model
{
public:

  typedef std::map<int,ThermalBC<T>*> ThermalBCMap;
  typedef std::map<int,ThermalVC<T>*> ThermalVCMap;
  class Impl;
  
  
  ThermalModel(const GeomFields& geomFields,
               ThermalFields& thermalFields, const MeshList& meshes);
  
  virtual ~ThermalModel();

  virtual void init();
  
  ThermalBCMap& getBCMap();
  ThermalVCMap& getVCMap();
  
  ThermalBC<T>& getBC(const int id);

  ThermalModelOptions<T>& getOptions();
  
  void computeIBFaceTemperature(const StorageSite& particles);

  T getHeatFluxIntegral(const Mesh& mesh, const int faceGroupId);
  
  void printBCs();
#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC))
  void dumpMatrix(const string fileBase);
#endif
  void advance(const int niter);

  void updateTime();
private:
  shared_ptr<Impl> _impl;
};

#endif

