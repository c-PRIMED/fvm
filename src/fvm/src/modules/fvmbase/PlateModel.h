// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _PLATEMODEL_H_
#define _PLATEMODEL_H_

#include "Model.h"

#include "GeomFields.h"
#include "PlateFields.h"

#include "Mesh.h"
#include "LinearSolver.h"

#include "PlateBC.h"

template<class T>
class PlateModel : public Model
{
public:

  typedef std::map<int,PlateBC<T>*> PlateBCMap;
  typedef std::map<int,PlateVC<T>*> PlateVCMap;
  
  class Impl;
  
  
  PlateModel(const GeomFields& geomFields,
               PlateFields& plateFields, const MeshList& meshes);
  
  virtual ~PlateModel();

  virtual void init();

  virtual map<string,shared_ptr<ArrayBase> >&
  getPersistenceData();

  virtual void restart(); 

  PlateBCMap& getBCMap();
  PlateVCMap& getVCMap();

  PlateModelOptions<T>& getOptions();

  void printBCs();

  // do the specified number of iterations, return true if converged 
  bool advance(const int niter);

  void updateTime();
  void getMoment(const Mesh& mesh);

  void dumpMatrix(const string fileBase);
private:
  shared_ptr<Impl> _impl;
};

#endif

