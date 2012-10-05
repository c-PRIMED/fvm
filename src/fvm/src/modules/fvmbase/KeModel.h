// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _KEMODEL_H_
#define _KEMODEL_H_

#include "Model.h"

#include "GeomFields.h"
#include "KeFields.h"
#include "FlowFields.h"
#include "Mesh.h"


#include "KeBC.h"

template<class T>
class KeModel : public Model
{
public:

  typedef std::map<int,KeBC<T>*> KeBCMap;
  typedef std::map<int,KeVC<T>*> KeVCMap;
  class Impl;


  KeModel(const GeomFields& geomFields,
          KeFields& keFields,
          FlowFields& flowFields,
          const MeshList& meshes);

  virtual ~KeModel();

  virtual void init();

  KeBCMap& getBCMap();
  KeVCMap& getVCMap();

  KeBC<T>& getBC(const int id);

  KeModelOptions<T>& getOptions();

  void getViscosity(const Mesh& mesh);
  void updateTimek();
  void updateTimee();
  void printBCs();
/*
#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC))
  void dumpMatrix(const string fileBase);
#endif
*/
  void advance(const int niter);

private:
  shared_ptr<Impl> _impl;
};

#endif

