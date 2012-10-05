// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _MODEL_H_
#define _MODEL_H_


#include "StorageSite.h"
#include "Mesh.h"
#include "ArrayBase.h"

class Model
{
public:
  Model(const MeshList& meshes);
  virtual ~Model();
    
  DEFINE_TYPENAME("Model");

  virtual void init() = 0;

  virtual map<string,shared_ptr<ArrayBase> >&
  getPersistenceData();

  virtual void restart();
  
protected:
  const MeshList _meshes;
  StorageSiteList _varSites;
  StorageSiteList _fluxSites;
  map<string,shared_ptr<ArrayBase> > _persistenceData;
};

#endif
