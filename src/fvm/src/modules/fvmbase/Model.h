#ifndef _MODEL_H_
#define _MODEL_H_


#include "StorageSite.h"
#include "Model.h"

class Model
{
public:
  Model(const MeshList& meshes);
  virtual ~Model();
    
  DEFINE_TYPENAME("Model");

protected:
  const MeshList& _meshes;
  StorageSiteList _varSites;
  StorageSiteList _fluxSites;  
};

#endif
