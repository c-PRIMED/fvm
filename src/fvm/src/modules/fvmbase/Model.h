#ifndef _MODEL_H_
#define _MODEL_H_


#include "StorageSite.h"
#include "Mesh.h"

class Model
{
public:
  Model(const MeshList& meshes);
  virtual ~Model();
    
  DEFINE_TYPENAME("Model");

  virtual void init() = 0;
  
protected:
  const MeshList _meshes;
  StorageSiteList _varSites;
  StorageSiteList _fluxSites;  
};

#endif
