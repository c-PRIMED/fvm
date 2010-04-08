
#include "Model.h"

Model::Model(const MeshList& meshes) :
  _meshes(meshes),
  _varSites(),
  _fluxSites()
{}

Model::~Model()
{}

map<string,shared_ptr<ArrayBase> >&
Model::getPersistenceData()
{
  return _persistenceData;
}

void
Model::restart()
{}
 
