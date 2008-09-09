
#include "Model.h"

Model::Model(const MeshList& meshes) :
  _meshes(meshes),
  _varSites(),
  _fluxSites()
{}

Model::~Model()
{}

