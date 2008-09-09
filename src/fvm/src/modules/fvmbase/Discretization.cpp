
#include "Discretization.h"

Discretization::Discretization(const Model& model, const MeshList& meshes):
  _model(model),
  _meshes(meshes),
  _coordField(getGlobalField(Mesh::coordinate)),
  _areaField(getGlobalField(Mesh::area)),
  _areaMagField(getGlobalField(Mesh::areaMag)),
  _volumeField(getGlobalField(Mesh::volume))
{}


Discretization::~Discretization()
{}
