#ifndef _DISCRETIZATION_H_
#define _DISCRETIZATION_H_

#include "misc.h"
#include "GlobalFields.h"
#include "Mesh.h"

class MultiFieldMatrix;
class MultiField;
class Model;

class Discretization
{
public:

  Discretization(const Model& model, const MeshList& meshes);

  virtual ~Discretization();

  virtual void discretize(const Mesh& mesh, MultiFieldMatrix& matrix,
                          MultiField& x, MultiField& r) = 0;

  DEFINE_TYPENAME("Discretization");
protected:
  const Model& _model;
  const MeshList& _meshes;
  const Field& _coordField;
  const Field& _areaField;
  const Field& _areaMagField;
  const Field& _volumeField;
};
#endif
