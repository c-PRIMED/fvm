// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _DISCRETIZATION_H_
#define _DISCRETIZATION_H_

#include "misc.h"
#include "Mesh.h"

class MultiFieldMatrix;
class MultiField;
class Model;

class Discretization
{
public:

  Discretization(const MeshList& meshes);

  virtual ~Discretization();

  virtual void discretize(const Mesh& mesh, MultiFieldMatrix& matrix,
                          MultiField& x, MultiField& r) = 0;

  DEFINE_TYPENAME("Discretization");
protected:
  const MeshList& _meshes;
};

typedef vector<shared_ptr<Discretization> > DiscrList;
#endif
