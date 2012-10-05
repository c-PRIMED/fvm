// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _LINEARIZER_H_
#define _LINEARIZER_H_

#include "misc.h"
class MultiFieldMatrix;
class MultiField;


#include "Mesh.h"
#include "Discretization.h"

class Linearizer
{
public:
  Linearizer();

  virtual void linearize(DiscrList& discretizations,
                         const MeshList& meshes, MultiFieldMatrix& matrix,
                         MultiField& x, MultiField& b);
};

#endif
