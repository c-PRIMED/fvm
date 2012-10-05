// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.


#include "GradientModel.h"

map<const Mesh*, shared_ptr<GradientMatrixBase> >
GradientModelBase::_gradientMatricesMap;

void
GradientModelBase::clearGradientMatrix(const Mesh& mesh)
{
  if (_gradientMatricesMap.find(&mesh) != _gradientMatricesMap.end())
    _gradientMatricesMap.erase(&mesh);
}
