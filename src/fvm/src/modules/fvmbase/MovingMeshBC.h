// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "misc.h"
#include "FloatVarDict.h"


template<class T>
struct MovingMeshModelOptions : public FloatVarDict<T>
{
  MovingMeshModelOptions()
  {
      this->defineVar("timeStep",T(0.1));
      this->defineVar("underrelaxation",T(1.0));

      this->absTolerance = 1e-4;
      this->nNodeDisplacementSweeps = 20;
      this->relativeTolerance = 1e-1;
      this->timeDiscretizationOrder = 1;
  }
  double absTolerance;
  double relativeTolerance;
  int nNodeDisplacementSweeps;
  int timeDiscretizationOrder;

#ifndef SWIG
#endif
};

