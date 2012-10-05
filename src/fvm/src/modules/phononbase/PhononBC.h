// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _PHONONBC_H_
#define _PHONONBC_H_

#include "misc.h"
#include "FloatVarDict.h"
#include "AMG.h"

template<class T>
struct PhononBC : public FloatVarDict<T>
{
  PhononBC()
  {
      this->defineVar("specifiedTemperature",T(300.0));
      this->defineVar("specifiedReflection",T(0.0));
      this->defineVar("specifiedHeatflux",T(0.0));
  }
  string bcType;
};

template<class T>
struct PhononModelOptions : public FloatVarDict<T>
{
  PhononModelOptions()
  {
    this->defineVar("timeStep",T(0.1));
    this->defineVar("initialTemperature",T(310.0));
    this->defineVar("Tref",T(299.0));
    this->printNormalizedResiduals = true;
    this->transient = false;
    this->timeDiscretizationOrder=1;
    this->PhononLinearSolver=0;
    this->absTolerance=1e-8;
    this->relTolerance=1-4;
    this->showResidual=5;
    this->defineVar("transmissivity0to1", T(0.0));
    this->defineVar("transmissivity1to0", T(0.0));	   
  }
  
  bool printNormalizedResiduals;
  bool transient;
  int timeDiscretizationOrder;
  LinearSolver* PhononLinearSolver;
  double absTolerance;
  double relTolerance;
  int showResidual;

#ifndef SWIG
  LinearSolver& getPhononLinearSolver()
  {
    if (this->PhononLinearSolver == 0)
    {
        LinearSolver* ls(new AMG());
        ls->relativeTolerance = 1e-5;
        ls->nMaxIterations = 100;
        ls->verbosity=0;
        this->PhononLinearSolver = ls;
    }
    return *this->PhononLinearSolver ;
  }
#endif

};


#endif
