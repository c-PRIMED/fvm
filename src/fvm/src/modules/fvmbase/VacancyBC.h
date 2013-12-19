// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "FloatVarDict.h"
#include "AMG.h"

template<class T>
struct VacancyBC : public FloatVarDict<T>
{
  VacancyBC()
  {
      this->defineVar("specifiedConcentration",T(300.0));
      this->defineVar("specifiedVacaFlux",T(0.0));
      this->defineVar("convectiveCoefficient", T(0.0));
      this->defineVar("farFieldConcentration", T(300.0));
  }
  string bcType;
};

template<class T>
struct VacancyVC : public FloatVarDict<T>
{
  VacancyVC()
  {
      this->defineVar("vacancyDiffusioncoefficient",T(1.0));
      this->defineVar("density", T(1.0));
      this->defineVar("specificVaca", T(1.0));
      	
  }
  string vcType;
};

template<class T>
struct VacancyModelOptions : public FloatVarDict<T>
{
  VacancyModelOptions()
  {
    this->defineVar("initialConcentration",T(300.0));
    this->defineVar("timeStep", T(1e-7));
    this->relativeTolerance=1e-8;
    this->absoluteTolerance=1e-16;
    this->linearSolver = 0;
    this->useCentralDifference=false;
    this->transient=false;
    this->timeDiscretizationOrder = 1;
  }
  double relativeTolerance;
  double absoluteTolerance;
  bool useCentralDifference;
  LinearSolver *linearSolver;
  bool transient;
  int timeDiscretizationOrder;
#ifndef SWIG
  LinearSolver& getLinearSolver()
  {
    if (this->linearSolver == 0)
    {
        LinearSolver* ls(new  AMG());
        ls->relativeTolerance = 1e-1;
        ls->nMaxIterations = 20;
        ls->verbosity=0;
        this->linearSolver = ls;
    }
    return *this->linearSolver ;
  }
#endif
};

