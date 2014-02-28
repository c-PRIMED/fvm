// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "FloatVarDict.h"
#include "AMG.h"

template<class T>
struct FractureBC : public FloatVarDict<T>
{
  FractureBC()
  {
      this->defineVar("specifiedPhaseFieldValue",T(1.0));
      this->defineVar("specifiedPhaseFieldFlux",T(0.0));
      //this->defineVar("convectiveCoefficient", T(0.0));
      //this->defineVar("farFieldTemperature", T(300.0));
  }
  string bcType;
};

template<class T>
struct FractureVC : public FloatVarDict<T>
{
  FractureVC()
  {
      this->defineVar("fractureConductivity",T(1.0));
      this->defineVar("fractureSource", T(0.0));
      this->defineVar("fractureSourceCoef", T(0.0));
      	
  }
  string vcType;
};

template<class T>
struct FractureModelOptions : public FloatVarDict<T>
{
  FractureModelOptions()
  {
    this->defineVar("initialPhaseFieldValue",T(1.0));
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

