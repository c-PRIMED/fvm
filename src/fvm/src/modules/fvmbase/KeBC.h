// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "FloatVarDict.h"
#include "AMG.h"

template<class T>
struct KeBC : public FloatVarDict<T>
{
  KeBC()
  {
      this->defineVar("Specifiedk",T(1.0));
      this->defineVar("specifiede",T(1.0));

  }
  string bcType;
};

template<class T>
struct KeVC : public FloatVarDict<T>
{
  KeVC()
  {
      this->defineVar("InitialEnergy",T(1.0));
      this->defineVar("InitialDissipation",T(1.0));
      this->defineVar("c1",T(1e-5));
      this->defineVar("c2",T(7.69e-06));


  }
  string vcType;
};

template<class T>
struct KeModelOptions : public FloatVarDict<T>
{
  KeModelOptions()
  {
    this->defineVar("timeStep",T(0.1));
    this->relativeTolerance=1e-7;
    this->absoluteTolerance=1e-7;
    this->defineVar("energyURF",T(0.8));
    this->defineVar("dissipationURF",T(0.8));
    this->linearSolver = 0;
    this->useCentralDifference=false;
    this->transient = false;
    this->timeDiscretizationOrder=1;
    this-> cmu=0.09;
    this-> c2mu=1.92;
    this-> c1mu=1.44;
    this-> vk=0.4187;
    this-> emp=9.793;
    this-> sigmak =1.0;
    this-> sigmae =1.3;
    this->printNormalizedResiduals = true;

 }     

  bool printNormalizedResiduals ;
  bool transient;
  int timeDiscretizationOrder;
  double relativeTolerance;
  double absoluteTolerance;
  double cmu;
  double c2mu;
  double c1mu;
  double vk;
  double emp;
  double sigmak;
  double sigmae;
  bool useCentralDifference;
  LinearSolver *linearSolver;

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

