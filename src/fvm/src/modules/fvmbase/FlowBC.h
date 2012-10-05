// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "misc.h"
#include "FloatVarDict.h"
#include "AMG.h"

template<class T>
struct FlowBC : public FloatVarDict<T>
{
  FlowBC()
  {
      this->defineVar("specifiedXVelocity",T(0.0));
      this->defineVar("specifiedYVelocity",T(0.0));
      this->defineVar("specifiedZVelocity",T(0.0));
      this->defineVar("specifiedPressure",T(0.0));
      this->defineVar("accomodationCoefficient",T(1.0));
  }
  string bcType;
};

template<class T>
struct FlowVC : public FloatVarDict<T>
{
  FlowVC()
  {
      this->defineVar("viscosity",T(1e-3));
      this->defineVar("density",T(1.0));
      this->defineVar("eddyviscosity",T(1e-5));
      this->defineVar("totalviscosity",T(2e-3));
  }
  string vcType;
};


template<class T>
struct FlowModelOptions : public FloatVarDict<T>
{
  FlowModelOptions()
  {
    this->defineVar("initialXVelocity",T(0.0));
    this->defineVar("initialYVelocity",T(0.0));
    this->defineVar("initialZVelocity",T(0.0));
    this->defineVar("initialPressure",T(0.0));
    this->defineVar("momentumURF",T(0.7));
    this->defineVar("velocityURF",T(1.0));
    this->defineVar("pressureURF",T(0.3));
    this->defineVar("timeStep",T(0.1));
    this->defineVar("operatingPressure",T(101325.0));
    this->defineVar("operatingTemperature",T(300.0));
    this->defineVar("molecularWeight",T(28.966));

    this->momentumTolerance=1e-3;
    this->continuityTolerance=1e-3;
    this->printNormalizedResiduals = true;
    this->transient = false;
    this->correctVelocity = true;
    this->timeDiscretizationOrder=1;
    this->momentumLinearSolver = 0;
    this->pressureLinearSolver = 0;
    this->coupledLinearSolver = 0;

    this->incompressible = true;
    this->turbulent = false;
    this->vk = 0.4187;
    this->emp = 9.793;
    this-> cmu = 0.09;
   

  }
  
  bool printNormalizedResiduals;
  double momentumTolerance;
  double continuityTolerance;
  bool transient;
  bool turbulent;
  bool correctVelocity;
  int timeDiscretizationOrder;
  LinearSolver *momentumLinearSolver;
  LinearSolver *pressureLinearSolver;
  LinearSolver *coupledLinearSolver;
  double cmu;
  double vk;
  double emp;
  bool incompressible;
#ifndef SWIG
  LinearSolver& getMomentumLinearSolver()
  {
    if (this->momentumLinearSolver == 0)
    {
        LinearSolver* ls(new AMG());
        ls->relativeTolerance = 1e-1;
        ls->nMaxIterations = 20;
        ls->verbosity=0;
        this->momentumLinearSolver = ls;
    }
    return *this->momentumLinearSolver ;
  }

  LinearSolver& getPressureLinearSolver()
  {
    if (this->pressureLinearSolver == 0)
    {
        LinearSolver* ls(new AMG());
        ls->relativeTolerance = 1e-1;
        ls->nMaxIterations = 20;
        ls->verbosity=0;
        this->pressureLinearSolver = ls;
    }
    return *this->pressureLinearSolver ;
  }
#endif
};

