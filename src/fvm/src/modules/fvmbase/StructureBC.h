// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "misc.h"
#include "FloatVarDict.h"
#include "TractionVal.h"
#include "AMG.h"

template<class T>
struct StructureBC : public FloatVarDict<T>
{
  StructureBC()
  {
      this->defineVar("specifiedXDeformation",T(0.0));
      this->defineVar("specifiedYDeformation",T(0.0));
      this->defineVar("specifiedZDeformation",T(0.0));
      this->defineVar("specifiedXXTraction",T(0.0));
      this->defineVar("specifiedXYTraction",T(0.0));
      this->defineVar("specifiedXZTraction",T(0.0));
      this->defineVar("specifiedYXTraction",T(0.0));
      this->defineVar("specifiedYYTraction",T(0.0));
      this->defineVar("specifiedYZTraction",T(0.0));
      this->defineVar("specifiedZXTraction",T(0.0));
      this->defineVar("specifiedZYTraction",T(0.0));
      this->defineVar("specifiedZZTraction",T(0.0));
      this->defineVar("specifiedXForce",T(0.0));
      this->defineVar("specifiedYForce",T(0.0));
      this->defineVar("specifiedZForce",T(0.0));
      this->defineVar("specifiedXDistForce",T(0.0));
      this->defineVar("specifiedYDistForce",T(0.0));
      this->defineVar("specifiedZDistForce",T(0.0));
  }
  string bcType;
};

template<class T>
struct StructureVC : public FloatVarDict<T>
{
  StructureVC()
  {
      this->defineVar("eta",T(1.0));
      this->defineVar("eta1",T(1.0));
      this->defineVar("density",T(1.0));
      this->defineVar("alpha",T(1.0));
  }
  string vcType;
};


template<class T>
struct StructureModelOptions : public FloatVarDict<T>
{
  StructureModelOptions()
  {
    this->defineVar("initialXDeformation",T(0.0));
    this->defineVar("initialYDeformation",T(0.0));
    this->defineVar("initialZDeformation",T(0.0));
    this->defineVar("deformationURF",T(0.7));
    this->defineVar("timeStep",T(0.1));
    this->defineVar("operatingTemperature",T(300.0));
    this->defineVar("residualXXStress",T(0.));
    this->defineVar("residualYYStress",T(0.));
    this->defineVar("residualZZStress",T(0.));

    this->deformationTolerance=1e-4;
    this->printNormalizedResiduals = true;
    this->transient = false;
    this->timeDiscretizationOrder=1;
    this->variableTimeStep = false;
    this->timeStepN1=0.1;
    this->timeStepN2=0.1;
    this->thermo = false;
    this->residualStress = false;
    this->deformationLinearSolver = 0;
    this->coupledLinearSolver = 0;
    this->creep = false;
    this->creepModel=1;
    this->A = 0.03/3600.;
    this->B = 1.8e8;
    this->m = 2.0;
    this->n = 2.0;
    this->Sy0 = 1.0e9;

    this->incompressible = true;
  }
  
  bool printNormalizedResiduals;
  double deformationTolerance;
  bool transient;
  int timeDiscretizationOrder;
  bool variableTimeStep;
  double timeStepN1;
  double timeStepN2;
  bool thermo;
  bool residualStress;
  LinearSolver *deformationLinearSolver;
  LinearSolver *coupledLinearSolver;
  bool creep;
  int creepModel;
  double A;
  double B;
  double m;
  double n;
  double Sy0;

  bool incompressible;
#ifndef SWIG
  LinearSolver& getDeformationLinearSolver()
  {
    if (this->deformationLinearSolver == 0)
    {
        LinearSolver* ls(new AMG());
        ls->relativeTolerance = 1e-1;
        ls->nMaxIterations = 20;
        ls->verbosity=0;
        this->deformationLinearSolver = ls;
    }
    return *this->deformationLinearSolver ;
  }
#endif
};

