// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "FloatVarDict.h"
#include "AMG.h"

template<class T>
struct BatterySpeciesBC : public FloatVarDict<T>
{
  BatterySpeciesBC()
  {
      this->defineVar("specifiedConcentration",T(0.0));
      this->defineVar("specifiedMassFlux",T(1.0));
      
  }
  string bcType;
};

template<class T>
struct BatterySpeciesVC : public FloatVarDict<T>
{
  BatterySpeciesVC()
  {
      this->defineVar("massDiffusivity",T(1.e-9));
      this->defineVar("initialConcentration",T(1.0));
  }
  string vcType;
};

template<class T>
struct BatteryPotentialBC : public FloatVarDict<T>
{
  BatteryPotentialBC()
  {
      this->defineVar("specifiedPotential",T(300.0));
      this->defineVar("specifiedPotentialFlux",T(0.0));
  }
  string bcType;
};

template<class T>
struct BatteryPotentialVC : public FloatVarDict<T>
{
  BatteryPotentialVC()
  {
      this->defineVar("conductivity",T(1.0));
      this->defineVar("initialPotential",T(0.0));
  }
  string vcType;
};

template<class T>
struct BatteryThermalBC : public FloatVarDict<T>
{
  BatteryThermalBC()
  {
      this->defineVar("specifiedTemperature",T(300.0));
      this->defineVar("specifiedHeatFlux",T(0.0));
  }
  string bcType;
};

template<class T>
struct  BatteryThermalVC : public FloatVarDict<T>
{
   BatteryThermalVC()
  {
      this->defineVar("thermalConductivity",T(1.0));
      this->defineVar("density", T(1.0));
      this->defineVar("specificHeat", T(1.0));
      this->defineVar("initialTemperature",T(300.0));   	
  }
  string vcType;
};

template<class T>
struct BatteryModelOptions : public FloatVarDict<T>
{
  BatteryModelOptions()
  {
    this->defineVar("ButlerVolmerRRConstant",T(5.0e-7));
    this->defineVar("ButlerVolmerAnodeShellMeshID", int(-1));
    this->defineVar("ButlerVolmerCathodeShellMeshID", int(-1));
    this->defineVar("BatteryElectrolyteMeshID", int(-1));
    this->defineVar("timeStep",T(0.1));
    this->defineVar("interfaceSpeciesUnderRelax",T(1.0));
    this->relativeTolerance=1e-8;
    this->absoluteTolerance=1e-16;
    this->relativeSpeciesTolerance=1e-8;
    this->absoluteSpeciesTolerance=1e-16;
    this->relativePotentialTolerance=1e-8;
    this->absolutePotentialTolerance=1e-16;
    this->relativeThermalTolerance=1e-8;
    this->absoluteThermalTolerance=1e-16;
    this->relativePCTolerance=1e-8;
    this->absolutePCTolerance=1e-16;
    this->linearSolver = 0;
    this->linearSolverSpecies = 0;
    this->linearSolverPotential = 0;
    this->linearSolverThermal = 0;
    this->linearSolverPC = 0;
    this->useCentralDifference=false;
    this->transient = true;
    this->ButlerVolmer = false;
    this->thermalModelPC = true;
    this->timeDiscretizationOrder=2;
    this->advanceVerbosity=1;
  }
  double relativeTolerance;
  double absoluteTolerance;
  double relativeSpeciesTolerance;
  double absoluteSpeciesTolerance;
  double relativePotentialTolerance;
  double absolutePotentialTolerance;
  double relativeThermalTolerance;
  double absoluteThermalTolerance;
  double relativePCTolerance;
  double absolutePCTolerance;
  bool useCentralDifference;
  LinearSolver *linearSolver;
  LinearSolver *linearSolverSpecies;
  LinearSolver *linearSolverPotential;
  LinearSolver *linearSolverThermal;
  LinearSolver *linearSolverPC;
  bool transient;
  bool ButlerVolmer;
  bool thermalModelPC;
  int timeDiscretizationOrder;
  int advanceVerbosity;

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

  LinearSolver& getLinearSolverSpecies()
  {
    if (this->linearSolverSpecies == 0)
    {
        LinearSolver* ls(new  AMG());
        ls->relativeTolerance = 1e-1;
        ls->nMaxIterations = 20;
        ls->verbosity=0;
        this->linearSolverSpecies = ls;
    }
    return *this->linearSolverSpecies ;
  }
LinearSolver& getLinearSolverPotential()
  {
    if (this->linearSolverPotential == 0)
    {
        LinearSolver* ls(new  AMG());
        ls->relativeTolerance = 1e-1;
        ls->nMaxIterations = 20;
        ls->verbosity=0;
        this->linearSolverPotential = ls;
    }
    return *this->linearSolverPotential ;
  }

LinearSolver& getLinearSolverThermal()
  {
    if (this->linearSolverThermal == 0)
    {
        LinearSolver* ls(new  AMG());
        ls->relativeTolerance = 1e-1;
        ls->nMaxIterations = 20;
        ls->verbosity=0;
        this->linearSolverThermal = ls;
    }
    return *this->linearSolverThermal ;
  }

LinearSolver& getLinearSolverPC()
  {
    if (this->linearSolverPC == 0)
    {
        LinearSolver* ls(new  AMG());
        ls->relativeTolerance = 1e-1;
        ls->nMaxIterations = 20;
        ls->verbosity=0;
        this->linearSolverPC = ls;
    }
    return *this->linearSolverPC ;
  }
#endif
};

