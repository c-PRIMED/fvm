// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _COMETBC_H_
#define _COMETBC_H_

#include "misc.h"
#include "FloatVarDict.h"
#include "AMG.h"

template<class T>
struct COMETBC : public FloatVarDict<T>
{
  COMETBC()
  {
      this->defineVar("specifiedXVelocity",T(0.00));
      this->defineVar("specifiedYVelocity",T(0.0));
      this->defineVar("specifiedZVelocity",T(0.0));
      this->defineVar("specifiedPressure",T(1.0));
      this->defineVar("specifiedDensity",T(1.0));
      this->defineVar("specifiedTemperature",T(1.0));
      this->defineVar("accommodationCoefficient",T(1.0));
      this->defineVar("specifiedMassFlowRate",T(1.0));
      this->defineVar("specifiedTauxx",T(1.0));
      this->defineVar("specifiedTauyy",T(1.0));
      this->defineVar("specifiedTauzz",T(1.0));
      this->defineVar("specifiedTauxy",T(0.0));
      this->defineVar("specifiedTauyz",T(0.0));
      this->defineVar("specifiedTauzx",T(0.0));
  }
  string bcType;
};

template<class T>
struct COMETVC : public FloatVarDict<T>
{
  COMETVC()
  {
      this->defineVar("viscosity",T(1e-3));
      this->defineVar("density",T(1.0));
  }
  string vcType;
};


template<class T>
struct COMETModelOptions : public FloatVarDict<T>
{
  COMETModelOptions()
  {
    this->defineVar("initialXVelocity",T(0.0));
    this->defineVar("initialYVelocity",T(0.0));
    this->defineVar("initialZVelocity",T(0.0));
    this->defineVar("initialPressure",T(0.0));
    this->defineVar("momentumURF",T(0.7));
    this->defineVar("velocityURF",T(1.0));
    this->defineVar("pressureURF",T(0.3));
    this->defineVar("timeStep",T(1E-3)); 

    this->defineVar("rho_init",T(9.28E-6));
    this->defineVar("T_init",T(273.15));
    this->defineVar("nonDimLt",T(1.0));
    this->defineVar("nonDimLx",T(1.0)); 
    this->defineVar("nonDimLy",T(1.0));
    this->defineVar("nonDimLz",T(1.0));


    this->defineVar("operatingPressure",T(101325.0));
    this->defineVar("operatingTemperature",T(300.0));

    this->defineVar("molecularWeight",T(40.0));
    this ->defineVar("Tmuref",T(273.15));
    this ->defineVar("muref",T(2.117e-5));  
    this ->defineVar("mu_w",T(0.81));  
    this->Prandtl=2.0/3.0; 
    this->SpHeatRatio=5.0/3.0; 
    // Argon  2.117e-5,0.81
    // Helium 1.865e-5,0.66
    // Air  1.7116e-5, 0.74
    // Nitrogen 1.781e-5,
    
    // Argon, Helium Pr=2/3, gamma=5/3
    // Air, Nitrogen Pr=3/4, gamma=7/5
    
    this->Tolerance=1e-4;
    this->printNormalizedResiduals = true;
    this->transient = false;
   
    this->fgamma=2;
   

    this->timeDiscretizationOrder=1;
    this->conOrder=1;
    this->showResidual=5;
    this->maxLevels=2;
    this->AgglomerationMethod="FaceArea";
    this->preSweeps=0;
    this->postSweeps=2;
    this->preCoarsestSweeps=1;
    this->postCoarsestSweeps=1;
    this->method=1;
    this->relaxDistribution=0;
    this->underRelaxation=1.0;
    this->minCells=1;
    this->CentralDifference=false;
    this->KineticLinearSolver = 0;
   
    this-> printCellNumber=0;
    this->defineVar("printDirectionNumber",545);
    
    this->NewtonsMethod_ktrial=50;
    this->relativeTolerance=1e-14;
    this->absoluteTolerance=1e-22; 

    this->BoltzmannConstant=1.38e-23;
    this->Planck=6.26068E-34;
    this->epsilon_ES=1e-50;
    this->pi=acos(-1.0);//3.14159;
 
    this->Knq_direction=0;
    //used in Newton's Method for Equilibrium distribution function
    this->defineVar("ToleranceX",T(1e-8));
    this->defineVar("ToleranceF",T(1e-16));
 }
  
  bool printNormalizedResiduals;
  double Tolerance;
  bool transient;
 
  double BoltzmannConstant;
  double Planck;

  double relativeTolerance;
  double absoluteTolerance;

  double epsilon_ES;
  double pi;
  
  int NewtonsMethod_ktrial;
  int timeDiscretizationOrder;
  int conOrder;
  double absTolerance;
  double relTolerance;
  int showResidual;
  int maxLevels;
  double underRelaxation;
  string AgglomerationMethod;
  int preSweeps;
  int postSweeps;
  int preCoarsestSweeps;
  int postCoarsestSweeps;
  int method;
  int relaxDistribution;
  int minCells;
  bool CentralDifference;

  LinearSolver *KineticLinearSolver;
  int printCellNumber;
  int fgamma;
  int Knq_direction;
  double Prandtl;
  double SpHeatRatio;

#ifndef SWIG
  LinearSolver& getKineticLinearSolver()
  {
    if (this->KineticLinearSolver == 0)
    {
        LinearSolver* ls(new AMG());
        ls->relativeTolerance = 1e-15;
        ls->nMaxIterations = 100;
        ls->verbosity=0;
        this->KineticLinearSolver = ls;
    }
    return *this->KineticLinearSolver ;
  }
#endif

};

#endif
