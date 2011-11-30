#ifndef _COMETBC_H_
#define _COMETBC_H_

#include "misc.h"
#include "FloatVarDict.h"

template<class T>
struct COMETBC : public FloatVarDict<T>
{
  COMETBC()
  {
      this->defineVar("specifiedTemperature",T(300.0));
      this->defineVar("specifiedReflection",T(0.0));
  }
  string bcType;
};

template<class T>
struct COMETModelOptions : public FloatVarDict<T>
{
  COMETModelOptions()
  {
    this->defineVar("timeStep",T(0.1));
    this->defineVar("initialTemperature",T(310.0));
    this->defineVar("Tref",T(299.0));
    this->printNormalizedResiduals = true;
    this->transient = false;
    this->timeDiscretizationOrder=1;
    this->absTolerance=1e-8;
    this->relTolerance=1-4;
    this->showResidual=5;
    this->maxLevels=2;
    this->AgglomerationMethod="FaceArea";
    this->preSweeps=0;
    this->postSweeps=2;
    this->relFactor=.1;
    this->withNormal=false;
  }
  
  bool printNormalizedResiduals;
  bool transient;
  int timeDiscretizationOrder;
  double absTolerance;
  double relTolerance;
  int showResidual;
  int maxLevels;
  string AgglomerationMethod;
  int preSweeps;
  int postSweeps;
  double relFactor;
  bool withNormal;

};

#endif
