#include "misc.h"
#include "FloatVarDict.h"
#include "AMG.h"

template<class T>
struct KineticBC : public FloatVarDict<T>
{
  KineticBC()
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
struct KineticVC : public FloatVarDict<T>
{
  KineticVC()
  {
      this->defineVar("viscosity",T(1e-3));
      this->defineVar("density",T(1.0));
  }
  string vcType;
};


template<class T>
struct KineticModelOptions : public FloatVarDict<T>
{
  KineticModelOptions()
  {
    this->defineVar("initialXVelocity",T(0.0));
    this->defineVar("initialYVelocity",T(0.0));
    this->defineVar("initialZVelocity",T(0.0));
    this->defineVar("initialPressure",T(0.0));
    this->defineVar("momentumURF",T(0.7));
    this->defineVar("velocityURF",T(1.0));
    this->defineVar("pressureURF",T(0.3));
    this->defineVar("timeStep",T(0.0001)); 
    this->defineVar("nonDimLength",T(1.0));
    this->defineVar("operatingPressure",T(101325.0));
    this->defineVar("operatingTemperature",T(300.0));
    this->defineVar("molecularWeight",T(28.966));

    this->Tolerance=1e-4;
    this->printNormalizedResiduals = true;
    this->transient = false;
    this->timeDiscretizationOrder=1;
    this->KineticLinearSolver = 0;
     
    this->relativeTolerance=1e-8;
    this->absoluteTolerance=1e-16;
  }
  
  bool printNormalizedResiduals;
  double Tolerance;
  bool transient;
  
  double relativeTolerance;
  double absoluteTolerance;

  int timeDiscretizationOrder;
  LinearSolver *KineticLinearSolver;
 
#ifndef SWIG
  LinearSolver& getKineticLinearSolver()
  {
    if (this->KineticLinearSolver == 0)
    {
        LinearSolver* ls(new AMG());
        ls->relativeTolerance = 1e-1;
        ls->nMaxIterations = 20;
        ls->verbosity=0;
        this->KineticLinearSolver = ls;
    }
    return *this->KineticLinearSolver ;
  }
#endif

};

