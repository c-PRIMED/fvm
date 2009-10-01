#include "FloatVarDict.h"
#include "AMG.h"

template<class T>
struct ElectronicBC : public FloatVarDict<T>
{
  ElectronicBC()
  {
      this->defineVar("specifiedPotential",T(300.0));
      this->defineVar("specifiedPotentialFlux",T(0.0));
      this->defineVar("specifiedCharge",T(300.0));
      this->defineVar("specifiedChargeFlux",T(0.0));
      this->defineVar("timeStep",T(0.1));
  
     
  }
  string bcType;
};


template<class T>
struct ElectronicModelOptions : public FloatVarDict<T>
{
  ElectronicModelOptions()
  {
    this->defineVar("initialCharge",T(300.0));
    this->defineVar("initialPotential",T(300.0));
    this->defineVar("initialTunnelingCharge", T(1.0));
    this->electrostaticsTolerance=1e-8;
    this->chargetransportTolerance=1e-8;
    this->electrostaticsLinearSolver = 0;
    this->chargetransportLinearSolver = 0;
    this->transient = false;
    this->ibm = false;
    this->electrostatics = true;
    this->chargetransport = false;
    this->tunneling = false;
    this->emission = false;
    this->capture = false;
    this->drift = false;
    this->diffusion = false;
  }
  bool printNormalizedResiduals;

  double electrostaticsTolerance;
  double chargetransportTolerance;
  double tunnelingtransportTolerance;
  
  bool ibm;
  bool transient;
  bool tunneling;
  bool electrostatics;
  bool chargetransport;
  bool emission;
  bool capture;
  bool drift;
  bool diffusion;
  

  int timeDiscretizationOrder;
  LinearSolver *electrostaticsLinearSolver;
  LinearSolver *chargetransportLinearSolver;

#ifndef SWIG
  LinearSolver& getElectroStaticsLinearSolver()
  {
    if (this->electrostaticsLinearSolver == 0)
    {
        LinearSolver* ls(new  AMG());
        ls->relativeTolerance = 1e-1;
        ls->nMaxIterations = 20;
        ls->verbosity=0;
        this->electrostaticsLinearSolver = ls;
    }
    return *this->electrostaticsLinearSolver ;
  }

  LinearSolver& getChargeTransportLinearSolver()
  {
    if (this->chargetransportLinearSolver == 0)
    {
        LinearSolver* ls(new  AMG());
        ls->relativeTolerance = 1e-1;
        ls->nMaxIterations = 20;
        ls->verbosity=0;
        this->chargetransportLinearSolver = ls;
    }
    return *this->chargetransportLinearSolver ;
  }
#endif
};

