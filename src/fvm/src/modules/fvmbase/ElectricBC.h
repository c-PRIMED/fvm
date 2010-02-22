#include "FloatVarDict.h"
#include "AMG.h"

template<class T>
struct ElectricBC : public FloatVarDict<T>
{
  ElectricBC()
  {
      this->defineVar("specifiedPotential",T(300.0));
      this->defineVar("specifiedXElecField",T(0.0));
      this->defineVar("specifiedYElecField",T(0.0));
      this->defineVar("specifiedZElecField",T(0.0));
      this->defineVar("specifiedCharge",T(0.0));
      this->defineVar("specifiedChargeFlux",T(0.0));
      this->defineVar("timeStep",T(0.1)); 
     
  }
  string bcType;
};

template<class T>
struct ElectricVC : public FloatVarDict<T>
{
  ElectricVC()
  {
      this->defineVar("dielectric_constant",T(1.0));
  }
  string vcType;
};

template<class T>
struct ElectricModelOptions : public FloatVarDict<T>
{
  ElectricModelOptions()
  {
    this->defineVar("initialCharge",T(0));
    this->defineVar("initialPotential",T(0.0));
    this->defineVar("initialTotalCharge", T(0.0));
    this->defineVar("initialTunnelingCharge", T(1.0));
    this->electrostaticsTolerance=1e-8;
    this->chargetransportTolerance=1e-8;
    this->electrostaticsLinearSolver = 0;
    this->chargetransportLinearSolver = 0;
    this->transient = false;
    this->ibm = false;
    this->electrostatics_enable = true;
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
  bool electrostatics_enable;
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

