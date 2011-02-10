#ifndef _KINETICBC_H_
#define _KINETICBC_H_

#include "misc.h"
#include "FloatVarDict.h"
#include "AMG.h"

template<class T>
struct KineticBC : public FloatVarDict<T>
{
  KineticBC()
  {
      this->defineVar("specifiedXVelocity",T(0.00));
      this->defineVar("specifiedYVelocity",T(0.0));
      this->defineVar("specifiedZVelocity",T(0.0));
      this->defineVar("specifiedPressure",T(1.0));
      this->defineVar("specifiedDensity",T(1.0));
      this->defineVar("specifiedTemperature",T(1.0));
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
    this->defineVar("timeStep",T(1E-3)); 

    this->defineVar("rho_init",T(9.28E-6));
    this->defineVar("T_init",T(273.15));
    this->defineVar("nonDimLength",T(1.0));

    this->defineVar("operatingPressure",T(101325.0));
    this->defineVar("operatingTemperature",T(300.0));
    this->defineVar("molecularWeight",T(40.0));

    this ->defineVar("Tmuref",T(273.15));
    this ->defineVar("mu_w",T(0.81));
    this ->defineVar("muref",T(2.117e-5)); //Argon
    
    
    this->defineVar("pi",T(3.1416));
    
    this->Tolerance=1e-4;
    this->printNormalizedResiduals = true;
    this->transient = true;
    //this->ESBGK_fgamma = false; 
    
    this->fgamma=2;
    this->Prandtl=2.0/3.0; 
    this->SpHeatRatio=5.0/3.0; 
    this->timeDiscretizationOrder=1;
    this->KineticLinearSolver = 0;
   
    this-> printCellNumber=0;
    this->defineVar("printDirectionNumber",545);
    
    this->NewtonsMethod_ktrial=50;
    this->relativeTolerance=1e-14;
    this->absoluteTolerance=1e-22; 

    this->BoltzmannConstant=1.38e-23;

    //used in Newton's Method for Equilibrium distribution function
    this->defineVar("ToleranceX",T(1e-8));
    this->defineVar("ToleranceF",T(1e-16));
  }
  
  bool printNormalizedResiduals;
  double Tolerance;
  bool transient;
  //bool ESBGK_fgamma;
  double BoltzmannConstant;

  double relativeTolerance;
  double absoluteTolerance;
  
  //double ToleranceX;
  //double ToleranceF;
  int NewtonsMethod_ktrial;
  int timeDiscretizationOrder;
  LinearSolver *KineticLinearSolver;
  int printCellNumber;
  //int printDirectionNumber;
  int fgamma;
  
  double Prandtl;
  double SpHeatRatio;

#ifndef SWIG
  LinearSolver& getKineticLinearSolver()
  {
    if (this->KineticLinearSolver == 0)
    {
        LinearSolver* ls(new AMG());
        ls->relativeTolerance = 1e-5;
        ls->nMaxIterations = 20;
        ls->verbosity=0;
        this->KineticLinearSolver = ls;
    }
    return *this->KineticLinearSolver ;
  }
#endif

};

#endif
