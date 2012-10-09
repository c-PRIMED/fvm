// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _ELECTRICBC_H_
#define _ELECTRICBC_H_


#include "FloatVarDict.h"
#include "AMG.h"

template<class T>
struct ElectricBC : public FloatVarDict<T>
{
  ElectricBC()
  {
      this->defineVar("specifiedPotential",T(300.0));
      this->defineVar("specifiedPotentialFlux",T(0.0));
      this->defineVar("specifiedXElecField",T(0.0));
      this->defineVar("specifiedYElecField",T(0.0));
      this->defineVar("specifiedZElecField",T(0.0));
      this->defineVar("specifiedCharge",T(0.0));
      this->defineVar("specifiedChargeFlux",T(0.0));
      this->defineVar("timeStep",T(1.0)); 
     
  }
  string bcType;
};

template<class T>
struct ElectricVC : public FloatVarDict<T>
{
  ElectricVC()
  {
      this->defineVar("dielectric_constant",T(7.9));
  }
  string vcType;
};

template<class T>
struct ElectricModelConstants : public FloatVarDict<T>
{
  ElectricModelConstants()
  {
    this->defineVar("dielectric_ionization", T(3.0));              //in EV
    this->defineVar("dielectric_bandgap", T (5.0));                 // in EV
    this->defineVar("optical_dielectric_constant", T(4.0));
    this->defineVar("dielectric_thickness", T(2.5e-7));          //in meter    
    this->defineVar("dielectric_constant", T(7.9));
    this->defineVar("electron_capture_cross", T(1e-17));         //in m^-2

    this->defineVar("membrane_workfunction", T (5.0));              //in EV

    this->defineVar("substrate_workfunction", T (5.0));             //in EV

    this->defineVar("membrane_voltage", T (0.0));

    this->defineVar("substrate_voltage", T(0.0));

    this->defineVar("OP_temperature", T(300.0));                   //in K
    this->defineVar("electron_effmass", T(0.5));
    this->defineVar("poole_frenkel_emission_frequency", T(1.0e+12));  //in s^-1
    this->defineVar("electron_mobility", T(50e4));                     // m^2 / Vs
    this->defineVar("electron_saturation_velocity", T(1e9));	      // m/s
    this->defineVar("voltage", T(100));                           // V

    this->defineVar("substrate_id", int(5));               //mesh id of substrate
    this->defineVar("membrane_id", int(4));                //mesh id of membrane
    this->defineVar("nLevel", int(0));                   //number of grid in normal direction
    this->defineVar("normal_direction", int(2));          //normal direction index; values: 0 1 2
    this->defineVar("nTrap", int(1));   
   
  }

 
  vector<T> electron_trapdensity;
  vector<T> electron_trapdepth;

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
    this->defineVar("timeStep",T(0.1));
    this->defineVar("Interface_A_coeff",T(1.0));
    this->defineVar("Interface_B_coeff",T(0.0));
    this->defineVar("ButlerVolmerRRConstant",T(5.0e-7));
    this->defineVar("ButlerVolmerAnodeShellMeshID", int(-1));
    this->defineVar("ButlerVolmerCathodeShellMeshID", int(-1));
    this->defineVar("BatteryElectrolyteMeshID", int(-1));
    this->ButlerVolmer = false;
    this->electrostaticsTolerance=1e-8;
    this->chargetransportTolerance=1e-8;
    this->electrostaticsLinearSolver = 0;
    this->chargetransportLinearSolver = 0;
    this->timeDiscretizationOrder = 1;
    this->transient_enable = true;
    this->ibm_enable = false;
    this->electrostatics_enable = true;
    this->chargetransport_enable = true;
    this->tunneling_enable = false;
    this->emission_enable = false;
    this->capture_enable = false;
    this->injection_enable = false;
    this->drift_enable = false;
    this->diffusion_enable = false;
    this->trapbandtunneling_enable = false;
    this->printNormalizedResiduals = true;
  }
  bool printNormalizedResiduals;

  double electrostaticsTolerance;
  double chargetransportTolerance;
  double tunnelingtransportTolerance;
  
  bool ibm_enable;
  bool transient_enable;
  bool tunneling_enable;
  bool emission_enable;
  bool electrostatics_enable;
  bool chargetransport_enable; 
  bool capture_enable;
  bool injection_enable;
  bool drift_enable;
  bool diffusion_enable;  
  bool trapbandtunneling_enable;
  bool ButlerVolmer;

  int timeDiscretizationOrder;
  LinearSolver *electrostaticsLinearSolver;
  LinearSolver *chargetransportLinearSolver;

#ifndef SWIG
  LinearSolver& getElectroStaticsLinearSolver()
  {
    if (this->electrostaticsLinearSolver == 0)
    {
        LinearSolver* ls(new  AMG());
        ls->relativeTolerance = 1e-3;
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
        ls->relativeTolerance = 1e-3;
        ls->nMaxIterations = 20;
        ls->verbosity=0;
        this->chargetransportLinearSolver = ls;
    }
    return *this->chargetransportLinearSolver ;
  }
#endif
};




#endif

