#include "FloatVarDict.h"

template<class T>
struct ThermalBC : public FloatVarDict<T>
{
  ThermalBC()
  {
      this->defineVar("specifiedTemperature",T(300.0));
      this->defineVar("specifiedHeatFlux",T(0.0));
      
  }
  string bcType;
};


template<class T>
struct ThermalModelOptions : public FloatVarDict<T>
{
  ThermalModelOptions()
  {
    this->defineVar("initialTemperature",T(300.0));
    this->relativeTolerance=1e-8;
    this->absoluteTolerance=1e-16;
  }
  double relativeTolerance;
  double absoluteTolerance;
};

