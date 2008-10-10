#include "FloatVarDict.h"

template<class T>
struct FlowBC : public FloatVarDict<T>
{
  FlowBC()
  {
      this->defineVar("specifiedXVelocity",T(0.0));
      this->defineVar("specifiedYVelocity",T(0.0));
      this->defineVar("specifiedZVelocity",T(0.0));
      this->defineVar("specifiedPressure",T(0.0));
  }
  string bcType;
};

template<class T>
struct FlowVC : public FloatVarDict<T>
{
  FlowVC()
  {
      this->defineVar("viscosity",T(1e-3));
      this->defineVar("density",T(1.0));
  }
  string vcType;
};


template<class T>
struct FlowModelOptions : public FloatVarDict<T>
{
  FlowModelOptions()
  {
    this->defineVar("initialXVelocity",T(0.0));
    this->defineVar("initialYVelocity",T(0.0));
    this->defineVar("initialZVelocity",T(0.0));
    this->defineVar("initialPressure",T(0.0));
    this->defineVar("momentumURF",T(0.7));
    this->defineVar("pressureURF",T(0.3));
    this->defineVar("momentumTolerance",T(1e-3));
    this->defineVar("continuityTolerance",T(1e-3));
    this->printNormalizedResiduals = true;
  }
  bool printNormalizedResiduals;
};

