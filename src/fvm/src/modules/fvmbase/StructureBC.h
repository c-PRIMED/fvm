#include "misc.h"
#include "FloatVarDict.h"
#include "TractionVal.h"
#include "AMG.h"

template<class T>
struct StructureBC : public FloatVarDict<T>
{
  StructureBC()
  {
      this->defineVar("specifiedXDeformation",T(0.0));
      this->defineVar("specifiedYDeformation",T(0.0));
      this->defineVar("specifiedZDeformation",T(0.0));
      this->defineVar("specifiedXXTraction",T(0.0));
      this->defineVar("specifiedXYTraction",T(0.0));
      this->defineVar("specifiedXZTraction",T(0.0));
      this->defineVar("specifiedYXTraction",T(0.0));
      this->defineVar("specifiedYYTraction",T(0.0));
      this->defineVar("specifiedYZTraction",T(0.0));
      this->defineVar("specifiedZXTraction",T(0.0));
      this->defineVar("specifiedZYTraction",T(0.0));
      this->defineVar("specifiedZZTraction",T(0.0));
      this->defineVar("specifiedXForce",T(0.0));
      this->defineVar("specifiedYForce",T(0.0));
      this->defineVar("specifiedZForce",T(0.0));
  }
  string bcType;
};

template<class T>
struct StructureVC : public FloatVarDict<T>
{
  StructureVC()
  {
      this->defineVar("eta",T(1.0));
      this->defineVar("eta1",T(1.0));
      this->defineVar("density",T(1.0));
  }
  string vcType;
};


template<class T>
struct StructureModelOptions : public FloatVarDict<T>
{
  StructureModelOptions()
  {
    this->defineVar("initialXDeformation",T(0.0));
    this->defineVar("initialYDeformation",T(0.0));
    this->defineVar("initialZDeformation",T(0.0));
    this->defineVar("deformationURF",T(0.7));
    this->defineVar("timeStep",T(0.1));
    this->defineVar("operatingTemperature",T(300.0));

    this->deformationTolerance=1e-4;
    this->printNormalizedResiduals = true;
    this->transient = false;
    this->deformationLinearSolver = 0;
    this->coupledLinearSolver = 0;

    this->incompressible = true;
  }
  
  bool printNormalizedResiduals;
  double deformationTolerance;
  bool transient;
  LinearSolver *deformationLinearSolver;
  LinearSolver *coupledLinearSolver;

  bool incompressible;
#ifndef SWIG
  LinearSolver& getDeformationLinearSolver()
  {
    if (this->deformationLinearSolver == 0)
    {
        LinearSolver* ls(new AMG());
        ls->relativeTolerance = 1e-1;
        ls->nMaxIterations = 20;
        ls->verbosity=0;
        this->deformationLinearSolver = ls;
    }
    return *this->deformationLinearSolver ;
  }
#endif
};

