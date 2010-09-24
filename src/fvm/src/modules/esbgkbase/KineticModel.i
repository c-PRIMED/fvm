%{
/* include C++ header files necessary to compile the interface */

#include "Quadrature.h"
#include "MacroFields.h"
#include "DistFunctFields.h"
#include "KineticModel.h"
#include "KineticBC.h"
#include "TimeDerivativeDiscretization_Kmodel.h"
#include "CollisionTermDiscretization.h"
#include "ConvectionDiscretization_Kmodel.h"
#include "KineticBoundaryConditions.h"
#include "AMG.h"
  %}

using namespace std;
using namespace boost;
//%include "FloatVarDict.i"
%include "KineticModel.h"

%import "Mesh.i"

%import "DistFunctFields.h"
%import "MacroFields.h"
%import "Model.i"
%import "KineticBC.h"

 /*
template <class T>
struct KineticBC : public FloatVarDict<T>
{
  string bcType;
};

template <class T>
struct KineticVC : public FloatVarDict<T>
{
};

template <class T>
struct KineticModelOptions : public FloatVarDict<T>
{
  bool printNormalizedResiduals;
  double Tolerance;
  bool transient;

  double relativeTolerance;
  double absoluteTolerance;

  int timeDiscretizationOrder;
    LinearSolver *KineticLinearSolver;
    };*/


%template(KineticBCD) KineticBC<double>;
%template(KineticBCList) std::vector<KineticBC<double>*>;
%template(KineticBCsMap) std::map<int,KineticBC<double>*>;

%template(KineticVCD) KineticVC<double >;
%template(KineticVCList) std::vector<KineticVC<double>*>;
%template(KineticVCsMap) std::map<int,KineticVC<double>*>;
%template(KineticModelOptionsD) KineticModelOptions<double>;



%template(KineticModelD) KineticModel<double>;
%template(vdfField) DistFunctFields<double>;
