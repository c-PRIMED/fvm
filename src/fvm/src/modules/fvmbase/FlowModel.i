%{
#include "FlowModel.h"
#include "AMG.h"
%}


using namespace std;
using namespace boost;

%include "FloatVarDict.i"

template <class T>
struct FlowBC : public FloatVarDict<T>
{
  string bcType;
}; 

template <class T>
struct FlowVC : public FloatVarDict<T>
{
}; 

template <class T>
struct FlowModelOptions : public FloatVarDict<T>
{
  bool printNormalizedResiduals;
  bool transient;
  bool correctVelocity;
  bool incompressible;
  double momentumTolerance;
  double continuityTolerance;
  int timeDiscretizationOrder;
  LinearSolver *momentumLinearSolver;
  LinearSolver *pressureLinearSolver;
  LinearSolver *coupledLinearSolver;
  double cmu;
  double vk;
  double emp;
  bool turbulent;
}; 

//%template(Vector3) Vector<ATYPE_STR,3>;

%template(FlowBCA) FlowBC< ATYPE_STR >;
%template(FlowBCList) std::vector<FlowBC< ATYPE_STR >* >;
%template(FlowBCsMap) std::map<int,FlowBC< ATYPE_STR >* >;

%template(FlowVCA) FlowVC< ATYPE_STR >;
%template(FlowVCList) std::vector<FlowVC< ATYPE_STR >* >;
%template(FlowVCsMap) std::map<int,FlowVC< ATYPE_STR >* >;

%template(FlowModelOptionsA) FlowModelOptions< ATYPE_STR >;



%include "FlowModel.h"


%template(FlowModelA) FlowModel< ATYPE_STR >;

