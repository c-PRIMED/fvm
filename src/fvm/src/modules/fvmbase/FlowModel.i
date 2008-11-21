%{
#include "FlowModel.h"
#include "AMG.h"
%}

%include "atype.i"
%include "FloatVarDict.i"
%include "GeomFields.h"
%include "FlowFields.h"
%include "LinearSolver.i"
%include "Vector.i"

using namespace std;
using namespace boost;

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
  double momentumTolerance;
  double continuityTolerance;
  int timeDiscretizationOrder;
  LinearSolver *momentumLinearSolver;
  LinearSolver *pressureLinearSolver;
  LinearSolver *coupledLinearSolver;
}; 

%template(Vector3) Vector<ATYPE_STR,3>;

%template(FlowBCA) FlowBC<ATYPE_STR>;
%template(FlowBCList) std::vector<FlowBC<ATYPE_STR>* >;
%template(FlowBCsMap) std::map<int,FlowBC<ATYPE_STR>* >;

%template(FlowVCA) FlowVC<ATYPE_STR>;
%template(FlowVCList) std::vector<FlowVC<ATYPE_STR>* >;
%template(FlowVCsMap) std::map<int,FlowVC<ATYPE_STR>* >;

%template(FlowModelOptionsA) FlowModelOptions<ATYPE_STR>;


%include "Model.i"

%include "FlowModel.h"


%template(FlowModelA) FlowModel<ATYPE_STR>;

