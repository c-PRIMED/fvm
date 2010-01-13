%{
#include "ThermalModel.h"
%}

using namespace std;

%include "FloatVarDict.i"

template <class T>
struct ThermalBC : public FloatVarDict<T>
{
  string bcType;
}; 

template <class T>
struct ThermalModelOptions : public FloatVarDict<T>
{
  double relativeTolerance;
  double absoluteTolerance;
  bool useCentralDifference;
  LinearSolver *linearSolver;
}; 


%template(ThermalBCA) ThermalBC<ATYPE_STR>;
%template(ThermalBCList) std::vector<ThermalBC<ATYPE_STR>* >;
%template(ThermalBCsMap) std::map<int,ThermalBC<ATYPE_STR>* >;

%template(ThermalModelOptionsA) ThermalModelOptions<ATYPE_STR>;


%import "Model.i"

%include "ThermalModel.h"


%template(ThermalModelA) ThermalModel<ATYPE_STR>;

