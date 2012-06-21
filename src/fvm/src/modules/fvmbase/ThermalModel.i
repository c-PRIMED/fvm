%{
#include "ThermalModel.h"
%}

using namespace std;

%include "FloatVarDict.i"
%include "Model.i"

template <class T>
struct ThermalBC : public FloatVarDict<T>
{
  string bcType;
}; 

template <class T>
struct ThermalVC : public FloatVarDict<T>
{
  string vcType;
}; 

template <class T>
struct ThermalModelOptions : public FloatVarDict<T>
{
  double relativeTolerance;
  double absoluteTolerance;
  LinearSolver *linearSolver;
  bool useCentralDifference;
  bool transient;
}; 


%template(ThermalBCA) ThermalBC< ATYPE_STR >;
%template(ThermalVCA) ThermalVC< ATYPE_STR >;
%template(ThermalBCList) std::vector<ThermalBC< ATYPE_STR >* >;
%template(ThermalBCsMap) std::map<int,ThermalBC< ATYPE_STR >* >;
%template(ThermalVCList) std::vector<ThermalVC< ATYPE_STR >* >;
%template(ThermalVCsMap) std::map<int,ThermalVC< ATYPE_STR >* >;

%template(ThermalModelOptionsA) ThermalModelOptions< ATYPE_STR >;


%import "Model.i"

%include "ThermalModel.h"


%template(ThermalModelA) ThermalModel< ATYPE_STR >;

