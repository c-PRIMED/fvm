%{
#include "ThermalModel.h"
%}

%include "atype.i"
%include "FloatVarDict.i"
%include "GeomFields.h"
%include "ThermalFields.h"

using namespace std;

template <class T>
struct ThermalBC : public FloatVarDict<T>
{
  string bcType;
}; 

template <class T>
struct ThermalModelOptions : public FloatVarDict<T>
{
}; 


%template(ThermalBCA) ThermalBC<ATYPE_STR>;
%template(ThermalBCList) std::vector<ThermalBC<ATYPE_STR>* >;
%template(ThermalBCsMap) std::map<int,ThermalBC<ATYPE_STR>* >;

%template(ThermalModelOptionsA) ThermalModelOptions<ATYPE_STR>;


%include "Model.i"

%include "ThermalModel.h"


%template(ThermalModelA) ThermalModel<ATYPE_STR>;

