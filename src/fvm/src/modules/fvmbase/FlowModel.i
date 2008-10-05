%{
#include "FlowModel.h"
%}

%include "atype.i"
%include "FloatVarDict.i"
%include "GeomFields.h"
%include "FlowFields.h"

using namespace std;

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
}; 


%template(FlowBCA) FlowBC<ATYPE_STR>;
%template(FlowBCList) std::vector<FlowBC<ATYPE_STR>* >;
%template(FlowBCsMap) std::map<int,FlowBC<ATYPE_STR>* >;

%template(FlowVCA) FlowVC<ATYPE_STR>;
%template(FlowVCList) std::vector<FlowVC<ATYPE_STR>* >;
%template(FlowVCsMap) std::map<int,FlowVC<ATYPE_STR>* >;

%template(FlowModelOptionsA) FlowModelOptions<ATYPE_STR>;


%include "Model.i"
%include "AMG.i"

%include "FlowModel.h"


%template(FlowModelA) FlowModel<ATYPE_STR>;

