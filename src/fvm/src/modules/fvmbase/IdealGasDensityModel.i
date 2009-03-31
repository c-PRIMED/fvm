%{
#include "IdealGasDensityModel.h"
#include "AMG.h"
%}

%include "atype.i"
%include "FloatVarDict.i"
%include "FlowFields.h"
%include "Vector.i"

using namespace std;
using namespace boost;

%include "Model.i"

%include "IdealGasDensityModel.h"

%template(IdealGasVCA) IdealGasVC<ATYPE_STR>;
%template(IdealGasVCList) std::vector<IdealGasVC<ATYPE_STR>* >;
%template(IdealGasVCsMap) std::map<int,IdealGasVC<ATYPE_STR>* >;




%template(IdealGasDensityModelA) IdealGasDensityModel<ATYPE_STR>;

