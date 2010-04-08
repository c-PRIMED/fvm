%{
#include "IdealGasDensityModel.h"
#include "AMG.h"
%}

using namespace std;
using namespace boost;

%include "FloatVarDict.i"

%include "IdealGasDensityModel.h"

%template(IdealGasVCA) IdealGasVC< ATYPE_STR >;
%template(IdealGasVCList) std::vector<IdealGasVC< ATYPE_STR >* >;
%template(IdealGasVCsMap) std::map<int,IdealGasVC< ATYPE_STR >* >;




%template(IdealGasDensityModelA) IdealGasDensityModel< ATYPE_STR >;

