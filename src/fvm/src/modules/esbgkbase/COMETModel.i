%{
/* include C++ header files necessary to compile the interface */

#include "Quadrature.h"
#include "MacroFields.h"
#include "DistFunctFields.h"
#include "COMETModel.h"
#include "COMETBC.h"
#include "COMETBoundaryConditions.h"
#include "AMG.h"
#include "COMETESBGKDiscretizer.h"
  %}

using namespace std;
using namespace boost;

%include "atype.i"
%import "Mesh.i"
%include "FloatVarDict.i"
%include "Model.i"

%include "COMETModel.h"


%import "DistFunctFields.h"
%import "MacroFields.h"
%include "COMETBC.h"
%include "COMETESBGKDiscretizer.h"


%template(COMETBCD) COMETBC<double>;
%template(COMETBCList) std::vector<COMETBC<double>*>;
%template(COMETBCsMap) std::map<int,COMETBC<double>*>;

%template(COMETVCD) COMETVC<double >;
%template(COMETVCList) std::vector<COMETVC<double>*>;
%template(COMETVCsMap) std::map<int,COMETVC<double>*>;
%template(COMETModelOptionsD) COMETModelOptions<double>;



%template(COMETModelD) COMETModel<double>;
%template(vdfField) DistFunctFields<double>;
