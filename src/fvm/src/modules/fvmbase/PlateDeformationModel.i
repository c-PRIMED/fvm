%{
#include "PlateDeformationModel.h"
%}


using namespace std;
using namespace boost;

%include "FloatVarDict.i"


//%template(Vector3) Vector<ATYPE_STR,3>;

%include "PlateDeformationModel.h"


%template(PlateDeformationModelA) PlateDeformationModel< ATYPE_STR >;

