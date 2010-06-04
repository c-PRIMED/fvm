%{
#include "StructureDeformationModel.h"
%}


using namespace std;
using namespace boost;

%include "FloatVarDict.i"


//%template(Vector3) Vector<ATYPE_STR,3>;

%include "StructureDeformationModel.h"


%template(StructureDeformationModelA) StructureDeformationModel< ATYPE_STR >;

