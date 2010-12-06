%{
#include "RosselandModel.h"
#include "AMG.h"
%}




using namespace std;
using namespace boost;

%include "FloatVarDict.i"


%include "RosselandModel.h"

%template(RosselandVCA) RosselandVC< ATYPE_STR >;
%template(RosselandVCList) std::vector<RosselandVC< ATYPE_STR >* >;
%template(RosselandVCsMap) std::map<int,RosselandVC< ATYPE_STR >* >;




%template(RosselandModelA) RosselandModel< ATYPE_STR >;

