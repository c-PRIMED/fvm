%{
#include "PhononConvectionDiscretization.h"
  %}

%include "PhononConvectionDiscretization.h"

%template(PhononConDiscA) PhononConvectionDiscretization<ATYPE_STR,ATYPE_STR,ATYPE_STR>;
