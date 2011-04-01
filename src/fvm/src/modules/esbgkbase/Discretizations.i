%{
#include "PhononConvectionDiscretization.h"
#include "PhononCollisionDiscretization.h"
  %}

%import "Discretization.h"
%include "PhononConvectionDiscretization.h"
%include "PhononCollisionDiscretization.h"

%template(PhononConDiscA) PhononConvectionDiscretization<ATYPE_STR,ATYPE_STR,ATYPE_STR>;
%template(PhononCollDiscA) PhononCollisionDiscretization<ATYPE_STR,ATYPE_STR,ATYPE_STR>;
