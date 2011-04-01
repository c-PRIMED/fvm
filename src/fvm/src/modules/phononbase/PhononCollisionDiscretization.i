%{
#include "PhononCollisionDiscretization.h"
  %}

%include "PhononCollisionDiscretization.h"

%template(PhononCollDiscA) PhononCollisionDiscretization<ATYPE_STR,ATYPE_STR,ATYPE_STR>;
