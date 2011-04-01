%{
#include "pmode.h"
  %}

%import "Field.i"
%include "pmode.h"

%template(pmodeA) pmode< ATYPE_STR >;
