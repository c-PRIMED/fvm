%{
#include "PhononModel.h"
  %}

%include "Model.i"
%include "PhononModel.h"

%template(PhononModelA) PhononModel< ATYPE_STR >;
