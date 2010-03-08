%{
/* include C++ header files necessary to compile the interface */

#include "quadrature.h"
#include "MacroParameters.h"
#include "DistFunctFields.h"
#include "KineticModel.h"
  %}

%import "Mesh.i"

%import "DistFunctFields.h"
%import "MacroParameters.h"
%import "KineticModel.h"


%template(KineticModelD) KineticModel<double>;
%template(vdfField) DistFunctFields<double>;
