%{
/* include C++ header files necessary to compile the interface */

#include "Quadrature.h"
#include "MacroFields.h"
#include "DistFunctFields.h"
#include "KineticModel.h"
  %}

%import "Mesh.i"

%import "DistFunctFields.h"
%import "MacroFields.h"
%import "KineticModel.h"


%template(KineticModelD) KineticModel<double>;
%template(vdfField) DistFunctFields<double>;
