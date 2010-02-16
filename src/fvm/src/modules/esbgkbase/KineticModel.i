%module KineticModel
%{
/* include C++ header files necessary to compile the interface */

#include "quadrature.h"
#include "MacroParameters.h"
#include "DistFunctFields.h"
#include "KineticModel.h"
  %}

%import "quadrature.h"
%import "DistFunctFields.h"
%import "MacroParameters.h"
%import "KineticModel.h"

%template(Quadrature) Quadrature<double>;
%template(KineticModel) KineticModel<double>;
%template(DistFunctFields) DistFunctFields<double>;
