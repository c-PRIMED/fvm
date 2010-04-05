%{
/* include C++ header files necessary to compile the interface */
#include "Quadrature.h"
%}

%import "Quadrature.h"
%template(QuadratureD) Quadrature<double>;
