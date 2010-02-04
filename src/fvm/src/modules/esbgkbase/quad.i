%module quad
%{
/* include C++ header files necessary to compile the interface */
#include "quadrature.h"
%}

%import "quadrature.h"
%template(Quadrature) Quadrature<double>;
