%module esbgk
%{
#include "Field.h"
%}

%include "std_string.i"
%include "std_vector.i"
%include "std_except.i"
%include "std_map.i"
%include "shared_ptr.i"
using namespace std;

%import "ArrayBase.i"
%import "Field.i"
%import "Mesh.i"
%import "LinearSolver.i"
%import "Vector.i"

%include "Quadrature.i"
%include "MacroFields.h"
%include "KineticModel.i"

