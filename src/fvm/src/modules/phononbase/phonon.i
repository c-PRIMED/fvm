%module phonon_atyped

%import "std_string.i"
%import "std_vector.i"
%import "std_except.i"
%import "std_map.i"
%import "shared_ptr.i"
%import "atype.i"
%import "Mesh.i"
%import "ArrayBase.i"
%import "Field.i"
%import "Vector.i"
 //%import "Discretization.h"

using namespace std;

%include "pmode.i"
%include "kvol.i"
%include "Kspace.i"
%include "PhononModel.i"
%include "PhononBoundary.i"
 //%include "PhononCollisionDiscretization.i"
 //%include "PhononConvectionDiscretization.i"