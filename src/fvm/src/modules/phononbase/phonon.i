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

using namespace std;

%include "Kspace.i"
%include "PhononModel.i"
%include "PhononBoundary.i"
%include "COMETModel.i"
%include "MatrixJML.i"
%include "ArrowHeadMatrix.i"
%include "SquareMatrix.i"
