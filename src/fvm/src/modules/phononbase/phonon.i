%module phonon_atyped

%include "std_string.i"
%include "std_vector.i"
%include "std_except.i"
%include "std_map.i"
%include "shared_ptr.i"

%include "atype.i"


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
