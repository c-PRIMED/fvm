
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

%include "GeomFields.h"
%include "FlowFields.h"
%include "ThermalFields.h"

%include "atype.i"
%include "Model.i"

%include "MeshMetricsCalculator.i"

%include "ThermalModel.i"

%include "FlowModel.i"

%include "IdealGasDensityModel.i"

%include "MovingMeshModel.i"

#ifdef USING_ATYPE_TANGENT

typedef Vector<Tangent,3> VecTangent3;
%template(VecTangent3) Vector<Tangent,3>;

#endif
