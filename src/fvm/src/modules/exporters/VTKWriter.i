%{
#include "VTKWriter.h"
#include "atype.h"
%}

%include "std_string.i"

%include "atype.i"
%import "Mesh.i"
%import "GeomFields.h"

using namespace std;

template<class T>
class VTKWriter
{
public:
  VTKWriter(const GeomFields& geomFields,
            const MeshList& meshes,
            const string fileName,
            const string comment,
            const bool binary,
            const int atypeComponent,
            const bool surfaceOnly=false);
  void init();
  void finish();
  void writeScalarField(const Field& field,
                        const string label);
  void writeVectorField(const Field& field,
                        const string label);
};

%template(VTKWriterA) VTKWriter< ATYPE_STR >;


