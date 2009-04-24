%{
#include "FluentDataExporter.h"
#include "atype.h"
%}

%include "std_string.i"

%include "atype.i"
%import "Mesh.i"

using namespace std;

template<class T>
class FluentDataExporter
{
public:
  FluentDataExporter(FluentReader& reader,
                     const string fileName,
                     const bool binary,
                     const int atypeComponent);
  void init();
  void finish();
  void writeScalarField(const Field& field,
                        const int fluentFieldId);
  void writeVectorField(const Field& field,
                        const int fluentFieldId);
};

%template(FluentDataExporterA) FluentDataExporter<ATYPE_STR>;

%include "NcDataWriter.i"

