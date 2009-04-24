%module NcReader
%{
#include "NcDataReader.h"
%}

%include "std_string.i"
%include "std_vector.i"


using namespace std;
class NcDataReader {

public :
  
    NcDataReader( const string& fname );

    MeshList    getMeshList();

  void createMappers(MeshList& globalMeshList);
};
