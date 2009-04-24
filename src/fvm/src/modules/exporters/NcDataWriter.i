%{
#include "NcDataWriter.h"
%}

%include "std_vector.i"
%include "std_string.i"

using namespace std;

class NcDataWriter {

public :

    NcDataWriter(const MeshList& meshes, const string& fname);
    void  record();
};

