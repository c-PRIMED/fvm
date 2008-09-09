%module importers
%{
#include "FluentReader.h"
%}

%include "std_string.i"

%include "Mesh.i"

using namespace std;

class FluentReader
{
public:
  FluentReader(const string& fileName);
  void readMesh();
  int getNumCells();
  MeshList getMeshList();
};
