%{
#include "Mesh.h"
%}

%include "std_vector.i"

using namespace std;

class Mesh
{
private:
  Mesh(const int i);
};

typedef vector<Mesh*> MeshList;


%template(MeshList) vector<Mesh*>;
