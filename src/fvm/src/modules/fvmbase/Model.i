%{
#include "Model.h"
%}

%include "std_vector.i"

%include "Mesh.i"

using namespace std;

class Model
{
public:
  Model(const MeshList& meshes);
  virtual void init()=0;
};

