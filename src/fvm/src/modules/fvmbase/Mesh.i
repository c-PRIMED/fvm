%{
#include "Mesh.h"
  %}

%include "std_vector.i"

using namespace std;

struct FaceGroup
{
  const int id;
  const string groupType;
private:
  FaceGroup();
};

class Mesh
{
public:
  int getFaceGroupCount() const;
  //  const FaceGroup&  getFaceGroup(const int i) const;
  
private:
  Mesh(const int i);
};

typedef vector<Mesh*> MeshList;


%template(MeshList) vector<Mesh*>;
