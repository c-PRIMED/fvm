%{
#include "Mesh.h"
  %}

%include "Vector.i"

%include "std_vector.i"

using namespace std;

class StorageSite
{
public:
  int getCount() const;
private:
  StorageSite();
};

struct FaceGroup
{
  const int id;
  const string groupType;
private:
  FaceGroup();
};

%template(VecD3) Vector<double,3>;

class Mesh
{
public:
  typedef Vector<double,3> VecD3;

  enum
    {
      IBTYPE_FLUID,
      IBTYPE_SOLID,
      IBTYPE_BOUNDARY
    };

  int getFaceGroupCount() const;
  int getIBTypeForCell(const int c) const;
  
  void setIBTypeForCell(const int c, const int type);

  VecD3 getCellCoordinate(const int c) const;
  const StorageSite& getCells() const;
  ArrayBase* getNodeCoordinates();
  
  %extend
  {
    int getNumberOfCells() {return self->getCells().getCount();}
    ArrayBase* getCoordinates()
    {
        cout << "mesh ext" << endl;
        return self->getNodeCoordinates();
    }
  }
  //  const FaceGroup&  getFaceGroup(const int i) const;
  
private:
  Mesh(const int i);
};

typedef vector<Mesh*> MeshList;


%template(MeshList) vector<Mesh*>;
