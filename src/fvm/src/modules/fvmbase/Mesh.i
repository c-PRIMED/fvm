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
  int getSelfCount() const;
private:
  StorageSite();
};

struct FaceGroup
{
  const int id;
  const string groupType;
  const StorageSite site;
private:
  FaceGroup();
};

typedef std::vector<FaceGroup*> FaceGroupVector;


%template(FaceGroupVector) std::vector<FaceGroup*>;

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
  int getIBTypeForCell(const int c) const;
  
  void setIBTypeForCell(const int c, const int type);

  VecD3 getCellCoordinate(const int c) const;
  const StorageSite& getCells() const;
  const StorageSite& getNodes() const;
  const StorageSite& getFaces() const;

  const StorageSite& getIBFaces() const;
  
  const FaceGroup& getInteriorFaceGroup() const;
  
  int getFaceGroupCount() const;
  int getBoundaryGroupCount() const;
  int getInterfaceGroupCount() const;
  
  ArrayBase* getNodeCoordinates();

  const CRConnectivity& getConnectivity(const StorageSite& from,
                                        const StorageSite& to) const;
  
  %extend
  {
    int getNumberOfCells() {return self->getCells().getCount();}
    ArrayBase* getCoordinates()
    {
        cout << "mesh ext" << endl;
        return self->getNodeCoordinates();
    }

    std::vector<FaceGroup*> getBoundaryGroups()
    {
      std::vector<FaceGroup*> v;
      const FaceGroupList& fgs = self->getBoundaryFaceGroups();
      foreach(FaceGroupPtr fg, fgs)
        v.push_back(fg.get());
      return v;
    }
  }
  //  const FaceGroup&  getFaceGroup(const int i) const;
  
private:
  Mesh(const int i);
};

typedef std::vector<Mesh*> MeshList;


%template(MeshList) std::vector<Mesh*>;
