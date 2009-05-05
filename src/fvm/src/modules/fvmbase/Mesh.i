%{
#include "Mesh.h"
  %}

//%include "Vector.i"

%include "std_vector.i"
%include "std_map.i"


using namespace std;


class StorageSite
{
  typedef map< const StorageSite*, boost::shared_ptr< Array<int> > > ScatterMap;
  typedef map< const StorageSite*, boost::shared_ptr< Array<int> > > GatherMap;

public:
  int getCount() const;
  int getSelfCount() const;
  const ScatterMap&  getScatterMap() const;
  const GatherMap&   getGatherMap()  const;
  
  
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


class Mesh
{
public:

  typedef map<int, boost::shared_ptr<StorageSite> > GhostCellSiteMap;

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
					
  const CRConnectivity& getAllFaceNodes() const;
  const CRConnectivity& getAllFaceCells() const;
  const CRConnectivity& getCellNodes() const;
  
  const CRConnectivity& getFaceCells(const StorageSite& site) const;
  const CRConnectivity& getFaceNodes(const StorageSite& site) const;
  const CRConnectivity& getCellFaces() const;
  const CRConnectivity& getCellCells() const;
					
  
  const GhostCellSiteMap& getGhostCellSiteMap() const;
  const StorageSite* getGhostCellSite( int id );

  %extend
  {
    int getNumberOfCells() {return self->getCells().getCount();}

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
