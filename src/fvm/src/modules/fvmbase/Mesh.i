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
  typedef map< const StorageSite*, boost::shared_ptr< Array<int> > > CommonMap;

public:
  int getCount() const;
  int getSelfCount() const;
  const ScatterMap&  getScatterMap() const;
  const GatherMap&   getGatherMap()  const;
  const CommonMap&   getCommonMap()  const;
  
  void clearGatherScatterMaps();
  
private:
  StorageSite();
};

struct FaceGroup
{
  const int id;
  string groupType;
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
      IBTYPE_BOUNDARY,
      IBTYPE_SOLID,
      IBTYPE_REALBOUNDARY,
      IBTYPE_UNKNOWN
    };
%extend{
  Mesh(const int dimension, const int id, const ArrayBase&  faceNodesCoord ) 
  {
     typedef Vector<double,3> Vec3D;
     typedef Array<Vec3D> Vec3DArray;
     const Vec3DArray& points(dynamic_cast<const Vec3DArray&>(faceNodesCoord));
     return new Mesh( dimension, id, points );
  }
}

  //VecD3 getCellCoordinate(const int c) const;
  const StorageSite& getCells() const;
  const StorageSite& getNodes() const;
  const StorageSite& getFaces() const;

  const StorageSite& getIBFaces() const;
  
  const FaceGroup& getInteriorFaceGroup() const;
  
  int getFaceGroupCount() const;
  int getBoundaryGroupCount() const;
  int getInterfaceGroupCount() const;
  int getID() const;
  shared_ptr<Array<int> > createAndGetBNglobalToLocal() const;
  const ArrayBase& getBNglobalToLocal() const;
  const StorageSite& getBoundaryNodes()const;
  
  const ArrayBase& getNodeCoordinates() const;

  const CRConnectivity& getConnectivity(const StorageSite& from,
                                        const StorageSite& to) const;
					
  const CRConnectivity& getAllFaceNodes() const;
  const CRConnectivity& getAllFaceCells() const;
  const CRConnectivity& getCellNodes() const;
  
  const CRConnectivity& getFaceCells(const StorageSite& site) const;
  const CRConnectivity& getFaceNodes(const StorageSite& site) const;
  const CRConnectivity& getCellFaces() const;
  const CRConnectivity& getCellCells() const;
					
  void findCommonNodes(Mesh& otherMesh);
  
  // const GhostCellSiteMap& getGhostCellSiteMap() const;
  // const StorageSite* getGhostCellSite( int id );

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
    std::vector<FaceGroup*> getAllFaceGroups()
    {
      std::vector<FaceGroup*> v;
      const FaceGroupList& fgs = self->getAllFaceGroups();
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
