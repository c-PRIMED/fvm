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

  %extend{
    boost::shared_ptr<ArrayBase> getCommonIndices(const StorageSite& other)
    {
      StorageSite::CommonMap::iterator pos = self->getCommonMap().find(&other);
      if (pos != self->getCommonMap().end())
      {
          return pos->second;
      }
      return boost::shared_ptr<ArrayBase>();
    }
  }
  
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
  Mesh(const int dimension, const ArrayBase&  faceNodesCoord ) 
  {
     typedef Vector<double,3> Vec3D;
     typedef Array<Vec3D> Vec3DArray;
     const Vec3DArray& points(dynamic_cast<const Vec3DArray&>(faceNodesCoord));
     return new Mesh( dimension, points );
  }

  Mesh(const int dimension,
       const int nCells,
       const ArrayBase&  nodesCoordB,
       const ArrayBase& faceCellIndicesB,
       const ArrayBase& faceNodeIndicesB,
       const ArrayBase& faceNodeCountB,
       const ArrayBase& faceGroupSizeB )
  {
     typedef Vector<double,3> Vec3D;
     typedef Array<Vec3D> Vec3DArray;
     typedef Array<int> IntArray;
     const Vec3DArray&
       nodesCoord(dynamic_cast<const Vec3DArray&>(nodesCoordB));
     
     const IntArray&
       faceCellIndices(dynamic_cast<const IntArray&>(faceCellIndicesB));
     const IntArray&
       faceNodeIndices(dynamic_cast<const IntArray&>(faceNodeIndicesB));
     const IntArray&
       faceNodeCount(dynamic_cast<const IntArray&>(faceNodeCountB));
     const IntArray&
       faceGroupSize(dynamic_cast<const IntArray&>(faceGroupSizeB));
      
     return new Mesh( dimension, nCells, nodesCoord, faceCellIndices,
                      faceNodeIndices, faceNodeCount, faceGroupSize);
      
      
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
  int getCellZoneID() const;
  
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
  void findCommonFaces(StorageSite& faces, StorageSite& otherFaces,
                       const GeomFields& geomFields);
  
  Mesh* extractBoundaryMesh();
  Mesh* extrude(int nz, double zmax);
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
    const FaceGroup* getFaceGroup(const int i)
    {
      const FaceGroupList& fgs = self->getAllFaceGroups();
      foreach(FaceGroupPtr fg, fgs)
        if (i == fg->id)
          return fg.get();
      return 0;
    }
  }
  //  const FaceGroup&  getFaceGroup(const int i) const;
  
private:
  Mesh(const int i);
};

typedef std::vector<Mesh*> MeshList;


%template(MeshList) std::vector<Mesh*>;
