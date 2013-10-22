%{
#include "Mesh.h"
  %}

//%include "Vector.i"

%include "std_vector.i"
%include "std_map.i"
%include "std_set.i"


using namespace std;


class StorageSite
{
  typedef map< const StorageSite*, boost::shared_ptr< Array<int> > > ScatterMap;
  typedef map< const StorageSite*, boost::shared_ptr< Array<int> > > GatherMap;
  typedef map< const StorageSite*, boost::shared_ptr< Array<int> > > CommonMap;
  typedef map<  const StorageSite*, map<int,int> >   ScatterIndex;  
public:
  int getCount() const;
  int getCountLevel1() const;
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

  %extend{
    std::map<int,int>   getGlobalToLocal( const StorageSite& other ){
         StorageSite::ScatterIndex scatterIndex = self->getScatterIndex();
         return scatterIndex.find(&other)->second;
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
  Mesh* extrude(int nz, double zmax, bool boundaryOnly=false);

  Mesh* createShell(const int fgId, Mesh& otherMesh, const int otherFgId);
  bool isShell() const;
  Mesh* createDoubleShell(const int fgId, Mesh& otherMesh, const int otherFgId, const bool connectedShell);
  bool isDoubleShell() const;
  bool isConnectedShell() const;
  const StorageSite& getParentFaceGroupSite() const;
  const StorageSite& getOtherFaceGroupSite() const;
  int getParentMeshID() const;
  int getOtherMeshID() const;
  
  // const GhostCellSiteMap& getGhostCellSiteMap() const;
  // const StorageSite* getGhostCellSite( int id );
  %extend{
    boost::shared_ptr<ArrayBase> getLocalToGlobalNodes()
    {
      boost::shared_ptr<ArrayBase> local = self->getLocalToGlobalNodesPtr();
      return self->getLocalToGlobalNodesPtr();

    }
  }
  const set<int>&  getBoundaryNodesSet() const;
  void setCommonFacesMap( const Mesh& bMesh );
  const map<int,int>& getCommonFacesMap() const;
  boost::shared_ptr<ArrayBase> getNodeCoordinatesPtr();   
  void setNodeRepeationArrayCoupling(const Mesh& bMesh);
  boost::shared_ptr< ArrayBase > getUpdatedNodesCoordCoupling(const GeomFields& geomField, const Mesh& bMesh);
  boost::shared_ptr< ArrayBase >   getLocalToGlobalPtr(){ return _localToGlobal;} 

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
    std::vector<FaceGroup*> getBoundaryFaceGroups()
    {
      std::vector<FaceGroup*> v;
      const FaceGroupList& fgs = self->getBoundaryFaceGroups();
      foreach(FaceGroupPtr fg, fgs)
        v.push_back(fg.get());
      return v;
    } 
	std::vector<FaceGroup*> getInterfaceGroups()
    {
      std::vector<FaceGroup*> v;
      const FaceGroupList& fgs = self->getInterfaceGroups();
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
 
  const CRConnectivity& getFaceCells2() const;
  void CRConnectivityPrint( const CRConnectivity& conn, int procID, const string& name );
  void CRConnectivityPrintFile(const CRConnectivity& conn, const string& name, const int procID) const;
  void InterfaceToBoundary();

private:
  Mesh(const int i);
};

typedef std::vector<Mesh*> MeshList;


%template(MeshList) std::vector<Mesh*>;
%template(MapInt)  std::map<int,int>;
