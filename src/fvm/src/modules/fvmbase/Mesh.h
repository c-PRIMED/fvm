// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _MESH_H_
#define _MESH_H_

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif


#include "Array.h"

#include "StorageSite.h"
#include "Vector.h"
#include "Field.h"
#include "FieldLabel.h"
#include <vector>
#include <map>
#include <set>

class CRConnectivity;
class GeomFields;

class MPM;

struct FaceGroup
{
  FaceGroup(const int count_,
            const int offset_,
            const StorageSite& parent_,
            const int id_,
            const string& groupType_) :
    site(count_,0,offset_,&parent_),
    id(id_),
    groupType(groupType_)
  {}
  
  StorageSite site;
  const int id;
  string groupType;
};


typedef shared_ptr<FaceGroup> FaceGroupPtr; 
typedef vector<FaceGroupPtr> FaceGroupList;  

class Mesh 
{
public:

  typedef Vector<double,3> VecD3;
  typedef Array<int> IntArray;

  typedef vector<int>  vecList;
  typedef multimap<int,int> multiMap;
  typedef map<int,int> mapInt;

  typedef pair<const StorageSite*, const StorageSite*> SSPair;
  typedef map<SSPair,shared_ptr<CRConnectivity> > ConnectivityMap;
  typedef pair<int,int>    PartIDMeshIDPair; //other Partition ID, other MeshID (always)
  typedef map< PartIDMeshIDPair, shared_ptr<StorageSite> > GhostCellSiteMap;
  typedef pair<const StorageSite*, const StorageSite*> EntryIndex;
  typedef map<EntryIndex, shared_ptr<ArrayBase> > GhostArrayMap;

  typedef map<int,int> PeriodicFacePairs;
  
  enum
    {
      CELL_BAR2,
      CELL_TRI3,
      CELL_QUAD4,
      CELL_TETRA4,
      CELL_HEX8,
      CELL_PYRAMID,
      CELL_PRISM,
      CELL_TYPE_MAX
    } CellType;
  

  enum
    {
      IBTYPE_FLUID=-1,
      IBTYPE_BOUNDARY=-2,
      IBTYPE_SOLID=-3,
      IBTYPE_REALBOUNDARY=-4,
      IBTYPE_UNKNOWN=-5
    };
  
  Mesh(const int dimension);
  Mesh(const int dimension, const Array<VecD3>&  faceNodesCoord ); 
  Mesh(const int dimension,
       const int nCells,
       const Array<VecD3>&  nodesCoord,
       const Array<int>& faceCellIndices,
       const Array<int>& faceNodeIndices,
       const Array<int>& faceNodeCount,
       const Array<int>& faceGroupSize );
  
  ~Mesh();

  DEFINE_TYPENAME("Mesh");

  int getDimension() const {return _dimension;}
  int getID() const {return _id;}
  
  const StorageSite& getFaces()   const {return _faces;}
  const StorageSite& getCells()   const {return _cells;}
  const StorageSite& getNodes()   const {return _nodes;}
  const StorageSite& getIBFaces() const {return _ibFaces;}

  const StorageSite* getGhostCellSiteScatter(  const PartIDMeshIDPair& id ) const
  { return _ghostCellSiteScatterMap.find(id)->second.get(); }

  const GhostCellSiteMap& getGhostCellSiteScatterMap() const
  { return _ghostCellSiteScatterMap; }

  GhostCellSiteMap& getGhostCellSiteScatterMap() 
  { return _ghostCellSiteScatterMap; }

  const StorageSite* getGhostCellSiteGather( const PartIDMeshIDPair& id ) const
  { return _ghostCellSiteGatherMap.find(id)->second.get(); }

  GhostCellSiteMap& getGhostCellSiteGatherMap() 
  { return _ghostCellSiteGatherMap; }

  const GhostCellSiteMap& getGhostCellSiteGatherMap() const
  { return _ghostCellSiteGatherMap; }






  const StorageSite* getGhostCellSiteScatterLevel1(  const PartIDMeshIDPair& id ) const
  { return _ghostCellSiteScatterMapLevel1.find(id)->second.get(); }

  const GhostCellSiteMap& getGhostCellSiteScatterMapLevel1() const
  { return _ghostCellSiteScatterMapLevel1; }

  GhostCellSiteMap& getGhostCellSiteScatterMapLevel1() 
  { return _ghostCellSiteScatterMapLevel1; }

  const StorageSite* getGhostCellSiteGatherLevel1( const PartIDMeshIDPair& id ) const
  { return _ghostCellSiteGatherMapLevel1.find(id)->second.get(); }

  GhostCellSiteMap& getGhostCellSiteGatherMapLevel1() 
  { return _ghostCellSiteGatherMapLevel1; }

  const GhostCellSiteMap& getGhostCellSiteGatherMapLevel1() const
  { return _ghostCellSiteGatherMapLevel1; }



  StorageSite& getFaces() {return _faces;}
  StorageSite& getCells() {return _cells;}
  StorageSite& getNodes() {return _nodes;}
  StorageSite& getIBFaces() {return _ibFaces;}
 
  // this should only be used when we know that the connectivity
  // exists, for connectivities that are computed on demand the
  // specific functions below should be used
  
  const CRConnectivity& getConnectivity(const StorageSite& from,
                                        const StorageSite& to) const;

  const CRConnectivity& getAllFaceNodes() const;
  const CRConnectivity& getAllFaceCells() const;
  const CRConnectivity& getCellNodes() const;
  
  const CRConnectivity& getFaceCells(const StorageSite& site) const;
  const CRConnectivity& getFaceNodes(const StorageSite& site) const;
  const CRConnectivity& getCellFaces() const;
  const CRConnectivity& getCellCells() const;
  const CRConnectivity& getCellCells2() const;
  const CRConnectivity& getFaceCells2() const;

  CRConnectivity& getAllFaceCells();
  
  const FaceGroup& getInteriorFaceGroup() const {return *_interiorFaceGroup;}
  
  int getFaceGroupCount() const {return _faceGroups.size();}
  int getBoundaryGroupCount() const {return _boundaryGroups.size();}
  int getInterfaceGroupCount() const {return _interfaceGroups.size();}

  const FaceGroupList& getBoundaryFaceGroups() const
  {return _boundaryGroups;}

  const FaceGroupList& getInterfaceGroups() const
  {return _interfaceGroups;}
  
  const FaceGroupList& getAllFaceGroups() const
  {return _faceGroups;}
  
  const FaceGroup& getFaceGroup(const int fgId) const;
  
  const StorageSite& createInteriorFaceGroup(const int size);
  const StorageSite& createInterfaceGroup(const int size,const int offset, 
                                    const int id);
  const StorageSite& createBoundaryFaceGroup(const int size, const int offset, 
                                       const int id, const string& boundaryType);

  void setCoordinates(shared_ptr<Array<VecD3> > x) {_coordinates = x;}
  void setFaceNodes(shared_ptr<CRConnectivity> faceNodes);
  void setFaceCells(shared_ptr<CRConnectivity> faceCells);
  
  void setConnectivity(const StorageSite& rowSite, const StorageSite& colSite,
		       shared_ptr<CRConnectivity> conn);
  void eraseConnectivity(const StorageSite& rowSite,
                         const StorageSite& colSite) const;
		    
  
  shared_ptr<Array<int> > createAndGetBNglobalToLocal() const;
  const ArrayBase& getBNglobalToLocal() const;
  const StorageSite& getBoundaryNodes() const;  

  const Array<VecD3>& getNodeCoordinates() const {return *_coordinates;}
  Array<VecD3>& getNodeCoordinates() {return *_coordinates;}
  shared_ptr<ArrayBase> getNodeCoordinatesPtr() {return _coordinates;}

  void setNumOfAssembleMesh( int nmesh ){ _numOfAssembleMesh = nmesh; }

  //VecD3 getCellCoordinate(const int c) const;

  const Array<int>& getIBFaceList() const;

  Array<int>& getCellColors() { return *_cellColor;}
  const Array<int>& getCellColors() const { return *_cellColor;}

  Array<int>& getCellColorsOther() { return *_cellColorOther;}
  const Array<int>& getCellColorsOther() const { return *_cellColorOther;}

  shared_ptr< Array<int> > getLocalToGlobalNodesPtr()  { return _localToGlobalNodes;}
  const map<int,int> getGlobalToLocalNodes() const { return _globalToLocalNodes;}
  map<int,int> getGlobalToLocalNodes()    { return _globalToLocalNodes;}


  shared_ptr< ArrayBase >        getLocalToGlobalPtr(){ return _localToGlobal;}
  Array<int>&        getLocalToGlobal(){ return *_localToGlobal;}
  const Array<int>&  getLocalToGlobal() const { return *_localToGlobal;}

  map<int,int>&        getGlobalToLocal(){ return _globalToLocal;}
  const map<int,int>&  getGlobalToLocal() const { return _globalToLocal;}


  multiMap& getCellCellsGlobal() { return _cellCellsGlobal;}
  const multiMap& getCellCellsGlobal() const { return _cellCellsGlobal;}

  bool isMergedMesh() const { return _isAssembleMesh;}
  int  getNumOfAssembleMesh() const { return _numOfAssembleMesh;}

  const set<int>&  getBoundaryNodesSet() const { return _boundaryNodesSet;}
  set<int>&  getBoundaryNodesSet() { return _boundaryNodesSet;}
  const map<int,int>& getCommonFacesMap() const { return _commonFacesMap;}
  const map<int,int>& getCommonFacesMapOther() const { return _commonFacesMapOther;}

  
  void createScatterGatherCountsBuffer();
  void recvScatterGatherCountsBufferLocal();
  void syncCounts();
  
  void createScatterGatherIndicesBuffer();
  void recvScatterGatherIndicesBufferLocal();
  void syncIndices();
  
  const CRConnectivity& getCellCellsGhostExt() const { 
#ifdef FVM_PARALLEL
      int numProcs = MPI::COMM_WORLD.Get_size();
      if ( numProcs > 1 )
          return *_cellCellsGhostExt;
      else
#endif      
        return this->getCellCells();
  }

  const ArrayBase&      getSendCounts ( const EntryIndex& e ) const { return *_sendCounts[e];}
  const ArrayBase&      getSendIndices( const EntryIndex& e ) const { return *_sendIndices[e];}

  const ArrayBase&      getRecvCounts ( const EntryIndex& e ) const { return *_recvCounts[e];}
  const ArrayBase&      getRecvIndices( const EntryIndex& e ) const { return *_recvIndices[e];}
  
  
 
  void  createCellCellsGhostExt();
   
  void uniqueFaceCells();


  void setIBFaces(shared_ptr<Array<int> > faceList) {_ibFaceList = faceList;}

  void createGhostCellSiteScatter( const PartIDMeshIDPair& id, shared_ptr<StorageSite> site ); 
  void createGhostCellSiteGather ( const PartIDMeshIDPair& id, shared_ptr<StorageSite> site ); 
  void createGhostCellSiteScatterLevel1( const PartIDMeshIDPair& id, shared_ptr<StorageSite> site ); 
  void createGhostCellSiteGatherLevel1 ( const PartIDMeshIDPair& id, shared_ptr<StorageSite> site ); 

  void createCellColor();
  void createLocalGlobalArray();
  void createLocalToGlobalNodesArray();

  void                      setNodeRepeationArrayCoupling(const Mesh& bMesh);
  shared_ptr< ArrayBase >   getUpdatedNodesCoordCoupling(const GeomFields& geomField, const Mesh& bMesh);

  void setCommonFacesMap( const Mesh& bMesh );


  void findCommonNodes(Mesh& other);
  void findCommonFaces(StorageSite& faces, StorageSite& otherFaces,
                       const GeomFields& geomFields);
  bool COMETfindCommonFaces(StorageSite& faces, StorageSite& otherFaces,
			    const GeomFields& geomFields);
  
  Mesh* extractBoundaryMesh();
  Mesh* extrude(int nz, double zmax, bool boundaryOnly=false);

  Mesh* createShell(const int fgId, Mesh& otherMesh, const int otherFgId);
  Mesh* createDoubleShell(const int fgId, Mesh& otherMesh, const int otherFgId, const bool connectedShell);

  int getCellZoneID() const { return _cellZoneID;}
  void setCellZoneID(const int id) {_cellZoneID = id;}
  void setID(const int id) {_id = id;}
  
  bool isShell() const {return _isShell;}
  bool isDoubleShell() const {return _isDoubleShell;}
  bool isConnectedShell() const {return _isConnectedShell;}
  int getParentMeshID() const {return _parentMeshID;}
  int getOtherMeshID() const {return _otherMeshID;}

  const StorageSite& getParentFaceGroupSite() const
  {return *_parentFaceGroupSite;}

  const StorageSite& getOtherFaceGroupSite() const
  {return *_otherFaceGroupSite;}

  ConnectivityMap&  getConnectivityMap() {return _connectivityMap;}

  PeriodicFacePairs& getPeriodicFacePairs() { return _periodicFacePairs;}
  const PeriodicFacePairs& getPeriodicFacePairs() const { return _periodicFacePairs;}
  void CRConnectivityPrint(const CRConnectivity& conn, int procID, const string& name);
  void CRConnectivityPrintFile(const CRConnectivity& conn, const string& name, const int procID) const;
  void InterfaceToBoundary();

protected:
   
  
  const int _dimension;

  // used for persistence etc. Each mesh we create has a unique ID
  int _id;

  // used to get bcs from Fluent case.
  int _cellZoneID;
  
  StorageSite _cells;
  StorageSite _faces;
  StorageSite _nodes;

  StorageSite _ibFaces;
  mutable StorageSite* _boundaryNodes;
 
  shared_ptr<FaceGroup> _interiorFaceGroup;
  FaceGroupList _faceGroups;
  FaceGroupList _boundaryGroups;
  FaceGroupList _interfaceGroups;
  mutable ConnectivityMap _connectivityMap;
  shared_ptr<Array<VecD3> > _coordinates;
  mutable shared_ptr<Array<int> > _boundaryNodeGlobalToLocalPtr;

  mutable shared_ptr<Array<int> > _ibFaceList;

  shared_ptr< Array<int> > _cellColor;
  shared_ptr< Array<int> > _cellColorOther;
  int  _numOfAssembleMesh;
  bool _isAssembleMesh;


  GhostCellSiteMap   _ghostCellSiteScatterMap;
  GhostCellSiteMap   _ghostCellSiteGatherMap;

  GhostCellSiteMap   _ghostCellSiteScatterMapLevel1;
  GhostCellSiteMap   _ghostCellSiteGatherMapLevel1;
  
  mutable GhostArrayMap   _sendCounts;
  mutable GhostArrayMap   _recvCounts;
  mutable GhostArrayMap   _sendIndices;
  mutable GhostArrayMap   _recvIndices;

  shared_ptr<StorageSite> _cellSiteGhostExt;
  shared_ptr<CRConnectivity> _cellCellsGhostExt;

  
  map<int,int>                 _commonFacesMap;    //key=localFaceID, value = corresponding boundaryMehsFace
  map<int,int>                 _commonFacesMapOther;//key=boundaryMeshFace(sitting in this process), value=localFaceID
  set<int>                     _boundaryNodesSet; //store local number of boundary nodes  
  shared_ptr< Array<int>  >    _repeatNodes;
  shared_ptr< Array<int>  >    _localToGlobalNodes;
  shared_ptr< Array<int>  >    _localToGlobal;
  mutable map <int,int>        _globalToLocal;
  mutable map <int,int>        _globalToLocalNodes;
  multiMap             _cellCellsGlobal; //this hold cellCells information in global numbering, key is local,
                                        // values are global neighbouring and itself(global again)
 
  //mutable Array<int> *_cellTypes;
  //mutable Array<int> *_cellTypeCount;
  mutable shared_ptr<CRConnectivity> _cellCells2;
  mutable shared_ptr<CRConnectivity> _faceCells2;

  bool _isShell;
  bool _isDoubleShell;
  bool _isConnectedShell;

  // used if this is a shell mesh, points to the face group site that
  // can be used to obtain the area of the faces that makes up the
  // cells
  
  const StorageSite* _parentFaceGroupSite;
  const StorageSite* _otherFaceGroupSite;
  int _parentMeshID;
  int _otherMeshID;
  
  static int _lastID;

  PeriodicFacePairs _periodicFacePairs;

private:
  void createRowColSiteCRConn();
  void countCRConn();
  void addCRConn();
  int  getNumBounCells();
  int  get_request_size();

      
  
  
  
  
};

typedef vector<Mesh*> MeshList;

#endif
