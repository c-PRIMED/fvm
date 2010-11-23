#include <set>
#include "Mesh.h"
#include "StorageSite.h"
#include "CRConnectivity.h"
#include "Cell.h"
#include <cassert>
#include "KSearchTree.h"
#include "GeomFields.h"

#define epsilon 1e-6

int Mesh::_lastID = 0;

Mesh::Mesh(const int dimension):
  _dimension(dimension),
  _id(_lastID++),
  _cellZoneID(-1),
  _cells(0),
  _faces(0),
  _nodes(0),
  _ibFaces(0),
  _boundaryNodes(0),
  _interiorFaceGroup(),
  _faceGroups(),
  _boundaryGroups(),
  _interfaceGroups(),
  _connectivityMap(),
  _coordinates(),
  _boundaryNodeGlobalToLocalPtr(),
  _ibFaceList(),
  _numOfAssembleMesh(1),
  _isAssembleMesh(false)
{
  logCtor();
}

Mesh::Mesh( const int dimension,
            const Array<VecD3>&  faceNodesCoord ): 
  _dimension(dimension),
  _id(_lastID++),
  _cellZoneID(-1),
  _cells(0),
  _faces(0),
  _nodes(0),
  _ibFaces(0),
  _boundaryNodes(0),
  _interiorFaceGroup(),
  _faceGroups(),
  _boundaryGroups(),
  _interfaceGroups(),
  _connectivityMap(),
  _coordinates(),
  _boundaryNodeGlobalToLocalPtr(),
  _ibFaceList(),
  _numOfAssembleMesh(1),
  _isAssembleMesh(false)
{
  int faceNodeCount = _dimension == 2 ? 2 :4;
  int totNodes      = faceNodesCoord.getLength();

  // counting duplicate nodes as well for 3d
  int totFaces      = totNodes / faceNodeCount;

  //check if this is corect integer division
  assert( (faceNodeCount*totFaces) == totNodes );
  //set sites
  StorageSite& faceSite = getFaces();
  StorageSite& nodeSite = getNodes();
  faceSite.setCount( totFaces );
  nodeSite.setCount( totNodes );
  //interior face group (we have only one interface for this
  createBoundaryFaceGroup(totFaces,0,0,"wall");
  
  //setting coordinates
  shared_ptr< Array<VecD3> > coord( new Array< VecD3 > ( totNodes ) );
  *coord = faceNodesCoord;
  setCoordinates( coord );
  
  //faceNodes constructor
  shared_ptr<CRConnectivity> faceNodes( new CRConnectivity(faceSite,
                                                           nodeSite) );

  faceNodes->initCount();
  //addCount
  for ( int i = 0; i < totFaces; i++ )
    faceNodes->addCount(i, faceNodeCount); 
  //finish count
  faceNodes->finishCount();
  //add operation 
  int face     = 0;
  int nodeIndx = 0;
  for( int i = 0; i < totFaces; i++ )
  {
      for( int j =0; j < faceNodeCount; j++ )
      {
          faceNodes->add(face, nodeIndx++);
      }
      face++;
  }
  //finish add
  faceNodes->finishAdd();
  
  //setting faceNodes
  SSPair key(&faceSite,&nodeSite);
  _connectivityMap[key] = faceNodes;

  logCtor();
}

Mesh::Mesh( const int dimension,
            const int nCells,
            const Array<VecD3>&  nodesCoord,
            const Array<int>& faceCellIndices,
            const Array<int>& faceNodeIndices,
            const Array<int>& faceNodeCount,
            const Array<int>& faceGroupSize
            ): 
  _dimension(dimension),
  _id(_lastID++),
  _cellZoneID(-1),
  _cells(0),
  _faces(0),
  _nodes(0),
  _ibFaces(0),
  _boundaryNodes(0),
  _interiorFaceGroup(),
  _faceGroups(),
  _boundaryGroups(),
  _interfaceGroups(),
  _connectivityMap(),
  _coordinates(),
  _boundaryNodeGlobalToLocalPtr(),
  _ibFaceList(),
  _numOfAssembleMesh(1),
  _isAssembleMesh(false)
{
  int nFaces = faceNodeCount.getLength();
  int nNodes      = nodesCoord.getLength();

  StorageSite& faceSite = getFaces();
  StorageSite& nodeSite = getNodes();
  StorageSite& cellSite = getCells();
  faceSite.setCount( nFaces );
  nodeSite.setCount( nNodes );

  //interior face group (we have only one interface for this
  createInteriorFaceGroup(faceGroupSize[0]);
  
  int nFaceGroups = faceGroupSize.getLength();

  int offset = faceGroupSize[0];
  int nBoundaryFaces =0;
  for(int nfg=1; nfg<nFaceGroups; nfg++)
  {
      createBoundaryFaceGroup(faceGroupSize[nfg], offset, nfg, "wall");
      offset += faceGroupSize[nfg];
      nBoundaryFaces += faceGroupSize[nfg];
  }
  
  cellSite.setCount( nCells, nBoundaryFaces );
  
  //setting coordinates
  shared_ptr< Array<VecD3> > coord( new Array< VecD3 > ( nNodes ) );
  *coord = nodesCoord;
  setCoordinates( coord );
  
  
  //faceNodes constructor
  shared_ptr<CRConnectivity> faceNodes( new CRConnectivity(faceSite,
                                                           nodeSite) );

  faceNodes->initCount();

  shared_ptr<CRConnectivity> faceCells( new CRConnectivity(faceSite,
                                                           cellSite) );

  faceCells->initCount();

  
  //addCount
  for ( int f = 0; f < nFaces; f++ )
  {
      faceNodes->addCount(f, faceNodeCount[f]);
      faceCells->addCount(f, 2);
  }
  
  //finish count
  faceNodes->finishCount();
  faceCells->finishCount();
  
  //add operation 
  int nfn=0;
  int nfc=0;

  for( int f = 0; f < nFaces; f++ )
  {
      for( int j =0; j < faceNodeCount[f]; j++ )
      {
          faceNodes->add(f, faceNodeIndices[nfn++]);
      }
      faceCells->add(f,faceCellIndices[nfc++]);
      faceCells->add(f,faceCellIndices[nfc++]);
  }
  //finish add
  faceNodes->finishAdd();
  faceCells->finishAdd();
  
  //setting faceNodes
  SSPair key(&faceSite,&nodeSite);
  _connectivityMap[key] = faceNodes;
  
  SSPair key2(&faceSite,&cellSite);
  _connectivityMap[key2] = faceCells;

  logCtor();
}


Mesh::~Mesh()
{
  logDtor();
}




const StorageSite&
Mesh::createInteriorFaceGroup(const int size)
{
  _interiorFaceGroup = shared_ptr<FaceGroup>(new FaceGroup(size,0,_faces,0,"interior"));
  _faceGroups.push_back(_interiorFaceGroup);
  return _interiorFaceGroup->site;
}


const StorageSite&
Mesh::createInterfaceGroup(const int size, const int offset, const int id)
{
  shared_ptr<FaceGroup> fg(new FaceGroup(size,offset,_faces,id,"interface"));
  _faceGroups.push_back(fg);
  _interfaceGroups.push_back(fg);
  return fg->site;
}


const StorageSite&
Mesh::createBoundaryFaceGroup(const int size, const int offset, const int id, const string& boundaryType)
{
  shared_ptr<FaceGroup> fg(new FaceGroup(size,offset,_faces,id,boundaryType));
  _faceGroups.push_back(fg);
  _boundaryGroups.push_back(fg);
  return fg->site;
}


shared_ptr<Array<int> >
Mesh::createAndGetBNglobalToLocal() const
{
  if(!_boundaryNodeGlobalToLocalPtr)
  {
      const int nNodes = _nodes.getCount();
      _boundaryNodeGlobalToLocalPtr = shared_ptr<Array<int> >(new Array<int>(nNodes));
      Array<int>& globalToLocal = *_boundaryNodeGlobalToLocalPtr;
      globalToLocal = -1;
      int BoundaryNodeCount=0;
      int nLocal=0;
      foreach(const FaceGroupPtr fgPtr, getAllFaceGroups())
      {
	  const FaceGroup& fg = *fgPtr;
	  if (fg.groupType != "interior")
	  {
	      const StorageSite& BoundaryFaces = fg.site;
	      const CRConnectivity& BoundaryFaceNodes = getFaceNodes(BoundaryFaces);
	      const Array<int>& BFArray = BoundaryFaceNodes.getRow();
	      const Array<int>& BNArray = BoundaryFaceNodes.getCol();
	      const int nBFaces = BoundaryFaceNodes.getRowDim();
	      for(int i=0;i<nBFaces;i++)
	      {
		  for(int ip=BFArray[i];ip<BFArray[i+1];ip++)
		  {
		      const int j = BNArray[ip];
		      if (globalToLocal[j] == -1)
			globalToLocal[j] = nLocal++;
		  }
	      }
	  }
      }
      BoundaryNodeCount = nLocal;
  }
  return _boundaryNodeGlobalToLocalPtr;   
}


const StorageSite& Mesh::getBoundaryNodes() const 
{
  if(!_boundaryNodes)
  {
      const int nNodes = _nodes.getCount();
      shared_ptr<Array<int> > GlobalToLocalPtr = createAndGetBNglobalToLocal();
      Array<int>& GlobalToLocal = *GlobalToLocalPtr;
      int BoundaryNodeCount = 0;
      int nLocal = 0;
      for(int i=0;i<nNodes;i++)
      {
          if(GlobalToLocal[i] != -1)
	    nLocal++;
      }
      BoundaryNodeCount = nLocal;
      _boundaryNodes = new StorageSite(BoundaryNodeCount,0,0,0);
  }
  return *_boundaryNodes;
}


const ArrayBase& Mesh::getBNglobalToLocal() const
{
  return *(createAndGetBNglobalToLocal());
}


void Mesh::setConnectivity(const StorageSite& rowSite, const StorageSite& colSite,
	                       shared_ptr<CRConnectivity> conn)
{
  SSPair key(&rowSite,&colSite);
  _connectivityMap[key] = conn;
}

void Mesh::eraseConnectivity(const StorageSite& rowSite,
                             const StorageSite& colSite) const
{
  SSPair key(&rowSite,&colSite);
  _connectivityMap.erase(key);
}


const CRConnectivity&
Mesh::getAllFaceNodes() const
{
  SSPair key(&_faces,&_nodes);
  ConnectivityMap::const_iterator pos = _connectivityMap.find(key);
  if (pos != _connectivityMap.end())
    return *pos->second;
  throw CException("face nodes not defined");
}

const CRConnectivity&
Mesh::getAllFaceCells() const
{
  SSPair key(&_faces,&_cells);
  ConnectivityMap::const_iterator pos = _connectivityMap.find(key);
  if (pos != _connectivityMap.end())
    return *pos->second;
  throw CException("face cells not defined");
}

const CRConnectivity&
Mesh::getFaceCells(const StorageSite& faces) const
{
  SSPair key(&faces,&_cells);
  ConnectivityMap::const_iterator pos = _connectivityMap.find(key);
  if (pos != _connectivityMap.end())
    return *pos->second;

  shared_ptr<CRConnectivity> thisFaceCells =
    getAllFaceCells().createOffset(faces,faces.getOffset(),faces.getCount());
  _connectivityMap[key] = thisFaceCells;
  return *thisFaceCells;
}

const CRConnectivity&
Mesh::getFaceNodes(const StorageSite& faces) const
{
  SSPair key(&faces,&_nodes);
  ConnectivityMap::const_iterator pos = _connectivityMap.find(key);
  if (pos != _connectivityMap.end())
    return *pos->second;

  shared_ptr<CRConnectivity> thisFaceNodes =
    getAllFaceNodes().createOffset(faces,faces.getOffset(),faces.getCount());
  _connectivityMap[key] = thisFaceNodes;
  return *thisFaceNodes;
}

const CRConnectivity&
Mesh::getConnectivity(const StorageSite& from, const StorageSite& to) const
{
  SSPair key(&from,&to);
  ConnectivityMap::const_iterator pos = _connectivityMap.find(key);
  if (pos != _connectivityMap.end())
    return *pos->second;
  throw CException("connectivity not defined");
}

const CRConnectivity&
Mesh::getCellNodes() const
{
  SSPair key(&_cells,&_nodes);
 


  ConnectivityMap::const_iterator pos = _connectivityMap.find(key);

  if (pos != _connectivityMap.end())
    return *pos->second;


  const CRConnectivity& faceCells = getAllFaceCells();
  const CRConnectivity& faceNodes = getAllFaceNodes();
  
  SSPair keycf(&_cells,&_faces);
  shared_ptr<CRConnectivity> cellFaces = faceCells.getTranspose();
  shared_ptr<CRConnectivity> cellNodes = cellFaces->multiply(faceNodes,false);

  _connectivityMap[keycf] = cellFaces;
  _connectivityMap[key] = cellNodes;
  
  orderCellFacesAndNodes(*cellFaces, *cellNodes, faceNodes,
                         faceCells, *_coordinates);
  return *cellNodes;
}

const CRConnectivity&
Mesh::getCellFaces() const
{
  SSPair key(&_cells,&_faces);
  ConnectivityMap::const_iterator pos = _connectivityMap.find(key);
  if (pos != _connectivityMap.end())
    return *pos->second;

  const CRConnectivity& faceCells = getAllFaceCells();
  shared_ptr<CRConnectivity> cellFaces = faceCells.getTranspose();

  _connectivityMap[key] = cellFaces;
  return *cellFaces;
}

CRConnectivity&
Mesh::getAllFaceCells()
{
  SSPair key(&_faces,&_cells);
  ConnectivityMap::const_iterator pos = _connectivityMap.find(key);
  if (pos != _connectivityMap.end())
    return *pos->second;
  throw CException("face cells not defined");
}


const CRConnectivity&
Mesh::getCellCells() const
{
  SSPair key(&_cells,&_cells);
  ConnectivityMap::const_iterator pos = _connectivityMap.find(key);
  if (pos != _connectivityMap.end())
    return *pos->second;

  const CRConnectivity& faceCells = getAllFaceCells();
  const CRConnectivity& cellFaces = getCellFaces();
  shared_ptr<CRConnectivity> cellCells = cellFaces.multiply(faceCells,true);
  _connectivityMap[key] = cellCells;
  return *cellCells;
}

const CRConnectivity&
Mesh::getCellCells2() const
{
  if (!_cellCells2)
  { 
#ifdef FVM_PARALLEL
      //CRConnectivity constructor
      _cellCells2 = shared_ptr<CRConnectivity> ( new CRConnectivity(this->getCells(), this->getCells()) );
      //initCount
      _cellCells2->initCount();

      const int ncells = this->getCells().getSelfCount();
      multimap<int,int>::const_iterator it0;
      multimap<int,int>::const_iterator it1;
      //loop over inner cells
      for ( int n = 0; n < ncells; n++ ){ 
         set<int> setCells;
         for ( it0 = _cellCellsGlobal.equal_range(n).first; it0 != _cellCellsGlobal.equal_range(n).second; ++it0 ){
             const int cellID1 = it0->second;
             const int localID = _globalToLocal.find(cellID1)->second;
             setCells.insert( cellID1 );
             for (  it1 = _cellCellsGlobal.equal_range(localID).first; it1 != _cellCellsGlobal.equal_range(localID).second; ++it1 ){
                 const int cellID2 = it1->second;
                 setCells.insert(cellID2);
             }
         }
         //erase itself
         setCells.erase( (*_localToGlobal)[n]);
         _cellCells2->addCount(n,setCells.size());
      }
      //finish count
      _cellCells2->finishCount();
      //add cellcells2
      //loop over inner cells
      for ( int n = 0; n < ncells; n++ ){ 
         set<int> setCells;
         for ( it0 = _cellCellsGlobal.equal_range(n).first; it0 != _cellCellsGlobal.equal_range(n).second; it0++ ){
             const int cellID1 = it0->second;
             const int localID = _globalToLocal.find(cellID1)->second;
             setCells.insert( cellID1 );
             for (  it1 = _cellCellsGlobal.equal_range(localID).first; it1 != _cellCellsGlobal.equal_range(localID).second; it1++ ){
                 const int cellID2 = it1->second;
                 setCells.insert(cellID2);
             }
         }
         //erase itself
         setCells.erase((*_localToGlobal)[n]);
         foreach ( const set<int>::value_type globalID, setCells ){
             const int localCellID = _globalToLocal.find(globalID)->second;
            _cellCells2->add(n, localCellID);
         }
      }
      //finish add
      _cellCells2->finishAdd();
#endif

#ifndef FVM_PARALLEL
       const CRConnectivity& cellCells = getCellCells();
       _cellCells2 = cellCells.multiply(cellCells, true);
#endif



  }
  return *_cellCells2;
}

const CRConnectivity&
Mesh::getFaceCells2() const
{
  if (!_faceCells2)
  {
      const CRConnectivity& cellCells = getCellCells();
      const CRConnectivity& faceCells = getAllFaceCells();
      _faceCells2 = faceCells.multiply(cellCells, false);
  }
  return *_faceCells2;
}


void
Mesh::setFaceNodes(shared_ptr<CRConnectivity> faceNodes)
{
  SSPair key(&_faces,&_nodes);
  _connectivityMap[key] = faceNodes;
}


void
Mesh::setFaceCells(shared_ptr<CRConnectivity> faceCells)
{
  SSPair key(&_faces,&_cells);
  _connectivityMap[key] = faceCells;
}


const Array<int>&
Mesh::getIBFaceList() const
{
  if (_ibFaceList)  return (*_ibFaceList);
  throw CException("ib face list not defined");
}


void 
Mesh::createGhostCellSiteScatter(  const PartIDMeshIDPair& id, shared_ptr<StorageSite> site )
{
  _ghostCellSiteScatterMap.insert( pair<PartIDMeshIDPair, shared_ptr<StorageSite> >( id, site ) );

}

void 
Mesh::createGhostCellSiteGather( const PartIDMeshIDPair& id, shared_ptr<StorageSite> site )
{
  _ghostCellSiteGatherMap.insert( pair<PartIDMeshIDPair, shared_ptr<StorageSite> >( id, site ) );

}

void 
Mesh::createGhostCellSiteScatterLevel1(  const PartIDMeshIDPair& id, shared_ptr<StorageSite> site )
{
  _ghostCellSiteScatterMapLevel1.insert( pair<PartIDMeshIDPair, shared_ptr<StorageSite> >( id, site ) );

}

void 
Mesh::createGhostCellSiteGatherLevel1( const PartIDMeshIDPair& id, shared_ptr<StorageSite> site )
{
  _ghostCellSiteGatherMapLevel1.insert( pair<PartIDMeshIDPair, shared_ptr<StorageSite> >( id, site ) );

}


//this 
void 
Mesh::createCellColor()
{
    //cellColor color ghost cells respect to self-inner cells
   _cellColor      = shared_ptr< Array<int> > ( new Array<int>( _cells.getCount() ) );
    //cellColorOther color ghost cells in respect to other partition,
    //if partition interface has aligned with mesh interface this will be different than _cellColor
   _cellColorOther = shared_ptr< Array<int> > ( new Array<int>( _cells.getCount() ) );

   *_cellColor = -1;
   *_cellColorOther = -1;
   _isAssembleMesh = true;
}

void
Mesh::createLocalGlobalArray()
{
   _localToGlobal  = shared_ptr< Array<int> > ( new Array<int>( _cells.getCount() ) );
   *_localToGlobal  = -1;
}


void
Mesh::findCommonNodes(Mesh& other)
{
  StorageSite& nodes = _nodes;
  StorageSite& otherNodes = other._nodes;
  const int nNodes = nodes.getCount();
  const int nOtherNodes = otherNodes.getCount();

  const Array<VecD3>& coords = getNodeCoordinates();
  const Array<VecD3>& otherCoords = other.getNodeCoordinates();


  KSearchTree bNodeTree;

  // add all boundary nodes of this mesh to the tree
  {
      Array<bool> nodeMark(nNodes);
      nodeMark = false;
      foreach(const FaceGroupPtr fgPtr, getAllFaceGroups())
      {
          const FaceGroup& fg = *fgPtr;
          const StorageSite& faces = fg.site;
          if (fg.groupType!="interior")
          {
              const int nFaces = faces.getCount();
              const CRConnectivity& faceNodes = getFaceNodes(faces);
              for(int f=0; f<nFaces; f++)
              {
                  const int nFaceNodes = faceNodes.getCount(f);
                  for(int nn=0; nn<nFaceNodes; nn++)
                  {
                      const int n=faceNodes(f,nn);
                      if (!nodeMark[n])
                      {
                          nodeMark[n] = true;
                          bNodeTree.insert(coords[n],n);
                      }
                  }
              }
          }
      }
  }
  

  // loop over all the boundary nodes of the other mesh to find possible common ones
  Array<bool> nodeMark(nOtherNodes);

  typedef map<int,int> CommonNodesMap;
  CommonNodesMap commonNodesMap;

  Array<int> closest(2);
  
  nodeMark = false;

  
  foreach(const FaceGroupPtr fgPtr, other.getAllFaceGroups())
  {
      const FaceGroup& fg = *fgPtr;
      const StorageSite& faces = fg.site;
      if (fg.groupType!="interior")
      {
          const int nFaces = faces.getCount();
          const CRConnectivity& faceNodes = other.getFaceNodes(faces);

          for(int f=0; f<nFaces; f++)
          {
              const int nFaceNodes = faceNodes.getCount(f);
              for(int nn=0; nn<nFaceNodes; nn++)
              {
                  const int n=faceNodes(f,nn);
                  if (!nodeMark[n])
                  {
                      bNodeTree.findNeighbors(otherCoords[n],2,closest);
                      
                      double dist0 = mag(otherCoords[n] - coords[closest[0]]);
                      
                      // distance between the two closest point used as scale
                      
                      double distScale = mag(coords[closest[0]] - coords[closest[1]]);
                      
                      if (dist0 < distScale*epsilon)
                      {
                          // if another node has already been found as the
                          // closest for this one we have something wrong
                          if (commonNodesMap.find(closest[0]) != commonNodesMap.end())
                          {
                              throw CException("duplicate nodes on the mesh ?");
                          }
                          commonNodesMap.insert(make_pair(closest[0], n));
                      }
                      nodeMark[n] = true;
                  }
              }
          }
      }
  }

  const int nCommon = commonNodesMap.size();

  if (nCommon == 0)
    return;
  cout << "found " << nCommon << " common nodes " << endl;
  
  shared_ptr<IntArray> myCommonNodes(new IntArray(nCommon));
  shared_ptr<IntArray> otherCommonNodes(new IntArray(nCommon));

  int nc=0;
  foreach(CommonNodesMap::value_type& pos, commonNodesMap)
  {
      (*myCommonNodes)[nc] = pos.first;
      (*otherCommonNodes)[nc] = pos.second;
      nc++;
  }
  
  nodes.getCommonMap()[&otherNodes] = myCommonNodes;
  otherNodes.getCommonMap()[&nodes] = otherCommonNodes;
  
}

void
Mesh::findCommonFaces(StorageSite& faces, StorageSite& otherFaces,
                      const GeomFields& geomFields)
{
  const int count(faces.getCount());
  if (count != otherFaces.getCount())
    throw CException("face groups are not of the same length");

  const Array<VecD3>& coords =
    dynamic_cast<const Array<VecD3>& >(geomFields.coordinate[faces]);
  
  const Array<VecD3>& otherCoords =
    dynamic_cast<const Array<VecD3>& >(geomFields.coordinate[otherFaces]);

  const Array<VecD3>& area =
    dynamic_cast<const Array<VecD3>& >(geomFields.area[faces]);
  
  const Array<VecD3>& otherArea =
    dynamic_cast<const Array<VecD3>& >(geomFields.area[otherFaces]);

  KSearchTree thisFacesTree(coords);
  
  Array<int> closest(2);

  shared_ptr<IntArray> myCommonFaces(new IntArray(count));
  shared_ptr<IntArray> otherCommonFaces(new IntArray(count));
  
  for(int f=0; f<count; f++)
  {
      thisFacesTree.findNeighbors(otherCoords[f],2,closest);
      
      const int closestFace = closest[0];
      double dist0 = mag(otherCoords[f] - coords[closestFace]);
      
      // distance between the two closest point used as scale
      
      double distScale = mag(coords[closest[0]] - coords[closest[1]]);
      
      if (dist0 < distScale*epsilon)
      {
          double crossProductMag(mag2(cross(otherArea[f],area[closestFace])));
          if (crossProductMag > mag2(otherArea[f])*epsilon)
            throw CException("cross product is not small");
          
          (*otherCommonFaces)[f] = closestFace;
          (*myCommonFaces)[closestFace] = f;
          
      }
  }
  
  faces.getCommonMap()[&otherFaces] = myCommonFaces;
  otherFaces.getCommonMap()[&faces] = otherCommonFaces;
  
}


Mesh*
Mesh::extractBoundaryMesh()
{
  StorageSite& nodes = _nodes;
  const Array<VecD3>& coords = getNodeCoordinates();

  const int nodeCount = nodes.getCount();
  Array<int> globalToLocalNodes(nodeCount);

  globalToLocalNodes = -1;
  int bMeshNodeCount=0;
  int bMeshFaceCount=0;
  foreach(const FaceGroupPtr fgPtr, getAllFaceGroups())
  {
      const FaceGroup& fg = *fgPtr;
      const StorageSite& faces = fg.site;
      if (fg.groupType!="interior")
      {
          const int nFaces = faces.getCount();
          const CRConnectivity& faceNodes = getFaceNodes(faces);
          for(int f=0; f<nFaces; f++)
          {
              const int nFaceNodes = faceNodes.getCount(f);
              for(int nn=0; nn<nFaceNodes; nn++)
              {
                  const int n=faceNodes(f,nn);
                  if (globalToLocalNodes[n] == -1)
                  {
                      globalToLocalNodes[n] = bMeshNodeCount++;
                  }
              }
          }
          bMeshFaceCount += nFaces;
      }
  }


  Mesh *bMesh = new Mesh(_dimension);


  StorageSite& bMeshFaces = bMesh->getFaces();
  StorageSite& bMeshNodes = bMesh->getNodes();
  bMeshFaces.setCount( bMeshFaceCount );
  bMeshNodes.setCount( bMeshNodeCount );
  
  bMesh->createBoundaryFaceGroup(bMeshFaceCount,0,0,"wall");
  
  //setting coordinates
  shared_ptr< Array<VecD3> > bMeshCoordPtr( new Array< VecD3 > ( bMeshNodeCount ) );

  shared_ptr<IntArray> myCommonNodes(new IntArray(bMeshNodeCount));
  shared_ptr<IntArray> otherCommonNodes(new IntArray(bMeshNodeCount));

  for(int n=0; n<nodeCount; n++)
  {
      const int nLocal = globalToLocalNodes[n];
      if (nLocal >=0)
      {
          (*bMeshCoordPtr)[nLocal] = coords[n];
          (*myCommonNodes)[nLocal] = nLocal;
          (*otherCommonNodes)[nLocal] = n;
      }
  }
  nodes.getCommonMap()[&bMeshNodes] = myCommonNodes;
  bMeshNodes.getCommonMap()[&nodes] = otherCommonNodes;
         
  bMesh->setCoordinates( bMeshCoordPtr );
  
  //faceNodes constructor
  shared_ptr<CRConnectivity> bFaceNodes( new CRConnectivity(bMeshFaces,
                                                            bMeshNodes) );
  
  bFaceNodes->initCount();

  bMeshFaceCount=0;
  
  foreach(FaceGroupPtr fgPtr, getAllFaceGroups())
  {
      FaceGroup& fg = *fgPtr;
      StorageSite& faces = const_cast<StorageSite&>(fg.site);
      if (fg.groupType!="interior")
      {
          const int nFaces = faces.getCount();
          const CRConnectivity& faceNodes = getFaceNodes(faces);

          shared_ptr<IntArray> myCommonFaces(new IntArray(nFaces));
          shared_ptr<IntArray> otherCommonFaces(new IntArray(nFaces));

          for(int f=0; f<nFaces; f++)
          {
              const int nFaceNodes = faceNodes.getCount(f);
              bFaceNodes->addCount(bMeshFaceCount,nFaceNodes);
              (*myCommonFaces)[f] = bMeshFaceCount;
              (*otherCommonFaces)[f] = f;
              bMeshFaceCount++;
          }

          faces.getCommonMap()[&bMeshFaces] = myCommonFaces;
          bMeshFaces.getCommonMap()[&faces] = otherCommonFaces;
      }
  }

  bFaceNodes->finishCount();
  bMeshFaceCount=0;

  foreach(const FaceGroupPtr fgPtr, getAllFaceGroups())
  {
      const FaceGroup& fg = *fgPtr;
      const StorageSite& faces = fg.site;
      if (fg.groupType!="interior")
      {
          const int nFaces = faces.getCount();
          const CRConnectivity& faceNodes = getFaceNodes(faces);
          for(int f=0; f<nFaces; f++)
          {
              const int nFaceNodes = faceNodes.getCount(f);
              for(int nn=0; nn<nFaceNodes; nn++)
              {
                  const int n=faceNodes(f,nn);
                  const int nLocal = globalToLocalNodes[n];
                  bFaceNodes->add(bMeshFaceCount,nLocal);
              }
              bMeshFaceCount++;
          }
      }
  }
  
  bFaceNodes->finishAdd();
  //setting faceNodes
  SSPair key(&bMeshFaces,&bMeshNodes);
  bMesh->_connectivityMap[key] = bFaceNodes;

  return bMesh;

}

Mesh*
Mesh::extrude(int nz, double zmax)
{
  if (_dimension != 2)
    throw CException("can only extrude two dimensional mesh");

  const int myNCells = _cells.getSelfCount();
  const int myNFaces = _faces.getSelfCount();
  const int myNInteriorFaces = _interiorFaceGroup->site.getCount();
  const int myNBoundaryFaces = myNFaces - myNInteriorFaces;
  const int myNNodes = _nodes.getSelfCount();
  const CRConnectivity& myCellNodes = getCellNodes();
 
  const int eNInteriorFaces_rib = nz*myNInteriorFaces;
  const int eNInteriorFaces_cap = (nz-1)*myNCells;
  
  const int eNInteriorFaces = eNInteriorFaces_rib + eNInteriorFaces_cap;
  const int eNBoundaryFaces = nz*myNBoundaryFaces + 2*myNCells;
  const int eNFaces = eNInteriorFaces + eNBoundaryFaces;
  const int eNInteriorCells = nz*myNCells;
  const int eNBoundaryCells = eNBoundaryFaces;
  
  const int eNCells = eNInteriorCells + eNBoundaryCells;
  const int eNNodes = (nz+1)*myNNodes;
  
  Mesh *eMesh = new Mesh(3);

  // set nodes
  StorageSite& eMeshNodes = eMesh->getNodes();
  eMeshNodes.setCount( eNNodes );

  const Array<VecD3>& myCoords = getNodeCoordinates();
  shared_ptr< Array<VecD3> > eMeshCoordPtr( new Array< VecD3 > ( eNNodes ) );
  Array<VecD3>& eCoords = *eMeshCoordPtr;
  const double dz = zmax/nz;

  const double z0 = -zmax/2.0;
  for(int k=0, en=0; k<=nz; k++)
  {
      const double z = z0 + dz*k;
      for(int n=0; n<myNNodes; n++)
      {
          eCoords[en][0] = myCoords[n][0];
          eCoords[en][1] = myCoords[n][1];
          eCoords[en][2] = z;
          en++;
      }
  }

  eMesh->setCoordinates( eMeshCoordPtr );
  
          
  // set cells
  
  StorageSite& eMeshCells = eMesh->getCells();
  eMeshCells.setCount( eNInteriorCells, eNBoundaryCells);
  

  // set faces
  
  StorageSite& eMeshFaces = eMesh->getFaces();
  eMeshFaces.setCount(eNFaces);


  shared_ptr<CRConnectivity> eFaceNodes( new CRConnectivity(eMeshFaces,
                                                            eMeshNodes) );

  shared_ptr<CRConnectivity> eFaceCells( new CRConnectivity(eMeshFaces,
                                                            eMeshCells) );

  // set counts for face cells and face nodes

  
  eFaceNodes->initCount();
  eFaceCells->initCount();

  int f = 0;

  // rib faces first
  for(; f<eNInteriorFaces_rib; f++)
  {
      eFaceNodes->addCount(f,4);
      eFaceCells->addCount(f,2);
  }

  // interior cap faces
  for(int k=1; k<nz; k++)
  {
      for(int c=0; c<myNCells; c++)
      {
          eFaceNodes->addCount(f, myCellNodes.getCount(c));
          eFaceCells->addCount(f, 2);
          f++;
      }
  }

  eMesh->createInteriorFaceGroup(eNInteriorFaces);

  foreach(const FaceGroupPtr fgPtr, getAllFaceGroups())
  {
      const FaceGroup& fg = *fgPtr;
      const StorageSite& faces = fg.site;
      if (fg.groupType!="interior")
      {
          const int nBFaces = faces.getCount();

          eMesh->createBoundaryFaceGroup(nBFaces*nz,  f, fg.id, fg.groupType);

          for(int k=0;  k<nz; k++)
            for(int bf=0; bf<nBFaces; bf++)
            {
                eFaceNodes->addCount(f,4);
                eFaceCells->addCount(f,2);
                f++;
            }
      }
  }
  

  // z = 0 faces
  eMesh->createBoundaryFaceGroup(myNCells,  f, 10000, "wall");
  for(int c=0; c<myNCells; c++)
  {
      eFaceNodes->addCount(f, myCellNodes.getCount(c));
      eFaceCells->addCount(f, 2);
      f++;
  }
  
  // z = zmax faces
  eMesh->createBoundaryFaceGroup(myNCells,  f, 10001, "wall");
  for(int c=0; c<myNCells; c++)
  {
      eFaceNodes->addCount(f, myCellNodes.getCount(c));
      eFaceCells->addCount(f, 2);
      f++;
  }
  
  eFaceNodes->finishCount();
  eFaceCells->finishCount();


  // now set the indices of face cells and nodes

  
  f = 0;

  
  // rib faces first

  const CRConnectivity& myFaceNodes = getAllFaceNodes();
  const CRConnectivity& myFaceCells = getAllFaceCells();
  
  for(int k=0; k<nz; k++)
  {
      const int eCellOffset = k*myNCells;
      const int eNodeOffset = k*myNNodes;
      
      for(int myf=0; myf<myNInteriorFaces; myf++)
      {
          const int myNode0 = myFaceNodes(myf,0);
          const int myNode1 = myFaceNodes(myf,1);
          
          const int myCell0 = myFaceCells(myf,0);
          const int myCell1 = myFaceCells(myf,1);
          
          
          eFaceNodes->add(f, myNode0 + eNodeOffset);
          eFaceNodes->add(f, myNode1 + eNodeOffset );
          eFaceNodes->add(f, myNode1 + eNodeOffset + myNNodes);
          eFaceNodes->add(f, myNode0 + eNodeOffset + myNNodes);
          
          eFaceCells->add(f, myCell0 + eCellOffset);
          eFaceCells->add(f, myCell1 + eCellOffset);
          
          f++;
      }
  }
  
  // interior cap faces
  for(int k=1; k<nz; k++)
  {
      const int eCellOffset = k*myNCells;
      const int eNodeOffset = k*myNNodes;

      for(int c=0; c<myNCells; c++)
      {
          const int nCellNodes = myCellNodes.getCount(c);
          for(int nnc=0; nnc<nCellNodes; nnc++)
          {
              eFaceNodes->add(f, myCellNodes(c,nnc) + eNodeOffset);
          }

          eFaceCells->add(f, c + eCellOffset - myNCells);
          eFaceCells->add(f, c + eCellOffset);
          f++;
      }
  }


  // counter for boundary faces 
  int ebf = 0;


  // rib boundary faces
  foreach(const FaceGroupPtr fgPtr, getAllFaceGroups())
  {
      const FaceGroup& fg = *fgPtr;
      const StorageSite& faces = fg.site;
      if (fg.groupType!="interior")
      {
          const int nBFaces = faces.getCount();

          const CRConnectivity& bFaceNodes = getFaceNodes(faces);
          const CRConnectivity& bFaceCells = getFaceCells(faces);
          
          for(int k=0;  k<nz; k++)
          {
              const int eCellOffset = k*myNCells;
              const int eNodeOffset = k*myNNodes;
              for(int bf=0; bf<nBFaces; bf++)
              {
                  const int myNode0 = bFaceNodes(bf,0);
                  const int myNode1 = bFaceNodes(bf,1);
                  
                  const int myCell0 = bFaceCells(bf,0);
          
                  
                  eFaceNodes->add(f, myNode0 + eNodeOffset);
                  eFaceNodes->add(f, myNode1 + eNodeOffset);
                  eFaceNodes->add(f, myNode1 + eNodeOffset + myNNodes);
                  eFaceNodes->add(f, myNode0 + eNodeOffset + myNNodes);
                  
                  eFaceCells->add(f, myCell0 + eCellOffset);
                  eFaceCells->add(f, ebf + eNInteriorCells);
                  
                  f++;
                  ebf++;
              }
          }
      }
      
  }


  // z = 0 faces
  {  
      const int eCellOffset = 0;
      const int eNodeOffset = 0;

      for(int c=0; c<myNCells; c++)
      {
          const int nCellNodes = myCellNodes.getCount(c);
          
          // reverse order of face nodes
          for(int nnc=nCellNodes-1; nnc>=0; nnc--)
          {
              eFaceNodes->add(f, myCellNodes(c,nnc) + eNodeOffset);
          }
          
          eFaceCells->add(f, c + eCellOffset);
          eFaceCells->add(f, ebf + eNInteriorCells);
          f++;
          ebf++;
      }
      
  }
  
  // z = zmax faces
  {  
      const int eCellOffset = (nz-1)*myNCells;
      const int eNodeOffset = nz*myNNodes;

      for(int c=0; c<myNCells; c++)
      {
          const int nCellNodes = myCellNodes.getCount(c);
          
          // reverse order of face nodes
          for(int nnc=0; nnc<nCellNodes; nnc++)
          {
              eFaceNodes->add(f, myCellNodes(c,nnc) + eNodeOffset);
          }
          
          eFaceCells->add(f, c + eCellOffset);
          eFaceCells->add(f, ebf + eNInteriorCells);
          f++;
          ebf++;
      }
  }
  
  eFaceNodes->finishAdd();
  eFaceCells->finishAdd();

 //setting faceNodes
  SSPair key1(&eMeshFaces,&eMeshNodes);
  eMesh->_connectivityMap[key1] = eFaceNodes;

  SSPair key2(&eMeshFaces,&eMeshCells);
  eMesh->_connectivityMap[key2] = eFaceCells;

  return eMesh;
}

  

