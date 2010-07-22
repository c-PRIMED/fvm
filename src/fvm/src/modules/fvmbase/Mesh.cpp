#include "Mesh.h"
#include "StorageSite.h"
#include "CRConnectivity.h"
#include "Cell.h"
#include <cassert>
#include "KSearchTree.h"
#include "GeomFields.h"

#define epsilon 1e-6

Mesh::Mesh(const int dimension, const int id):
  _dimension(dimension),
  _id(id),
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

Mesh::Mesh( const int dimension, const int id,
            const Array<VecD3>&  faceNodesCoord ): 
  _dimension(dimension),
  _id(id),
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
      const CRConnectivity& cellCells = getCellCells();
      _cellCells2 = cellCells.multiply(cellCells, true);
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


  Mesh *bMesh = new Mesh(_dimension,0);


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
