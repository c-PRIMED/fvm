#include "Mesh.h"
#include "StorageSite.h"
#include "CRConnectivity.h"
#include "Cell.h"
#include <cassert>

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

  Mesh::Mesh( const int dimension, const int id, const Array<VecD3>&  faceNodesCoord ): 
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
   int faceNodeCount = 4;
   int totNodes      = faceNodesCoord.getLength(); // counting duplicate nodes as well
   int totFaces      = totNodes / faceNodeCount;
   //check if this is corect integer division
   assert( (faceNodeCount*totFaces) == totNodes );
   //set sites
   StorageSite& faceSite = getFaces();
   StorageSite& nodeSite = getNodes();
   faceSite.setCount( totFaces );
   nodeSite.setCount( totNodes );
   //interior face group (we have only one interface for this
   createInteriorFaceGroup(1);

   //setting coordinates
   shared_ptr< Array<VecD3> > coord( new Array< VecD3 > ( totNodes ) );
   *coord = faceNodesCoord;
   setCoordinates( coord );

   //faceNodes constructor
   shared_ptr<CRConnectivity> faceNodes( new CRConnectivity(faceSite, nodeSite) );
   //addCount
   for ( int i = 0; i < totFaces; i++ )
        faceNodes->addCount(i, faceNodeCount); 
   //finish count
   faceNodes->finishCount();
   //add operation 
   int face     = 0;
   int nodeIndx = 0;
   for( int i = 0; i < totFaces; i++ ){
      for( int j =0; j < faceNodeCount; j++ ){
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
  
  orderCellFacesAndNodes(*cellFaces, *cellNodes, faceNodes, faceCells);
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

