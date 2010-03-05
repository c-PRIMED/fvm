
#include "Mesh.h"
#include "StorageSite.h"
#include "CRConnectivity.h"
#include "Cell.h"

Mesh::Mesh(const int dimension, const int id):
  _dimension(dimension),
  _id(id),
  _cells(0),
  _faces(0),
  _nodes(0),
  _boundaryNodes(0),
  _ibFaces(0),
  _boundaryNodeGlobalToLocalPtr(),
  _interiorFaceGroup(),
  _faceGroups(),
  _boundaryGroups(),
  _interfaceGroups(),
  _connectivityMap(),
  _coordinates(),
  _ibType(),
  _ibFaceList(),
  _numOfAssembleMesh(1),
  _isAssembleMesh(false)
{
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
      foreach(const FaceGroupPtr fgPtr, getBoundaryFaceGroups())
      {
	  const FaceGroup& fg = *fgPtr;
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
Mesh::getIBType() const
{
  return getOrCreateIBType();
}

void
Mesh::setIBTypeForCell(const int c, const int type)
{
  Array<int>& ibType = getOrCreateIBType();
  ibType[c]=type;
}

int
Mesh::getIBTypeForCell(const int c) const
{
  Array<int>& ibType = getOrCreateIBType();
  return ibType[c];
}

Array<int>&
Mesh::getOrCreateIBType() const
{
  if (!_ibType)
  {
      _ibType = shared_ptr<Array<int> >(new Array<int>(_cells.getCount()));
      *_ibType = IBTYPE_FLUID;
  }
  return *_ibType;
}


#if 0
Mesh::VecD3
Mesh::getCellCoordinate(const int c) const
{
  return (*_coordinates)[c];
}
#endif

void
Mesh::createIBFaceList(const int size) const
{
  _ibFaceList = shared_ptr<Array<int> >(new Array<int> (size));
}

const Array<int>&
Mesh::getIBFaceList() const
{
  if (_ibFaceList)  return (*_ibFaceList);
  throw CException("ib face list not defined");
}

void
Mesh::addIBFace(const int i, const int c)
{
  (*_ibFaceList)[i]=c;
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

