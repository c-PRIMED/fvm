#include "misc.h"
#include "FluentReader.h"
#include "CRConnectivity.h"
#include "StorageSite.h"
#include "OneToOneIndexMap.h"
#include "Mesh.h"
#include "OneToOneIndexMap.h"

#include "Cell.h"

enum
  {
    READ_SIZES,
    READ_COUNTS,
    READ_MESH,
    READ_DATA
  };


FluentReader::FluentReader(const string& fileName) :
  SchemeReader(fileName),
  _dimension(0),
  _numCells(0),
  _numFaces(0),
  _numNodes(0),
  _numBoundaryFaces(0),
  _cells(0),
  _faces(0),
  _nodes(0),
  _faceNodes(),
  _faceCells(),
  _cellFaces(),
  _cellNodes(),
  _nodeCells(),
  _faceZones(),
  _cellZones(),
  _coords(0),
  _rpVarStringLength(0)
{}

FluentReader::~FluentReader()
{}

void FluentReader::readVectorData(Array<Vec3>& a,
                                  const int iBeg,
                                  const int iEnd,
                                  const bool isBinary, const bool isDP)
{
  moveToListOpen();
  
  const int count = iEnd-iBeg+1;
  const int nread = _dimension*count;
  if (isDP)
  {
      double *buff;
      if (_dimension == 3)
        buff = (double*) a.getData() + (iBeg-1)*_dimension;
      else
        buff = new double[nread];
      
      if (isBinary)
      {
          if (nread != (int)fread(buff,sizeof(double),nread,_fp))
            cerr << "error reading dp binary nodes" << endl;
      }
      else
      {
          double *bp = buff;
          for(int i=0; i<nread; i++)
            if (fscanf(_fp,"%le",bp++) != 1)
              cerr << "error reading dp formatted nodes" << endl;
      }
      
      if (_dimension == 2)
      {
          double *bp = buff;
          for(int i=iBeg; i<=iEnd; i++)
          {
              a[i-1][0] = *bp++;
              a[i-1][1] = *bp++;
          }
          
          delete [] buff;
      }
  }
  else
  {
      float *buff = new float[nread];
      if (isBinary)
      {
          if (nread != (int)fread(buff,sizeof(float),nread,_fp))
            cerr << "error reading sp binary nodes" << endl;
      }
      else
      {
          float *bp = buff;
          for(int i=0; i<nread; i++)
            if (fscanf(_fp,"%e",bp++) != 1)
              cerr << "error reading sp formatted nodes" << endl;
      }
      
      float *bp = buff;
      for(int i=iBeg; i<=iEnd; i++)
        for(int d=0; d<_dimension; d++)
          a[i-1][d] = *bp++;
      
      delete [] buff;
  }
}
                                  
void
FluentReader::readNodes(const int pass, const bool isBinary,
                        const bool isDP, const int sectionID)
{
  int threadId, iBeg, iEnd, type, dummy;
  readHeader(threadId,iBeg,iEnd,type,dummy);
  if (pass == READ_SIZES)
  {
      if (threadId == 0)
      {
          _numNodes=iEnd;
      }
      else 
      {
          moveToListOpen();
      }
      if (isBinary)
        closeSectionBinary(sectionID);
      else
        closeSection();
  }
  else if (pass == READ_COUNTS)
  {
      if (threadId != 0)
      {
          moveToListOpen();
      }
      if (isBinary)
        closeSectionBinary(sectionID);
      else
        closeSection();
  }
  else if (pass == READ_MESH)
  {
      if (threadId != 0)
      {
          readVectorData(_coords,iBeg,iEnd,isBinary,isDP);
      }
      if (isBinary)
        closeSectionBinary(sectionID);
      else
        closeSection();
  }
  
  return;
}

void
FluentReader::readCells(const int pass, const bool isBinary,
                        const int sectionID)
{
  int threadId, iBeg, iEnd, type, dummy;
  readHeader(threadId,iBeg,iEnd,type,dummy);
  if (pass == READ_SIZES)
  {
      if (threadId == 0)
      {
         _numCells=iEnd;
      }
      else if ((type == 1) || (type == 17))
      {
          FluentCellZone *cz = new FluentCellZone();
          cz->ID=threadId;
          cz->iBeg=iBeg-1;
          cz->iEnd=iEnd-1;
          cz->threadType=type;
          _cellZones[threadId]=cz;
      }
      else if (type == 32)
      {
          _numCells -= (iEnd-iBeg+1);
      }
      else
      {
          throw CException("cell thread type not handled"); 
      }
  }
  if (isBinary)
    closeSectionBinary(sectionID);
  else
    closeSection();
}


void
FluentReader::readFaces(const int pass, const bool isBinary,
                        const int sectionID)
{
  int threadId, iBeg, iEnd, type, shape;
  readHeader(threadId,iBeg,iEnd,type,shape);

  //  cerr <<  " thread id " << threadId << endl;
  if (pass == READ_SIZES)
  {
      if (threadId == 0)
      {
          _numFaces=iEnd;
      }
      else 
      {
          if ((type != 0) && (type != 31))
          {
              FluentFaceZone *fz = new FluentFaceZone();
              fz->ID=threadId;
              fz->iBeg=iBeg-1;
              fz->iEnd=iEnd-1;
              fz->threadType=type;
              _faceZones[threadId]=fz;
              
              
              moveToListOpen();
              for(int f=iBeg; f<=iEnd; f++)
              {
                  int numNodes = shape;
                  if (shape == 0 || shape == 5)
                    numNodes = readInt(isBinary);
                  
                  skipInt(numNodes,isBinary);
                  int c0 = readInt(isBinary);
                  int c1 = readInt(isBinary);

                  if (c0 == 0 || c1 == 0)
                    _numBoundaryFaces++;
              }
          }
          else
          {
              _numFaces -= (iEnd-iBeg+1);
          }
      }
      if (isBinary)
        closeSectionBinary(sectionID);
      else
        closeSection();
  }
  else if (pass == READ_COUNTS)
  {
      if  (threadId != 0)
        moveToListOpen();

      if ((threadId != 0) && (type != 0) && (type != 31))
      {
          CRConnectivity& faceNodes = *_faceNodes;
          CRConnectivity& faceCells = *_faceCells;
          if (shape < 0) shape = _dimension;
          
          for(int f=iBeg; f<=iEnd; f++)
          {
              int numNodes = shape;
              if (shape == 0 || shape == 5)
                numNodes = readInt(isBinary);
                  
                  faceNodes.addCount(f-1,numNodes);

                  skipInt(numNodes+2,isBinary);
                  //int c0 = readInt(isBinary);
                  //int c1 = readInt(isBinary);

                  //                  if (c0 == 0 || c1 == 0)
                  //  faceCells.addCount(f-1,1);
                  //else
                    faceCells.addCount(f-1,2);
          }
      }

      if (isBinary)
        closeSectionBinary(sectionID);
      else
        closeSection();
      
        //      if  (threadId != 0)
        //  moveToListClose();
  }
  else if (pass == READ_MESH)
  {
      if  (threadId != 0)
        moveToListOpen();

      if ((threadId != 0) && (type != 0) && (type != 31))
      {
          CRConnectivity& faceNodes = *_faceNodes;
          CRConnectivity& faceCells = *_faceCells;
          if (shape < 0) shape = _dimension;

          const int maxNodes = 100;
          int fnodes[maxNodes];
          for(int f=iBeg; f<=iEnd; f++)
          {
              int numNodes = shape;
              if (shape == 0 || shape == 5)
                numNodes = readInt(isBinary);

              // msvc doesn't like this
              //int fnodes[numNodes];
              bool reverseNodes = _dimension == 3;
              
              for(int i=0; i<numNodes; i++)
                fnodes[i] = readInt(isBinary)-1;

              int c0 = readInt(isBinary);
              int c1 = readInt(isBinary);

              if (c0 == 0) reverseNodes = !reverseNodes;
              
              if (c0 != 0)
                faceCells.add(f-1,c0-1);
              if (c1 != 0)
                faceCells.add(f-1,c1-1);

              if (c0 ==0 || c1==0)
              {
                  faceCells.add(f-1,_numCells+_numBoundaryFaces);
                  _numBoundaryFaces++;
              }
              
              if (reverseNodes)
                for(int i=0; i<numNodes; i++)
                  faceNodes.add(f-1,fnodes[numNodes-i-1]);
              else
                for(int i=0; i<numNodes; i++)
                  faceNodes.add(f-1,fnodes[i]);
                                
          }
      }

      //     if  (threadId != 0)
      //  moveToListClose();
        if (isBinary)
    closeSectionBinary(sectionID);
  else
    closeSection();

  }
  
#if 0
  if (isBinary)
    closeSectionBinary(sectionID);
  else
    closeSection();
#endif
  return;
}

void
FluentReader::readFacePairs(const int pass, const bool isBinary,
                            const int sectionID)
{
#if 0
  int iBeg, iEnd, leftID, rightID, dummy;
  readHeader(iBeg,iEnd,leftID,rightID,dummy);

  //  cerr <<  " thread id " << threadId << endl;
  if (pass == READ_SIZES)
  {
      const int count = iEnd-iBeg+1;
      Array<int>* leftFaces = new Array<int>(count);
      Array<int>* rightFaces = new Array<int>(count);
      
              
      moveToListOpen();
      for(int n=0; n<count; n++)
      {
          int leftF = readInt(isBinary);
          int rightF = readInt(isBinary);

          (*leftFaces)[n]=leftF;
          (*rightFaces)[n]=rightF;
      }
      if (isBinary)
        closeSectionBinary(sectionID);
      else
        closeSection();

      FluentFacePairs *fp = new FluentFacePairs();
      fp->setVar("count",count);
      fp->setVar("leftID",leftID);
      fp->setVar("leftFaces",*leftFaces);
      fp->setVar("rightID",rightID);
      fp->setVar("rightFaces",*rightFaces);
      

      _facePairs.add(make_pair(leftID,rightID),fp);
  }
  else if ((pass == READ_COUNTS) || (pass == READ_MESH))
  {
      moveToListOpen();
      if (isBinary)
        closeSectionBinary(sectionID);
      else
        closeSection();
  }
#endif
}

void FluentReader::read(const int pass)
{
  int id;
  while ((id = getNextSection()) != EOF)
  {
      bool isBinary =  (id > 1000);
      bool isDP = (id > 3000);

      //      cerr << "reading section " << id << endl;
      
      switch (id %1000)
      {
      case 0:
      case 1:
        moveToListClose();
        break;

      case 2:
        if (pass == READ_SIZES)
        {
            fscanf(_fp,"%d",&_dimension);
        }
        else
          moveToListClose();
        break;

      case 37:
        if (pass == READ_SIZES)
        {
            _rpVarStringLength = readListLength();
        }
        else if (pass == READ_COUNTS)
        {
            char* rpVarString = new char[_rpVarStringLength];
            readList(rpVarString);
            _rpVars=string(rpVarString,_rpVarStringLength);
            delete [] rpVarString;
        }
        else 
        {
            closeSection();
        }
        break;
        
      case 39:
      {
          int zoneId;
          char zoneName[256], zoneType[80];
          moveToListOpen();
          fscanf(_fp,"%d%s",&zoneId,zoneType);

          // need to read zoneName this way because it may or may not
          // be followed by an int
          {
              int n=0;
              char c;
              while ((c = getc(_fp)) != EOF)
              {
                  if (isspace(c))
                  {
                      if (n>0)
                      {
                          moveToListClose();
                          break;
                      }
                  }
                  else if (c == ')')
                    break;
                  else
                  {
                      if (c == '-') c='_';
                      zoneName[n++]=c;
                  }
              }
              zoneName[n]='\0';
          }
          
          if (pass == READ_SIZES)
          {
              if (_cellZones.find(zoneId) != _cellZones.end())
              {
                  FluentCellZone *z = _cellZones[zoneId];
                  z->zoneName=string(zoneName,strlen(zoneName));
                  z->zoneType=string(zoneType,strlen(zoneType));
              }
              else if (_faceZones.find(zoneId) != _faceZones.end())
              {
                  FluentFaceZone *z = _faceZones[zoneId];
                  z->zoneName=string(zoneName,strlen(zoneName));
                  z->zoneType=string(zoneType,strlen(zoneType));
              }
              
              _zoneVarStringLength[zoneId] = readListLength();
          }
          else if (pass == READ_COUNTS)
          {
              char* zoneVarString = new char[_zoneVarStringLength[zoneId]];
              readList(zoneVarString);
              if (_cellZones.find(zoneId) != _cellZones.end())
              {
                  FluentCellZone *z = _cellZones[zoneId];
                  z->zoneVars = 
                    string(zoneVarString,_zoneVarStringLength[zoneId]);
              }
              else if (_faceZones.find(zoneId) != _faceZones.end())
              {
                  FluentFaceZone *z = _faceZones[zoneId];
                  z->zoneVars = 
                    string(zoneVarString,_zoneVarStringLength[zoneId]);
              }
              delete [] zoneVarString;
          }
          else 
          {
              closeSection();
          }
          break;
      }
      
      case 10:
        readNodes(pass,isBinary,isDP,id);
        break;
        
      case 12:
        readCells(pass,isBinary,id);
        break;
        
      case 13:
        readFaces(pass,isBinary,id);
        break;

#if 0
      case 18:
        readFacePairs(pass,isBinary,id);
        break;
#endif
        
      default:
        if (isBinary)
          closeSectionBinary(id);
        else
          closeSection();
        break;
      }
  }
}


void
FluentReader::readMesh()
{
  // first pass to find the sizes of cell zones etc
  read(READ_SIZES);

  _cells.setCount(_numCells,_numBoundaryFaces);
  _faces.setCount(_numFaces);
  _nodes.setCount(_numNodes);

  _faceNodes = shared_ptr<CRConnectivity>(new CRConnectivity(_faces,_nodes));
  _faceCells = shared_ptr<CRConnectivity>(new CRConnectivity(_faces,_cells));

  _faceNodes->initCount();
  _faceCells->initCount();

  // rewind and read to determine counts for connectivities
  resetFilePtr();
  read(READ_COUNTS);

  _faceCells->finishCount();
  _faceNodes->finishCount();
  
  _coords.resize(_numNodes);
  _coords.zero();
  
  _numBoundaryFaces = 0;

  // finally rewind once again and read in all the info
  
  resetFilePtr();
  read(READ_MESH);

  _faceNodes->finishAdd();
  _faceCells->finishAdd();
}

int
FluentReader::getCellZoneID(const int c) const
{
  foreach(const CellZonesMap::value_type& pos, _cellZones)
  {
      const FluentCellZone& cz = *(pos.second);
      if (c >= cz.iBeg && c <= cz.iEnd)
        return cz.ID;
  }
  if (c < _numCells + _numBoundaryFaces)
    return 0;
  throw CException("getCellZoneID: invalid cell id");
}

const CRConnectivity&
FluentReader::getCellFaces()
{
  if (!_cellFaces)
  {
      _cellFaces = _faceCells->getTranspose();
  }
  return *_cellFaces;
}

const CRConnectivity&
FluentReader::getCellNodes()
{
  if (!_cellNodes)
  {
      const CRConnectivity& cellFaces = getCellFaces();
      _cellNodes = cellFaces.multiply(*_faceNodes,false);
  }
  return *_cellNodes;
}

const CRConnectivity&
FluentReader::getNodeCells()
{
  if (!_nodeCells)
  {
      const CRConnectivity& cellNodes = getCellNodes();
      _nodeCells = cellNodes.getTranspose();
  }
  return *_nodeCells;
}

void
FluentReader::buildZones()
{
  const CRConnectivity& faceCells  = *_faceCells;
  // determine the left and right cell zones of each face zone
  foreach(FaceZonesMap::value_type& pos, _faceZones)
  {
      FluentFaceZone& fz = *(pos.second);

      const int c0 = faceCells(fz.iBeg,0);
      fz.leftCellZoneId = getCellZoneID(c0);
      FluentCellZone *lcz = _cellZones[fz.leftCellZoneId];

      const int c1 = faceCells(fz.iBeg,1);
      fz.rightCellZoneId = getCellZoneID(c1);

      if (fz.rightCellZoneId == fz.leftCellZoneId)
      {
          lcz->interiorZoneIds.push_back(fz.ID);
      }
      else if (fz.rightCellZoneId > 0)
      {
          lcz->interfaceZoneIds.push_back(fz.ID);
          FluentCellZone *rcz = _cellZones[fz.rightCellZoneId];
          rcz->interfaceZoneIds.push_back(fz.ID);
      }
      else
        lcz->boundaryZoneIds.push_back(fz.ID);
  }
}
      
Mesh*
FluentReader::createMesh(const int cellZoneID)
{
  const CRConnectivity& faceCells  = *_faceCells;

  const Array<int>& fcRow = faceCells.getRow();
  const Array<int>& fcCol = faceCells.getCol();
  
  Mesh *mesh = new Mesh(_dimension);
  
  FluentCellZone& cz = *(_cellZones[cellZoneID]);


  vector<int> allFaceList;

  // determine interior faces
  int faceOffset = 0;
  foreach(int fzId, cz.interiorZoneIds)
  {
      const FluentFaceZone& fz = *(_faceZones[fzId]);
      for(int i=fz.iBeg; i<=fz.iEnd; i++) allFaceList.push_back(i);
  }

  faceOffset = allFaceList.size();
  mesh->createInteriorFaceGroup(faceOffset);

  vector<int> interfaceFaceList;
  foreach(int fzId, cz.interfaceZoneIds)
  {
      const FluentFaceZone& fz = *(_faceZones[fzId]);
      for(int i=fz.iBeg; i<=fz.iEnd; i++)
      {
          allFaceList.push_back(i);
          interfaceFaceList.push_back(i);
      }

      const int thisFZCount = fz.iEnd-fz.iBeg+1;
      mesh->createInterfaceGroup(thisFZCount,faceOffset,fzId);
      faceOffset += thisFZCount;
  }
  
  vector<int> boundaryCells;
  foreach(int fzId, cz.boundaryZoneIds)
  {
      const FluentFaceZone& fz = *_faceZones[fzId];
      for(int i=fz.iBeg; i<=fz.iEnd; i++) allFaceList.push_back(i);

      const int thisFZCount = fz.iEnd-fz.iBeg+1;

      mesh->createBoundaryFaceGroup(thisFZCount,faceOffset,fzId,fz.zoneType);
      faceOffset += thisFZCount;

      for(int j=fcRow[fz.iBeg]; j<fcRow[fz.iEnd+1]; j++)
        if (fcCol[j] > _numCells)
          boundaryCells.push_back(fcCol[j]);
  }
  
  mesh->getFaces().setCount(allFaceList.size());

  Array<int> allFaceArray(allFaceList.size());
  for(unsigned int i=0; i<allFaceList.size(); i++) allFaceArray[i] = allFaceList[i];
  
  shared_ptr<CRConnectivity> mFaceCells(_faceCells->getSubset(mesh->getFaces(),allFaceArray));
  shared_ptr<CRConnectivity> mFaceNodes(_faceNodes->getSubset(mesh->getFaces(),allFaceArray));
  

  int numMeshCells = cz.iEnd-cz.iBeg+1;

  Array<int> globalToLocalCellMap(_numCells+_numBoundaryFaces);
  globalToLocalCellMap = -1;

  vector<int> interfaceCellList;

  const int interfaceFaceCount = interfaceFaceList.size();
  if (interfaceFaceCount > 0)
  {
      Array<int> ifFacesArray(interfaceFaceCount);
      StorageSite ifFacesSite(interfaceFaceCount,0);
      StorageSite ifNodesSite(0,0);
      StorageSite ifCellsSite(0,0);

      shared_ptr<CRConnectivity> ifFaceNodes(_faceNodes->getLocalizedSubset(ifFacesSite,
                                                                            ifNodesSite,
                                                                            ifFacesArray));

      const Array<int>& interfaceNodes = ifFaceNodes->getLocalToGlobalMap();

      
      shared_ptr<CRConnectivity> ifNodeCells(getNodeCells().getLocalizedSubset(ifNodesSite,
                                                                               ifCellsSite,
                                                                               interfaceNodes));

      const Array<int>& interfaceAllCells = ifNodeCells->getLocalToGlobalMap();

      for(int i=0; i<interfaceAllCells.getLength(); i++)
      {
          const int c = interfaceAllCells[i];
          if ((c < cz.iBeg || c >cz.iEnd) &&
              (c < _numCells))
            interfaceCellList.push_back(c);
      }
  }

  const int numGhostCells = interfaceCellList.size();
  const int numBoundaryCells = boundaryCells.size();

  const int nTotalCells = numMeshCells+numGhostCells+numBoundaryCells;
  Array<int> allCellList(nTotalCells);
  
  int nc=0;
  for(int i=cz.iBeg; i<=cz.iEnd; i++)
  {
      allCellList[nc]=i;
      globalToLocalCellMap[i]=nc++;
  }
  
  foreach(int i, interfaceCellList)
  {
      allCellList[nc]=i;
      globalToLocalCellMap[i]=nc++;
  }

  foreach(int i, boundaryCells)
  {
      allCellList[nc]=i;
      globalToLocalCellMap[i]=nc++;
  }

  mesh->getCells().setCount(numMeshCells, numGhostCells+numBoundaryCells);

  StorageSite tempNodesSite(0,0);

  shared_ptr<CRConnectivity> czAllCellNodes(getCellNodes().getLocalizedSubset(mesh->getCells(),
                                                                              tempNodesSite,
                                                                              allCellList));
  shared_ptr<Array<Vec3> > coords = _coords.getSubset(czAllCellNodes->getLocalToGlobalMap());

  StorageSite& nodes = mesh->getNodes();
  nodes.setCount(coords->getLength());
  

  cz.mesh = mesh;

  foreach(const CellZonesMap::value_type& pos, _cellZones)
  {
      const FluentCellZone& ocz = *(pos.second);
      if (&ocz != &cz)
      {
          OneToOneIndexMap *im = getGhostCellMap(ocz,allCellList);
          if (im != 0)
            cz.ghostCellMaps[ocz.ID] = im;
      }
  }

  return mesh;
}

MeshList
FluentReader::getMeshList()
{
  MeshList meshes;

  foreach(const CellZonesMap::value_type& pos, _cellZones)
  {
      const FluentCellZone& cz = *(pos.second);
      meshes.push_back(createMesh(cz.ID));
  }

  foreach(const CellZonesMap::value_type& pos, _cellZones)
  {
      const FluentCellZone& cz = *(pos.second);
      Mesh *mesh = cz.mesh;

      foreach(const GhostCellMapsMap::value_type& pos2, cz.ghostCellMaps)
      {
          const FluentCellZone& ocz = *_cellZones[pos2.first];
          Mesh *omesh = ocz.mesh;
          mesh->getCells().getMappers()[&omesh->getCells()] = pos2.second;
      }
  }
  return meshes;
}

OneToOneIndexMap*
FluentReader::getGhostCellMap(const FluentCellZone& cz, const Array<int>& indices)
{
  const int iBeg = cz.iBeg;
  const int iEnd = cz.iEnd;

  int thisZoneCells=0;
  for(int ii=0; ii<indices.getLength(); ii++)
  {
      const int c = indices[ii];
      if (c >= iBeg && c<=iEnd) thisZoneCells++;
  }

  if (thisZoneCells == 0) return 0;
  
  Array<int> *fromPtr = new Array<int>(thisZoneCells);
  Array<int> *toPtr = new Array<int>(thisZoneCells);

  Array<int>& from = *fromPtr;
  Array<int>& to = *toPtr;

  thisZoneCells = 0;
  for(int ii=0; ii<indices.getLength(); ii++)
  {
      const int c = indices[ii];
      if (c >= iBeg && c<=iEnd)
      {
          to[thisZoneCells] = ii;
          from[thisZoneCells] = c-iBeg;
          thisZoneCells++;
      }
  }

  OneToOneIndexMap *imap = new OneToOneIndexMap(fromPtr,toPtr);
  return imap;
}


#if 0
               
PyRet
FluentReader::orderCellFacesAndNodes(PyArgsIn args)
{
  CRConnectivity& cellFaces = args.getRef<CRConnectivity>(0); 
  CRConnectivity& cellNodes = args.getRef<CRConnectivity>(1); 
  const CRConnectivity& faceNodes = args.getRef<CRConnectivity>(2); 
  const CRConnectivity& faceCells = args.getRef<CRConnectivity>(3);

  ::orderCellFacesAndNodes(cellFaces,cellNodes,faceNodes,faceCells);
  return PyReturn();
}

PyReturn
FluentReader::coupleCells(PyArgsIn args)
{
  const State fp = args.getState(0);
  const int count = fp.getInt("count");
  const Array<int>& leftFaces = fp.getChildRef<Array<int> >("leftFaces");
  const Array<int>& rightFaces = fp.getChildRef<Array<int> >("rightFaces");

  CRConnectivity& faceCells = getChildRef<CRConnectivity>("faceCells");

  for(int n=0; n<count; n++)
  {
      const int lf = leftFaces[n];
      const int rf = rightFaces[n];
      faceCells(lf,1) = faceCells(rf,0);
      faceCells(rf,1) = faceCells(lf,0);
  }
  return PyReturn();
}


PyReturn
FluentReader::_getCellFaces(PyArgsIn)
{
  CRConnectivity& faceCells = getChildRef<CRConnectivity>("faceCells");
  CRConnectivity& faceNodes = getChildRef<CRConnectivity>("faceNodes");

  CRConnectivity *cellFacesPtr = new CRConnectivity(faceCells.getColSite(),faceCells.getRowSite());
  CRConnectivity& cellFaces = *cellFacesPtr;
  
  cellFaces.initCount();

  const int nFaces = faceCells.getRowDim();

  for (int f=0; f<nFaces; f++)
  {
      if (faceNodes.getCount(f) > 0)
      {
          for(int nfc =0; nfc<faceCells.getCount(f); nfc++)
            cellFaces.addCount(faceCells(f,nfc),1);
      }
  }
  
  cellFaces.finishCount();

  for (int f=0; f<nFaces; f++)
  {
      if (faceNodes.getCount(f) > 0)
      {
          for(int nfc =0; nfc<faceCells.getCount(f); nfc++)
            cellFaces.add(faceCells(f,nfc),f);
      }
  } 
  cellFaces.finishAdd();
  return cellFacesPtr;
}
#endif
