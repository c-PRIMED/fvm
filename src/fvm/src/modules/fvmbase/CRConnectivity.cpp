// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "CRConnectivity.h"
#include "CException.h"
#include "StorageSite.h"
#include <map>
#include <set>

CRConnectivity::CRConnectivity(const StorageSite& rowSite,
                               const StorageSite& colSite) :
  _rowSite(&rowSite),
  _colSite(&colSite),
  _rowDim(_rowSite->getCountLevel1()),
  _colDim(_colSite->getCountLevel1()),
  _row(),
  _col(),
  _globalToLocalMap(),
  _localToGlobalMap(),
  _pairToColMappings(),
  _connType(CELLCELL1)
{
  logCtorVerbose("from rowSite and colSite of dimensions (%dx%d)" ,
                 _rowDim,_colDim);
}


CRConnectivity::~CRConnectivity()
{
  for(map<const CRConnectivity*,PairToColMapping*>::iterator p
        = _pairToColMappings.begin();
      p != _pairToColMappings.end();
      ++p)
    delete p->second;

  logDtorVerbose("of dimensions (%dx%d)" , _rowDim,_colDim);
}


void CRConnectivity::initCount()
{
  _row = shared_ptr<Array<int> >(new Array<int>(_rowDim+1));
  *_row = 0;
}

void CRConnectivity::finishCount()
{
  Array<int>& row = *_row;

  for(int i=0; i<_rowDim; i++)
    row[i+1] += row[i];
  
  const int colSize = (_rowDim == 0) ? 0 : row[_rowDim-1];

    for(int i=_rowDim; i>0; i--)
    row[i] = row[i-1];
  row[0] = 0;

  _col = shared_ptr<Array<int> >(new Array<int>(colSize));
}


void CRConnectivity::finishAdd()
{
  Array<int>& row = *_row;
  for(int i=_rowDim; i>0; i--)
    row[i] = row[i-1];
  row[0] = 0;
}


shared_ptr<CRConnectivity>
CRConnectivity::getTranspose() const
{
  shared_ptr<CRConnectivity> trPtr(new CRConnectivity(*_colSite,*_rowSite));
  CRConnectivity& tr = *trPtr;
  
  tr.initCount();

  const Array<int>& myRow = *_row;
  const Array<int>& myCol = *_col;
  for(int i=0; i<_rowDim; i++){
    for(int j=myRow[i]; j<myRow[i+1]; j++)
      tr.addCount(myCol[j],1);
  }
  tr.finishCount();

  for(int i=0; i<_rowDim; i++)
    for(int j=myRow[i]; j<myRow[i+1]; j++)
      tr.add(myCol[j],i);
  tr.finishAdd();
  return trPtr;
}

// returns a transpose connectivity for a "flattened" matrix, ie. a
// connectivity matrix where each row is made of varSize number of
// separate rows

shared_ptr<CRConnectivity>
CRConnectivity::getMultiTranspose(const int varSize) const
{
  shared_ptr<CRConnectivity> trPtr(new CRConnectivity(*_colSite,*_rowSite));
  CRConnectivity& tr = *trPtr;

  tr._rowDim *= varSize;
  tr._colDim *= varSize;
  
  tr.initCount();

  const Array<int>& myRow = *_row;
  const Array<int>& myCol = *_col;

  for(int i=0; i<_rowDim; i++)
  {
      for(int ndr=0; ndr<varSize; ndr++)
      {
          const int nfr = i*varSize + ndr;
          tr.addCount(nfr, varSize);
      }
      
      for(int jp=myRow[i]; jp<myRow[i+1]; jp++)
      {
          const int j = myCol[jp];

          for(int ndr=0; ndr<varSize; ndr++)
            for(int ndc=0; ndc<varSize; ndc++)
            {
                //const int nfr = i*varSize + ndr;
                const int nfc = j*varSize + ndc;
                tr.addCount(nfc, 1);
            }
      }
  }
  
  tr.finishCount();

  for(int i=0; i<_rowDim; i++)
  {
      for(int ndr=0; ndr<varSize; ndr++)
        for(int ndc=0; ndc<varSize; ndc++)
        {
            const int nfr = i*varSize + ndr;
            const int nfc = i*varSize + ndc;
            tr.add(nfc, nfr);
        }
      
      for(int jp=myRow[i]; jp<myRow[i+1]; jp++)
      {
          const int j = myCol[jp];

          for(int ndr=0; ndr<varSize; ndr++)
            for(int ndc=0; ndc<varSize; ndc++)
            {
                const int nfr = i*varSize + ndr;
                const int nfc = j*varSize + ndc;
                tr.add(nfc, nfr);
            }
      }
  }

  tr.finishAdd();
  return trPtr;
}


shared_ptr<CRConnectivity>
CRConnectivity::multiply(const CRConnectivity& b, const bool implicitDiagonal) const
{
  const bool isSquared = (&b == this);
  
  if (_colDim != b._rowDim)
    cerr << "invalid connectivity multiplication" << endl;

  shared_ptr<CRConnectivity> prPtr(new CRConnectivity(*_rowSite,*b._colSite));
  CRConnectivity& pr = *prPtr;

  pr.initCount();

  const Array<int>& myRow = *_row;
  const Array<int>& myCol = *_col;


  const Array<int>& bRow = *b._row;
  const Array<int>& bCol = *b._col;

  Array<bool> marker(b._colDim);
  Array<int> marked(b._colDim);

  for(int i=0; i<b._colDim; i++)
  {
      marker[i] = false;
      marked[i] = 0;
  }
  for(int i=0; i<_rowDim; i++)
  {
      int nMarked = 0;
      for(int ir = myRow[i]; ir<myRow[i+1]; ir++)
      {
          const int ja = myCol[ir];

          if (isSquared)
          {
              if (!marker[ja])
              {
                  marker[ja] = true;
                  if (ja  != i || !implicitDiagonal)
                    pr.addCount(i,1);
                  marked[nMarked++] = ja;
              }
              
          }
          
          for (int rb = bRow[ja]; rb<bRow[ja+1];  rb++)
          {
              const int jb = bCol[rb];
              if (!marker[jb])
              {
                  marker[jb] = true;
                  if (jb  != i || !implicitDiagonal)
                    pr.addCount(i,1);
                  marked[nMarked++] = jb;
              }
          }
      }
      for(int n=0; n<nMarked; n++)
        marker[marked[n]] = false;
  }

  // for(int i=0; i<pr._rowDim;i++)
  //  cerr << i << " " << (*pr._row)[i] << endl;
  
  pr.finishCount();

  for(int i=0; i<_rowDim; i++)
  {
      int nMarked = 0;
      for(int ir = myRow[i]; ir<myRow[i+1]; ir++)
      {
          const int ja = myCol[ir];
          if (isSquared)
          {
              if (!marker[ja])
              {
                  marker[ja] = true;
                  if (ja  != i || !implicitDiagonal)
                    pr.add(i,ja);
                  marked[nMarked++] = ja;
              }
          }
          
          for (int rb = bRow[ja]; rb<bRow[ja+1];  rb++)
          {
              const int jb = bCol[rb];
              if (!marker[jb])
              {
                  marker[jb] = true;
                  if (jb  != i || !implicitDiagonal)
                    pr.add(i,jb);
                  marked[nMarked++] = jb;
              }
          }
      }
      for(int n=0; n<nMarked; n++)
        marker[marked[n]] = false;
  }

  pr.finishAdd();
  return prPtr;
}


shared_ptr<CRConnectivity>
CRConnectivity::getSubset(const StorageSite& site,
                          const Array<int>& indices) const
{
  shared_ptr<CRConnectivity> subPtr(new CRConnectivity(site,*_colSite));

  CRConnectivity& sub = *subPtr;
  const Array<int>& myRow = *_row;
  const Array<int>& myCol = *_col;
  
  sub.initCount();

  for(int ii=0;ii<indices.getLength();ii++)
  {
      const int i = indices[ii];
      sub.addCount(ii,getCount(i));
  }

  sub.finishCount();
  for(int ii=0;ii<indices.getLength();ii++)
  {
      const int i = indices[ii];
      for(int ip=myRow[i];ip<myRow[i+1]; ip++)
      {
          const int j = myCol[ip];
          sub.add(ii,j);
      }
  }

  sub.finishAdd();
  return subPtr;
}



//this is  specialized subset for PartMesh class, V type faceCells siutation for interface ghost cells sharing common global cell  result
//in wrong number of local element counts and mapping (Gazi).
shared_ptr<CRConnectivity>
CRConnectivity::getLocalizedSubsetOfFaceCells( const StorageSite& newRowSite,  StorageSite& newColSite,
                                               const Array<int>& indices, const CRConnectivity& faceCells, const CRConnectivity& cellCells ) const
{
  const Array<int>& myRow = *_row;
  const Array<int>& myCol = *_col;
 
  map<int, int>  localToGlobalMap;
  map<int, vector<int> > faceToGlobalCellsMap; 
  map<int, vector<int> > faceToLocalCellsMap;

  shared_ptr<Array<int> > globalToLocalPtr(new Array<int>(_colDim));

  Array<int>& globalToLocal = *globalToLocalPtr;

  globalToLocal = -1;
 
   int max_sur_cells = 0;
   for ( int elem = 0; elem < cellCells.getRowDim(); elem++ )
        max_sur_cells = std::max( max_sur_cells, cellCells.getCount(elem) );


  int nLocal=0;
//loop over inner elements
  for(int ii=0;ii<indices.getLength();ii++)
  {
      bool inner_face  = true;
      const int i = indices[ii];
      vector<int> localConn ( myRow[i+1] - myRow[i], -1 );
      vector<int> globalConn( myRow[i+1] - myRow[i], -1 );
      int indx = 0;
      for(int ip=myRow[i];ip<myRow[i+1]; ip++)
      {
          const int j = myCol[ip];
          const int local_elem = faceCells.getGlobalToLocalMap()[j];
          globalConn.at(indx) = j;

         if ( cellCells.getCount(local_elem) != max_sur_cells )
               inner_face = false;
 
         if ( cellCells.getCount(local_elem)  == max_sur_cells  ){ //if it is a inner element
            if( globalToLocal[j] == -1  ){
             localToGlobalMap.insert( pair<int,int>(nLocal, j) ); 
             localConn.at(indx)  = nLocal;
             globalToLocal[j] = nLocal++;
            } else {
             localConn.at(indx) = globalToLocal[j];
            }
         } 
          indx++;

      }

      if ( inner_face ){ 
         faceToGlobalCellsMap[ii] = globalConn;
         faceToLocalCellsMap[ii]  = localConn;
      }
  }


  //outer faces (interface  + boundary)
  for(int ii=0;ii<indices.getLength();ii++)
  {
      bool outer_face = false;
      const int i = indices[ii];
      vector<int> localConn ( myRow[i+1] - myRow[i], -1 );
      vector<int> globalConn( myRow[i+1] - myRow[i], -1 );
      int indx = 0;
      for(int ip=myRow[i];ip<myRow[i+1]; ip++)
      {
          const int j = myCol[ip];
          const int local_elem = faceCells.getGlobalToLocalMap()[j];
           globalConn.at(indx) = j;
 
          if ( cellCells.getCount(local_elem) != max_sur_cells )
               outer_face = true;

          if ( cellCells.getCount(local_elem) != max_sur_cells ){ //if it is a boundary or interface cells
              localToGlobalMap.insert( pair<int,int>(nLocal, j) );
              localConn.at(indx)  = nLocal;
              globalToLocal[j]    = nLocal++;  //we still use globalToLocal for only inner and boundary (one-to-one mapping)
          } else {
               assert( globalToLocal[j] != -1 );
               localConn.at(indx) = globalToLocal[j];
          }
          indx++;

      }

      if ( outer_face ){ 
         faceToGlobalCellsMap[ii] = globalConn;
         faceToLocalCellsMap[ii]  = localConn;
      }

  }

  shared_ptr<Array<int> > localToGlobalPtr(new Array<int>(nLocal));

  Array<int>& localToGlobal = *localToGlobalPtr;

  newColSite.setCount(nLocal);
  shared_ptr<CRConnectivity> subPtr(new CRConnectivity(newRowSite,newColSite));
  CRConnectivity& sub = *subPtr;
  
  sub.initCount();

  for(int ii=0;ii<indices.getLength();ii++)
  {
      const int i = indices[ii];
      sub.addCount(ii,getCount(i));
  }

  sub.finishCount();
  

  map<int,int>::const_iterator it;
  for ( it = localToGlobalMap.begin(); it != localToGlobalMap.end(); it++ ){
     int local_id  = it->first;
     int global_id = it->second;
     localToGlobal[ local_id ] =  global_id;
   }

 for(int ii=0;ii<indices.getLength();ii++)
  {
      const int i = indices[ii];
      int indx = 0;
      for(int ip=myRow[i];ip<myRow[i+1]; ip++)
      {
          const int jLocal = faceToLocalCellsMap[ii].at(indx);
          sub.add(ii,jLocal);
          indx++;
      }
  }


  sub.finishAdd();
  sub._globalToLocalMap=globalToLocalPtr;
  sub._localToGlobalMap=localToGlobalPtr;

  return subPtr;
}


//this is  specialized subset for PartMesh class, V type faceCells siutation for interface ghost cells sharing common global cell  result
//in wrong number of local element counts and mapping (Gazi). This is more robust version then
//above implementation
shared_ptr<CRConnectivity>
CRConnectivity::getLocalizedSubsetOfFaceCells( const StorageSite& newRowSite,  StorageSite& newColSite,
                                               const Array<int>& indices,  const CRConnectivity& cellParts, const int partID ) const
{
  const Array<int>& myRow = *_row;
  const Array<int>& myCol = *_col;
 
  map<int, int>  localToGlobalMap;
  map<int, vector<int> > faceToLocalCellsMap;

  const StorageSite& rowSite = cellParts.getRowSite();
  const int cellSelfCount = rowSite.getSelfCount();
  const int cellCount     = rowSite.getCount();
  shared_ptr<Array<int> > globalToLocalPtr(new Array<int>(cellCount));
  Array<int>& globalToLocal = *globalToLocalPtr;

  globalToLocal = -1;
 
  int nLocal=0;
//loop over faces (sub)
  for(int ii=0;ii<indices.getLength();ii++)
  {
      bool inner_face  = true;
      const int i = indices[ii];
      //check if it is interface, pick one cell partition id, then compare connected cells if they are anot equal
      const int compID =  cellParts( myCol[myRow[i]], 0 );
      for(int ip=myRow[i];ip<myRow[i+1]; ip++){
          const int j = myCol[ip];
          if ( cellParts(j,0) != compID ) //interface if it is true
             inner_face = false;
      }

      //first check if it is boundary cell, then it is boundary interface
      for(int ip=myRow[i];ip<myRow[i+1]; ip++){
          const int j = myCol[ip];
          if ( j  >= cellSelfCount ) //either interface or boundary 
             inner_face = false;
      }

      vector<int> localConn ( myRow[i+1] - myRow[i], -1 );
      //looping over cell around faces
      if ( inner_face ){ //if it is a inner element
          int indx = 0;
          for(int ip=myRow[i];ip<myRow[i+1]; ip++)
          {
             const int j = myCol[ip];
             if( globalToLocal[j] == -1  ){
                localToGlobalMap.insert( pair<int,int>(nLocal, j) ); 
                localConn.at(indx)  = nLocal;
                globalToLocal[j] = nLocal++;
             } else {
                localConn.at(indx) = globalToLocal[j];
             }
             indx++;
         } 
      }
      if ( inner_face ){ 
         faceToLocalCellsMap[ii]  = localConn;
      }
  }


  //now there might be cells which is inner cell but surrounded by all boundary/interfaces faces, we have to renumber it before ghost cells
  for(int ii=0;ii<indices.getLength();ii++){
      const int i = indices[ii];
      for(int ip=myRow[i];ip<myRow[i+1]; ip++){
          const int j = myCol[ip];
          // these cells should be inner cell (j<..), never visited (glb..=-1), and staying in current partition (cellParts... == partID)
          if ( j  < cellSelfCount && globalToLocal[j] == -1 &&  cellParts(j,0) == partID ){
             localToGlobalMap.insert( pair<int,int>(nLocal, j) ); 
             globalToLocal[j] = nLocal++;
          }
      }
  }


  //outer faces (interface  + boundary)
  for(int ii=0;ii<indices.getLength();ii++)
  {
      bool outer_face = false;
      const int i = indices[ii];
      //check if it is interface, pick one cell partition id, then compare connected cells if they are anot equal
      const int compID =  cellParts( myCol[myRow[i]], 0 );
      for(int ip=myRow[i];ip<myRow[i+1]; ip++){
          const int j = myCol[ip];
          if ( cellParts(j,0) != compID ) //either interface or boundary 
             outer_face = true;
      }

      //first check if it is boundary cell, then it is boundary interface
      for(int ip=myRow[i];ip<myRow[i+1]; ip++){
          const int j = myCol[ip];
          if ( j  >= cellSelfCount ) //either interface or boundary 
             outer_face = true;
      }
      vector<int> localConn ( myRow[i+1] - myRow[i], -1 );
      if ( outer_face ){ //if it is a boundary or interface cells          
         int indx = 0;
         for(int ip=myRow[i];ip<myRow[i+1]; ip++){
             const int j = myCol[ip];
             if ( globalToLocal[j] == -1 ){
               localToGlobalMap.insert( pair<int,int>(nLocal, j) );
               localConn.at(indx) = nLocal++;
             } else { 
               localConn.at(indx) = globalToLocal[j];
             }
             indx++;
         }
      }

      if ( outer_face ){ 
         faceToLocalCellsMap[ii]  = localConn;
      }

  }

  shared_ptr<Array<int> > localToGlobalPtr(new Array<int>(nLocal));

  Array<int>& localToGlobal = *localToGlobalPtr;

  newColSite.setCount(nLocal);
  shared_ptr<CRConnectivity> subPtr(new CRConnectivity(newRowSite,newColSite));
  CRConnectivity& sub = *subPtr;
  
  sub.initCount();

  for(int ii=0;ii<indices.getLength();ii++)
  {
      const int i = indices[ii];
      sub.addCount(ii,getCount(i));
  }

  sub.finishCount();
  

  map<int,int>::const_iterator it;
  for ( it = localToGlobalMap.begin(); it != localToGlobalMap.end(); it++ ){
     int local_id  = it->first;
     int global_id = it->second;
     localToGlobal[ local_id ] =  global_id;
   }

 for(int ii=0;ii<indices.getLength();ii++)
  {
      const int i = indices[ii];
      int indx = 0;
      for(int ip=myRow[i];ip<myRow[i+1]; ip++)
      {
          const int jLocal = faceToLocalCellsMap[ii].at(indx);
          sub.add(ii,jLocal);
          indx++;
      }
  }

  sub.finishAdd();
  sub._globalToLocalMap=globalToLocalPtr;
  sub._localToGlobalMap=localToGlobalPtr;
  return subPtr;
}





shared_ptr<CRConnectivity>
CRConnectivity::getLocalizedSubset(const StorageSite& newRowSite,
                                   StorageSite& newColSite,
                                   const Array<int>& indices) const
{
  const Array<int>& myRow = *_row;
  const Array<int>& myCol = *_col;
  
  shared_ptr<Array<int> > globalToLocalPtr(new Array<int>(_colDim));

  Array<int>& globalToLocal = *globalToLocalPtr;

  globalToLocal = -1;

  int nLocal=0;
  for(int ii=0;ii<indices.getLength();ii++)
  {
      const int i = indices[ii];
      for(int ip=myRow[i];ip<myRow[i+1]; ip++)
      {
          const int j = myCol[ip];
          if (globalToLocal[j] == -1)
            globalToLocal[j] = nLocal++;
      }
  }

  shared_ptr<Array<int> > localToGlobalPtr(new Array<int>(nLocal));

  Array<int>& localToGlobal = *localToGlobalPtr;

  newColSite.setCount(nLocal);
  shared_ptr<CRConnectivity> subPtr(new CRConnectivity(newRowSite,newColSite));
  CRConnectivity& sub = *subPtr;
  
  sub.initCount();

  for(int ii=0;ii<indices.getLength();ii++)
  {
      const int i = indices[ii];
      sub.addCount(ii,getCount(i));
  }

  sub.finishCount();
  for(int ii=0;ii<indices.getLength();ii++)
  {
      const int i = indices[ii];
      for(int ip=myRow[i];ip<myRow[i+1]; ip++)
      {
          const int j = myCol[ip];
          const int jLocal = globalToLocal[j];
          sub.add(ii,jLocal);
      }
  }

  for(int i=0; i<_colDim; i++)
    if (globalToLocal[i] != -1)
      localToGlobal[globalToLocal[i]] = i;
  
  sub.finishAdd();
  sub._globalToLocalMap=globalToLocalPtr;
  sub._localToGlobalMap=localToGlobalPtr;

  return subPtr;
}


void
CRConnectivity::reorder( const Array<int>& indices ) 
{
  Array<int>& myCol       = *_col;

  for ( int i = 0; i < myCol.getLength(); i++ )
       myCol[i]   = indices[ myCol[i] ];
}


void
CRConnectivity::localize(  const Array<int>& globalToLocal,
                           const StorageSite& newColSite)
{
  const int newColDim = newColSite.getCount();

  if (globalToLocal.getLength() != _colDim)
  {
      ostringstream e;
      e << "global to local mapping is of wrong size ("  <<
        globalToLocal.getLength() << "); should be " << _colDim;
      throw CException(e.str());
  }

  Array<int>& myCol = *_col;
  for(int i=0; i<myCol.getLength(); i++)
  {
      const int localIndex = globalToLocal[myCol[i]];
      if (localIndex >=0 && localIndex < newColDim)
        myCol[i] = localIndex;
      else
      {
          ostringstream e;
          e << "invalid global to local mapping at index "  <<
            myCol[i];
          throw CException(e.str());
      }
  }
  _colDim = newColDim;
  _colSite = &newColSite;
}


const Array<Vector<int,2> >&
CRConnectivity::getPairToColMapping(const CRConnectivity& pairs) const
{
  if (_pairToColMappings.find(&pairs) == _pairToColMappings.end())
  {
      const Array<int>& crRow = getRow();
      const Array<int>& crCol = getCol();
  
      const int pairCount = pairs.getRowDim();

      Array<Vector<int,2> > *pairToCol = new Array<Vector<int,2> >(pairCount);

      _pairToColMappings[&pairs] = pairToCol;
  
      for(int np=0; np<pairCount; np++)
      {
          if (pairs.getCount(np) == 2)
          {
              int n0 = pairs(np,0);
              int n1 = pairs(np,1);

              bool found = false;
              for(int r0=crRow[n0]; r0<crRow[n0+1]; r0++)
                if (crCol[r0] == n1)
                {
                    (*pairToCol)[np][0] = r0;
                    found = true;
                    break;
                }

              ostringstream e;
              if (!found)
              {
                  e << "did not find col pos for ( " << n0 << "," << n1 << ")";
                  throw CException(e.str());
              }

              found = false;
              for(int r1=crRow[n1]; r1<crRow[n1+1]; r1++)
                if (crCol[r1] == n0)
                {
                    (*pairToCol)[np][1] = r1;
                    found = true;
                    break;
                }

              if (!found)
              {
                  e << "did not find col pos for ( " << n1 << "," << n0 << ")";
                  throw CException(e.str());
              }
          
          }
          else
          {
              ostringstream e;
              e << "invalid count for pairs connectivity" << endl;
              throw CException(e.str());
          }
      }      
  }
  
  return *_pairToColMappings[&pairs];
}


shared_ptr<CRConnectivity>
CRConnectivity::createOffset(const StorageSite& newRowSite,
                             const int offset, const int size) const
{
  shared_ptr<CRConnectivity> offPtr(new CRConnectivity(newRowSite,*_colSite));
  shared_ptr<ArrayBase> rowOffset(_row->createOffsetArray(offset,size+1));
  offPtr->_row = dynamic_pointer_cast<Array<int> >(rowOffset);
  offPtr->_col = _col;
  return offPtr;
}

void
CRConnectivity::clearPairToColMapping(const CRConnectivity& pairs) const
{
  if (_pairToColMappings.find(&pairs) != _pairToColMappings.end())
  {
      Array<Vector<int,2> > *pairToCol = _pairToColMappings[&pairs];
      delete pairToCol;
      
      _pairToColMappings.erase(&pairs);
  }
}
