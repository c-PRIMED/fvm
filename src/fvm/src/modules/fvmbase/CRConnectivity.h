// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _CRCONNECTIVITY_H_
#define _CRCONNECTIVITY_H_

#include "misc.h"
#include "Array.h"
#include "Vector.h"
#include "StorageSite.h"

/**
 * Stores connectivity information using a compressed row
 * format. Basically this means that we have a (integer) row array
 * that stores the indices into col array (also integer) where the non
 * zero entries for a particular row begin. Thus col[row[i]] through
 * col[row[i+1]] store the column indices j where a sparse matrix
 * based on this connectivity would have non zero entries.
 *
 * Objects of this class can be used to store connectivities between
 * mesh elements such as the faces of a cell or the cells surrounding
 * a node etc. It can also be used to represent the fill pattern of an
 * arbitrary sparse matrix (see CRMatrix).
 *
 * This class provides support for constructing the connectivity in a
 * two step fashion but is not usable for dynamically changing the
 * connectivity. Typical usage is
 *
 * -- construct by specifying rowSite and colSite. Thus the dimensions
 * are known but not the number of nonzero entries
 *
 * -- call initCount to start the process of counting non zeros
 *
 * -- call addCount to specify the number of non zeros for each row i
 * (this can be in any order and multiple calls for same row are
 * allowed)
 *
 * -- call finishCount - this allocates the col array and resets row
 * array so that it is ready for entries to be added
 *
 * -- call add to specify the nonzero location j for each row i (this
 * can also be in any order but each i,j must be specified only once)
 *
 * Since we treat this connectivity as a matrix we can perform
 * operations like transpose (to get node - cell connectvity from cell
 * - node connectivity, for example) and multiplication (multiplying
 * cell - faces with face - nodes to get cell nodes, for example)
 * 
 */

class CRConnectivity
{
public:
  typedef Array<Vector<int,2> > PairToColMapping;

  enum CRTYPE
    {
      CELLCELL1 = 1,
      CELLCELL2 = 2
    };
  
  CRConnectivity(const StorageSite& rowSite, const StorageSite& colSite);
  ~CRConnectivity();

  DEFINE_TYPENAME("CRConnectivity");
  
  int getRowDim() const {return _rowDim;}
  int getColDim() const {return _colDim;}

  const StorageSite& getRowSite() const {return *_rowSite;}
  const StorageSite& getColSite() const {return *_colSite;}

  shared_ptr<CRConnectivity> getTranspose() const;
  shared_ptr<CRConnectivity> getMultiTranspose(const int varSize) const;

  shared_ptr<CRConnectivity> multiply(const CRConnectivity& b,
                                      const bool implicitDiagonal) const;
  
  shared_ptr<CRConnectivity>
  createOffset(const StorageSite& newRowSite,
               const int offset, const int size) const;
    
  shared_ptr<CRConnectivity>
  getSubset(const StorageSite& site,const Array<int>& indices) const;

  void reorder( const Array<int>& indices );
  
  shared_ptr<CRConnectivity>
  getLocalizedSubset(const StorageSite& newRowSite,
                     StorageSite& newColSite,
                     const Array<int>& indices) const;

  shared_ptr<CRConnectivity>
  getLocalizedSubsetOfFaceCells(const StorageSite& newRowSite,  StorageSite& newColSite,
                                const Array<int>& indices, const CRConnectivity& faceCells, 
                                const CRConnectivity& cellCells) const;
  //this is more robust method  
  shared_ptr<CRConnectivity>
  getLocalizedSubsetOfFaceCells(const StorageSite& newRowSite,  StorageSite& newColSite,
                                const Array<int>& indices, const CRConnectivity& cellParts, const int partID) const;


  void
  localize(const Array<int>& globalToLocal,const StorageSite& newColSite);
 
  //@{
  
  /**
   * see comments above
   * 
   */

  void initCount();

  void addCount(const int index, const int count)
  {
    (*_row)[index] += count;
  }

  void finishCount();

  int add(const int index, const int val)
  {
    int pos = (*_row)[index];
    (*_col)[pos] = val;
    (*_row)[index]++;
    return pos;
  }

  void finishAdd();

  //@}
  
  int operator()(const int i, const int j) const
  {
    return (*_col)[(*_row)[i]+j];
  }

  /**
   * returns the index of the j'th non zero column for row i
   * 
   */

  int& operator()(const int i, const int j)
  {
    return (*_col)[(*_row)[i]+j];
  }

  int getCoeffPosition(const int i,  const int j) const
  {
    for (int nnb = (*_row)[i]; nnb<(*_row)[i+1]; nnb++)
    {
        if ((*_col)[nnb] == j)
          return nnb;
    }
    throw CException("invalid indices");
  }

  /**
   * returns the number of non zeroes for row i
   * 
   */

  int getCount(const int i) const
  {
    return ((*_row)[i+1] - (*_row)[i]);
  }

  //@{

  /**
   * direct access to row and col array
   * 
   */

  void resizeLocalToGlobalMap( int size ) 
  {
     _localToGlobalMap = shared_ptr< Array<int> > ( new Array<int>(size) );
  }

  const Array<int>& getRow() const {return *_row;}
  const Array<int>& getCol() const {return *_col;}
  Array<int>& getCol() {return *_col;}

  //@}


  const Array<int>& getGlobalToLocalMap() const {return *_globalToLocalMap;}
  const Array<int>& getLocalToGlobalMap() const {return *_localToGlobalMap;}

  const map<int,int>& getGlobalToLocalMapper() const {return _globalToLocalMapper;}
  map<int,int>& getGlobalToLocalMapper()           {return _globalToLocalMapper;}
  
  shared_ptr<Array<int> > getGlobalToLocalMapPtr() {return _globalToLocalMap;}
  shared_ptr<Array<int> > getLocalToGlobalMapPtr() {return _localToGlobalMap;}

  const PairToColMapping&
  getPairToColMapping(const CRConnectivity& pairs) const;
  
  void
  clearPairToColMapping(const CRConnectivity& pairs) const;

  void setConnType( CRTYPE type ){
      _connType = type;
  }

  const CRTYPE& getConnType() const { 
      return _connType;
  }
  

private:
  CRConnectivity(const CRConnectivity&);
  mutable StorageSite const *_rowSite;
  mutable StorageSite const *_colSite;
  mutable int _rowDim;
  mutable int _colDim;
  shared_ptr<Array<int> > _row;
  shared_ptr<Array<int> > _col;
  shared_ptr<Array<int> > _globalToLocalMap;
  map<int,int>            _globalToLocalMapper;
  shared_ptr<Array<int> > _localToGlobalMap;
  mutable map<const CRConnectivity*, PairToColMapping*> _pairToColMappings;
  CRTYPE  _connType;
};

#endif
