// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _SpikeStorage_H_
#define _SpikeStorage_H_

#include <vector>
/**
 *Spike data storage
 * 
 */
class CRConnectivity;

class SpikeStorage 
{
public:
 //conn: cellcells connectivtiy 
  SpikeStorage(const CRConnectivity& conn, int semi_bandwidth);
  ~SpikeStorage();
  //interior
  const vector<int>&  getLSPKInterior() const {return _LSPK_INTERIOR;}  
  const vector<int>&  getRSPKInterior() const {return _RSPK_INTERIOR;}  
  //ghost 
  const vector<int>&  getLSPKGhost() const {return _LSPK_GHOST;}  
  const vector<int>&  getRSPKGhost() const {return _RSPK_GHOST;}  
  //offdiag_ptr
  const vector<int>&  getLSPKOffDiagPtr() const {return _LSPK_OFFD_PTR;}  
  const vector<int>&  getRSPKOffDiagPtr() const {return _RSPK_OFFD_PTR;} 
  //get I and J indices (local)
  const vector<int>&  getLSPKIndexI() const { return _LSPK_I;}
  const vector<int>&  getLSPKIndexJ() const { return _LSPK_J;}
  const vector<int>&  getRSPKIndexI() const { return _RSPK_I;}
  const vector<int>&  getRSPKIndexJ() const { return _RSPK_J;}
  //get GhostCount
  const vector<int>&  getLSPKCountGhost() const { return _LSPKCountGhost;}
  const vector<int>&  getRSPKCountGhost() const { return _RSPKCountGhost;}

  int getBandWidth() const { return _bandwidth;}

  
private:
   void init();
   void gatherCellSizes();
   void syncCellIDs();
   void setGlobalIndices();
   void setOffDiagPtr();

   const CRConnectivity& _conn;
   int  _bandwidth;

   int _procID;
   int _localCellSelfCount;
   map<int,int>    _ghostMap; //map before syncLocal and after syncLocal
   
   vector<int>   _cellSelfCounts; //size=comm.size(), store each partition cell size 
   vector<int>  _LSPK_INTERIOR;
   vector<int>  _LSPK_GHOST;
   vector<int>  _LSPK_OFFD_PTR;
   vector<int>  _RSPK_INTERIOR;
   vector<int>  _RSPK_GHOST;
   vector<int>  _RSPK_OFFD_PTR;
   vector<int>  _LSPK_I;
   vector<int>  _LSPK_J;
   vector<int>  _RSPK_I;
   vector<int>  _RSPK_J;
   vector<int>  _glblIndices;
   vector<int>  _LSPKCountGhost;
   vector<int>  _RSPKCountGhost;

};

#endif

