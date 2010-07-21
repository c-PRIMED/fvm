#ifndef _SPIKEMATRIX_H_
#define _SPIKEMATRIX_H_

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

#include "Matrix.h"
#include "CRConnectivity.h"
#include "Array.h"
#include "Array2D.h"
#include "SpikeStorage.h"

template<typename T_Diag, typename T_OffDiag,  typename X>
class SpikeMatrix : public Matrix
{
public:

   //friend class LinearSystemMerger;

  typedef T_Diag Diag;
  typedef T_OffDiag OffDiag;
  typedef Array<Diag> DiagArray;
  typedef Array<OffDiag> OffDiagArray;
  typedef Array<X> XArray;
  
  SpikeMatrix( const CRConnectivity& conn,  const Array<T_Diag>& diag, const Array<T_OffDiag>& off_diag, const SpikeStorage& spike_storage ) :
   Matrix(),
   _conn(conn),
   _diag(diag),
   _offDiag(off_diag),
   _spikeStorage(spike_storage),
   _bandwidth(spike_storage.getBandWidth()),
   _ncells(conn.getRowSite().getSelfCount()),
   _A(2*_bandwidth+1,_ncells), 
   _LSPIKE(_ncells,_bandwidth),
   _RSPIKE(_ncells,_bandwidth)
  {
    initAssembly();
    logCtor();
  }


  virtual ~SpikeMatrix()
  {
    logDtor();
  }


private:
  void initAssembly()
  {
     setMatrix();
  }

  void setMatrix()
  {
     //forming A
     const int  nr = 2 * _bandwidth + 1;
     const int  nc = _ncells;
     //diagonal filling
     for ( int i = 0; i < nc; i++ )
         _A(_bandwidth,i) = _diag[i];
     const Array<int>& row = _conn.getRow();
     const Array<int>& col = _conn.getCol();
     for ( int i = 0; i < _ncells; i++ ){
         for (int n = row[i]; n < row[i+1]; n++ ){
	     int j =  col[n];
	     //check if it is in bandwidth and inner coefficient
	     if ( abs(j-i) <= _bandwidth && j < _ncells ){
	       _A(_bandwidth-(j-i),j) = _offDiag[ _conn(i,j) ];
	     }
         }
     }
     
     _A.print(cout);
  } 
  //void setLSpikeMtrx();
  /*void setRSpikeMtrx();*/


  const CRConnectivity& _conn;
  const Array<Diag>& _diag;
  const Array<OffDiag>& _offDiag;
  const SpikeStorage& _spikeStorage;
  const	int _bandwidth;
  const int _ncells;
  Array2D<Diag>  _A;
  Array2D<Diag>  _LSPIKE;
  Array2D<Diag>  _RSPIKE;
  
 
};


#endif
