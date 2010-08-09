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
   _L(_bandwidth,_bandwidth),
   _R(_bandwidth,_bandwidth),
   _LSpikeT(_bandwidth,_bandwidth),
   _RSpikeB(_bandwidth,_bandwidth)
  {
    initAssembly();
    logCtor();
  }


  virtual ~SpikeMatrix()
  {
    logDtor();
  }
  //initial assamely 
  void initAssembly()
  {
     setMatrix();
     setLSpikeMtrx();
     setRSpikeMtrx();
  }
  //solve 
  void lu()
  {
     const int b = _bandwidth;
     for ( int i = 0; i < _ncells-1; i++ ){
        const Diag pivot = _A(b,i);
	for ( int j = i+1;  j <= min(_ncells-1,i+b); j++ ){
	   const int j2 = b+j-i;
	   const Diag m = _A(j2,i) / pivot; //Division(/) should be defined as _A * pivot^-1 for matrix ops.
	   _A(j2,i) = m;
	   for ( int k = i+1; k <= min(_ncells-1,i+b); k++ ){
	      const int j2 = b+j-k;
	      const int i2 = b+i-k;
	      _A(j2,k) -= m * _A(i2,k);
           }
       }
    }  
    _A.print(cout);
  }



private:

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
  //left spike matrix
  void setLSpikeMtrx()
  {
     const vector<int>& vecI = _spikeStorage.getLSPKIndexI();
     const vector<int>& vecJ = _spikeStorage.getLSPKIndexJ();
     const vector<int>& offDiagPtr = _spikeStorage.getLSPKOffDiagPtr(); 
     const vector<int>& countGhost = _spikeStorage.getLSPKCountGhost();
     int indx = 0;
     for ( int n = 0; n < _bandwidth; n++ ){
        int ncount = countGhost[n];
	for ( int c = 0; c < ncount; c++){ 
          int i = vecI[indx];
	  int j = vecJ[indx];
	  _L(i,j) = _offDiag[offDiagPtr[indx++]];
	}
     }
     _L.print(cout);
  }
  //right spike matrix
  void setRSpikeMtrx()
  {
     const vector<int>& vecI = _spikeStorage.getRSPKIndexI();
     const vector<int>& vecJ = _spikeStorage.getRSPKIndexJ();
     const vector<int>& offDiagPtr = _spikeStorage.getRSPKOffDiagPtr();
     const vector<int>& countGhost = _spikeStorage.getRSPKCountGhost();
     int indx = 0;
     for ( int n = _ncells-_bandwidth; n < _ncells; n++ ){
        int ncount = countGhost[n];
	for ( int c = 0; c < ncount; c++){ 
          int i = _bandwidth - _ncells +  vecI[indx];
	  int j = vecJ[indx];
	  _R(i,j) = _offDiag[offDiagPtr[indx++]];
	}
     }
     _R.print(cout);
  }
   

  const CRConnectivity& _conn;
  const Array<Diag>& _diag;
  const Array<OffDiag>& _offDiag;
  const SpikeStorage& _spikeStorage;
  const	int _bandwidth;
  const int _ncells;
  Array2D<Diag>  _A;
  Array2D<Diag>  _L;     
  Array2D<Diag>  _R;
  Array2D<Diag>  _LSpikeT; //left-top spike matrix
  Array2D<Diag>  _RSpikeB; //right-top spike matrix
 
};


#endif
