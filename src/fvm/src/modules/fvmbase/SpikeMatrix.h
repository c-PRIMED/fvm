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
   _procID(0),
   _nprocs(1),
   _A(2*_bandwidth+1,_ncells), 
   _LL(_ncells, _bandwidth),
   _L(_bandwidth,_bandwidth),
   _RR(_ncells, _bandwidth),
   _R(_bandwidth,_bandwidth),
   _LSpike (_ncells,_bandwidth),
   _LSpikeT(_bandwidth,_bandwidth),
   _RSpike (_ncells,_bandwidth),
   _RSpikeB(_bandwidth,_bandwidth),
   _JokerSpikeT(_bandwidth,_bandwidth),
   _JokerSpikeB(_bandwidth,_bandwidth),
   _g(_bandwidth,_bandwidth),
   _yL(_ncells,_bandwidth),
   _yR(_bandwidth,_bandwidth)
  {
    initAssembly();
    logCtor();
  }

  void solve()
  {

  }

  virtual ~SpikeMatrix()
  {
    logDtor();
  }


private:
  //initial assamely 
  void initAssembly()
  {
#ifdef FVM_PARALLEL
     _procID = MPI::COMM_WORLD.Get_rank();
     _nprocs = MPI::COMM_WORLD.Get_size();
#endif
     setMatrix();
     setLMtrx();
     setRMtrx();
     lu();
     setLSpikeMtrx();
     setRSpikeMtrx();
     exchangeSpikeMtrx();
     //setRSpikeMtrxFull();
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
     
     //_A.print(cout);
  } 
  //left matrix
  void setLMtrx()
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
     //_L.print(cout);
  }
  
  //lu 
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
    //_A.print(cout);
  }

  //right matrix
  void setRMtrx()
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
     /*_R.print(cout);*/
  }
  
  //left Spike Mtrx
  void setLSpikeMtrx()
  {
     //copy _L to LSpike
     _LL.partialCopyFrom(_L);
     //_LL.print(cout);
     //zeros yL
     _yL.zeros();
     const int b = _bandwidth;
     for ( int n = 0; n < b; n++ ){
        _yL(0,n) = _LL(0,n);
	for ( int i = 1; i < _ncells; i++ ){
            Diag yi = _LL(i,n);
	    for ( int j = max(0,i-b); j <= i-1; j++ ){
	        const int i2 = b+i-j;
	        yi -=  _A(i2,j) * _yL(j,n);
	    }
            _yL(i,n) = yi;
	}
    }
     //_yL.print(cout);
    //backward solve
    for ( int n = 0; n < b; n++ ){
       _LSpike(_ncells-1,n) = _yL(_ncells-1,n) / _A(b,_ncells-1);
       for ( int i = _ncells-2; i >= 0; i-- ){
          Diag soli = _yL(i,n);
	  for ( int j = i+1; j <= min(_ncells-1,i+b); j++ ){
	     const int i2 = b+i-j;
              soli -=  _A(i2,j)*_LSpike(j,n);
	  }
	  _LSpike(i,n) = soli / _A(b,i);
      }
   }
    //_LSpike.print(cout);
   _LSpike.partialCopyTo(_LSpikeT);
   _LSpikeT.print(cout);

  }
  //Right Spike Mtrx
  void setRSpikeMtrxFull()
  {
     //copy _R to RSpike
     _RR.partialCopyFrom(_R);
     //_RR.print(cout);
     //zeros yR
     _yR.zeros();
     const int b = _bandwidth;
     for ( int n = 0; n < b; n++ ){
        _yR(0,n) = _RR(0,n);
	for ( int i = 1; i < _ncells; i++ ){
            Diag yi = _RR(i,n);
	    for ( int j = max(0,i-b); j <= i-1; j++ ){
	        const int i2 = b+i-j;
	        yi -=  _A(i2,j) * _yR(j,n);
	    }
            _yR(i,n) = yi;
	}
    }
    //_yR.print(cout);
    //backward solve
    for ( int n = 0; n < b; n++ ){
       _RSpike(_ncells-1,n) = _yR(_ncells-1,n) / _A(b,_ncells-1);
       for ( int i = _ncells-2; i >= 0; i-- ){
          Diag soli = _yR(i,n);
	  for ( int j = i+1; j <= min(_ncells-1,i+b); j++ ){
	     const int i2 = b+i-j;
              soli -=  _A(i2,j)*_RSpike(j,n);
	  }
	  _RSpike(i,n) = soli / _A(b,i);
      }
   }
   //_RSpike.print(cout);
   _RSpike.partialCopyTo(_RSpikeB);
   _RSpikeB.print(cout);

  }


  //right Spike Mtrx
  void setRSpikeMtrx()
  {
     //forward solve 
     //zeros yR
     _yR.zeros();
     const int b = _bandwidth;
     for ( int n = 0; n < b; n++ ){
        _yR(0,n) = _R(0,n);
	for ( int i = 1; i < b; i++ ){
	    const int ii =  _ncells-b+i;
            Diag yi = _R(i,n);
	    for ( int j = 0; j < i; j++ ){
	        const int jj = _ncells-b + j;
	        const int i2 = b+ii-jj;
	        yi -=  _A(i2,jj) * _yR(j,n);
	    }
            _yR(i,n) = yi;
	}
    }
    //_yR.print(cout);
    //backward solve
    for ( int n = 0; n < b; n++ ){
       _RSpikeB(b-1,n) = _yR(b-1,n) / _A(b,_ncells-1);
       for ( int i = b-2; i >= 0; i-- ){
          Diag soli = _yR(i,n);
	  const int ii = _ncells-b+i;
	  for ( int j = b-1; j > i; j-- ){
	     const int jj = _ncells-b+j;
	     const int i2 = b+ii-jj;
              soli -=  _A(i2,jj)*_RSpikeB(j,n);
	  }
	  _RSpikeB(i,n) = soli / _A(b,ii);
      }
   }
   _RSpikeB.print(cout);
 }

 //exchanging LSpikeT (to rank-1) and RSpikeB(to rank+1), will be stored _JokerSpike Mtrx 
 void   exchangeSpikeMtrx()
 {
 #ifdef FVM_PARALLEL
   //send-recv single call since each process is involving send and recv
   MPI::Status status;  
   if ( _procID != _nprocs-1)
   MPI::COMM_WORLD.Sendrecv(_RSpikeB.getData()    , _RSpikeB.getDataSize()    , MPI::BYTE, _procID+1, 1199,
                            _JokerSpikeT.getData(), _JokerSpikeT.getDataSize(), MPI::BYTE, _procID+1, 2199, status );
   if ( _procID != 0 )
   MPI::COMM_WORLD.Sendrecv(_LSpikeT.getData()    , _LSpikeT.getDataSize()    , MPI::BYTE, _procID-1, 2199,
		   _JokerSpikeB.getData(), _JokerSpikeB.getDataSize(), MPI::BYTE, _procID-1, 1199, status );
   _JokerSpikeB.print(cout);
   _JokerSpikeT.print(cout);
#endif
				   

 }

  const CRConnectivity& _conn;
  const Array<Diag>& _diag;
  const Array<OffDiag>& _offDiag;
  const SpikeStorage& _spikeStorage;
  const	int _bandwidth;
  const int _ncells;
  int  _procID;
  int  _nprocs;
  Array2D<Diag>  _A;
  Array2D<Diag>  _LL;     //C matrix full Nxb 
  Array2D<Diag>  _L;      //C matrix only bxb
  Array2D<Diag>  _RR;     //B matrix full Nxb
  Array2D<Diag>  _R;      //B matrix bxb
  Array2D<Diag>  _LSpike;  //full left spike matrix
  Array2D<Diag>  _LSpikeT; //left-top spike matrix
  Array2D<Diag>  _RSpike;  //full right spike matrix
  Array2D<Diag>  _RSpikeB; //right-top spike matrix
  Array2D<Diag>  _JokerSpikeT; // bxb matrix come from rank+1(source) 
  Array2D<Diag>  _JokerSpikeB; // bxb matrix come from rank-1(source)  
  Array2D<Diag>  _g;       //g matrix LU g = f
  Array2D<Diag>  _yL; 
  Array2D<Diag>  _yR;       //joker vector for intermiediate steps, 
                           //LU x = f,  Ux = y first solve L y = f and then solve U x = y 
};


#endif
