// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

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
#include "Array.h"

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
   _reducedA1(_bandwidth,_bandwidth),
   _reducedA2(_bandwidth,_bandwidth),
   _reducedRHS1(_bandwidth),
   _reducedRHS2(_bandwidth),
   _JokerZ1(_bandwidth),
   _JokerZ2(_bandwidth),
   _RHS(_ncells),
   _g(_ncells),
   _gB(_bandwidth),
   _gT(_bandwidth),
   _JokergB(_bandwidth),
   _JokergT(_bandwidth),
   _yL(_ncells,_bandwidth),
   _yR(_bandwidth,_bandwidth),
   _y(_ncells),
   _pp1(_bandwidth),
   _pp2(_bandwidth)
  {
    initAssembly();
    logCtor();
  }

  void solve( const XArray& f, XArray& x )
  {
      //negate rhs (f) for x
      luSolver(f, _g, true);
      //set gB, gT from 
      setgBgT();
      //communicate gB and gT to neighbourhood processors
      exchange_gTgB();
      //setting rhs for reduced system
      setReducedRHS();
      //solve reduced system
      solveReducedSystem();
      //exchange reduced system solutions
      exchange_reducedSol(); 
      //setting final system rhs 
      setRHS(f, true);
      //get solution
      luSolver( _RHS, x,  false);
      
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
     setReducedMtrx();
     //setRSpikeMtrxFull();
     if (_procID != _nprocs-1)
        denseMtrxLU (_reducedA1, _pp1);
     if ( _procID != 0 )
        denseMtrxLU (_reducedA2, _pp2);
  }
 
  void setMatrix()
  {
     //forming A
     //const int  nr = 2 * _bandwidth + 1;
     const int  nc = _ncells;
     //diagonal filling
     for ( int i = 0; i < nc; i++ )
         _A(_bandwidth,i) = _diag[i];
     const Array<int>& row = _conn.getRow();
     const Array<int>& col = _conn.getCol();
     for ( int i = 0; i < nc; i++ ){
         for (int n = row[i]; n < row[i+1]; n++ ){
	     int j =  col[n];
	     //check if it is in bandwidth and inner coefficient
	     if ( abs(j-i) <= _bandwidth && j < nc ){
	       _A(_bandwidth-(j-i),j) = _offDiag[n];
	     }
         }
     }
     
     /*_A.print(cout);*/
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
     /*_L.print(cout);*/
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
     /*_A.print(cout);*/
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
     /*_LL.print(cout);*/
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
     /*_yL.print(cout);*/
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
    /*_LSpike.print(cout);*/
   _LSpike.partialCopyTo(_LSpikeT);
   /*_LSpikeT.print(cout);*/

  }
  //Right Spike Mtrx
  void setRSpikeMtrxFull()
  {
     //copy _R to RSpike
     _RR.partialCopyFrom(_R);
     /*_RR.print(cout);*/
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
     /*_yR.print(cout);*/
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
    /*_RSpike.print(cout);*/
   _RSpike.partialCopyTo(_RSpikeB);
   /*_RSpikeB.print(cout);*/

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
     /*_yR.print(cout);*/
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
    /*_RSpikeB.print(cout);*/
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
   /*_JokerSpikeB.print(cout);*/
   /*_JokerSpikeT.print(cout);*/
#endif
				   

 }

//settting reduced matrix 
void  setReducedMtrx()
{
    //system above the processor line 
    if ( _procID != _nprocs-1 ){
	_reducedA1.setIdentity();
	for ( int i = 0; i < _bandwidth; i++ ){
	   for ( int j = 0; j < _bandwidth; j++ ){
	      for ( int k = 0; k < _bandwidth; k++ ){
	         _reducedA1(i,j) -=  _RSpikeB(i,k) * _JokerSpikeT(k,j);
	      }
	   }
	}
   }
    /*_reducedA1.print(cout);*/

    //system above the processor line 
    if ( _procID != 0 ){
	_reducedA2.setIdentity();
	for ( int i = 0; i < _bandwidth; i++ ){
	   for ( int j = 0; j < _bandwidth; j++ ){
	      for ( int k = 0; k < _bandwidth; k++ ){
	         _reducedA2(i,j) -=  _LSpikeT(i,k) * _JokerSpikeB(k,j);
	      }
	   }
	}
   }
    /*_reducedA2.print(cout);*/
}

//getting LU decomposigion for dense matrix, the LU decomposition is stored in 
void denseMtrxLU ( Array2D<Diag>&  A,  Array<int>& pp)
{
	/*A.print(cout);*/
    //fill permutation
    for( int i = 0; i < pp.getLength(); i++ )
       pp[i] = i;

    const int n = A.getRow();
    for ( int i = 0; i < n-1; i++ ){
       //doing paritial (column) pivoting
       Diag pivot = A(i,i);
       int jj = i;
       for ( int j = i+1; j < n; j++){
           if ( NumTypeTraits<Diag>::doubleMeasure( A(j,i) ) > NumTypeTraits<Diag>::doubleMeasure( pivot ) ){
	       jj = j;
	       pivot = A(j,i);
	   }
       }
       //explicitly swap A(i,i:n-1) and A(jj,i:n-1)
       if ( jj > i ){
         for ( int j = i; j < n; j++){ 
            const Diag tmp = A(i,j);
            A(i,j) = A(jj,j);
	    A(jj,j) = tmp;
	 }
	 pp[i] = jj;
       }
       ///////////
       pivot = A(i,i);
       for ( int j = i+1; j < n; j++ ){
           const Diag m = A(j,i) / pivot;
	   A(j,i) = m;
	   for ( int k = i+1; k < n; k++ ){
	       A(j,k) -= m * A(i,k);
	   }
       }
    }
    /*pp.print(cout);*/
    /*A.print(cout);*/
      

}


// generalized Lu solve Lu x = f as Vector
  void luSolver(const Array<X>& f, Array<X>& x, bool negate_rhs=false)
  {
      x.zero();
     //zeros y
     _y.zero();
      const int b = _bandwidth;
     if ( negate_rhs ){
        _y[0] = -f[0];
        for ( int i = 1; i < _ncells; i++ ){
           X yi = -f[i];
           for ( int j = max(0,i-b); j <= i-1; j++ ){
              const int i2 = b+i-j;
              yi -=  _A(i2,j) * _y[j];
           }
           _y[i] = yi;
        }
    } else {
        _y[0] = f[0];
        for ( int i = 1; i < _ncells; i++ ){
           X yi = f[i];
           for ( int j = max(0,i-b); j <= i-1; j++ ){
              const int i2 = b+i-j;
              yi -=  _A(i2,j) * _y[j];
           }
           _y[i] = yi;
        }
    }
     /*_y.print(cout);*/
    //backward solve
    x.zero();
    x[_ncells-1] = _y[_ncells-1] / _A(b,_ncells-1);
    for ( int i = _ncells-2; i >= 0; i-- ){
        X soli = _y[i];
        for ( int j = i+1; j <= min(_ncells-1,i+b); j++ ){
            const int i2 = b+i-j;
            soli -=  _A(i2,j)*x[j];
        }
	x[i] = soli / _A(b,i);
    }
    /*f.print(cout);*/
    /*cout << endl;*/
    /*_y.print(cout);*/
    /*cout << endl;*/
    /*x.print(cout);*/
    /*cout << endl;*/
  }

  //setting gB and gT from g
  void setgBgT()
  {
      //top part of g
      _gT.zero();
      for ( int i = 0; i < _bandwidth; i++ )
          _gT[i] = _g[i];
      //bottom part of g
      _gB.zero();
      int indx = 0;
      for ( int i = _ncells-_bandwidth; i < _ncells; i++ )
          _gB[indx++] = _g[i];
      /*_g.print(cout);*/
      /*_gT.print(cout);*/
      /*_gB.print(cout);*/

  }
  //exhange gT and gB to temporary arrays
  void exchange_gTgB()
  {
     _JokergB.zero();
     _JokergT.zero();
     #ifdef FVM_PARALLEL
     //send-recv single call since each process is involving send and recv
     MPI::Status status;  
     if (_procID != _nprocs-1)
     MPI::COMM_WORLD.Sendrecv(_gB.getData()     , _gB.getDataSize()     , MPI::BYTE, _procID+1, 3199,
                              _JokergT.getData(), _JokergT.getDataSize(), MPI::BYTE, _procID+1, 4199, status );
     if ( _procID != 0 )
     MPI::COMM_WORLD.Sendrecv(_gT.getData()    , _gT.getDataSize()      , MPI::BYTE, _procID-1, 4199,
                 	      _JokergB.getData(), _JokergB.getDataSize(), MPI::BYTE, _procID-1, 3199, status );
     #endif
     /*_JokergB.print(cout);*/
     /*_JokergT.print(cout);*/
	
  }
  
  //setting rhs for reduced system
  void setReducedRHS()
  {
      //rhs1
      _reducedRHS1.zero();
      for ( int i = 0; i < _bandwidth; i++ ){
          X dot_product = NumTypeTraits<X>::getZero();
	  for ( int j = 0; j < _bandwidth; j++ ){
	       dot_product +=   _RSpikeB(i,j) * _JokergT[j];
          }
	  _reducedRHS1[i] = _gB[i] - dot_product;
      }
      /*_reducedRHS1.print(cout);*/
      //rhs2
      _reducedRHS2.zero();
      for ( int i = 0; i < _bandwidth; i++ ){
          X dot_product = NumTypeTraits<X>::getZero();
	  for ( int j = 0; j < _bandwidth; j++ ){
	       dot_product += _LSpikeT(i,j) * _JokergB[j];
          }
	  _reducedRHS2[i] = _gT[i] - dot_product;
      }
      /*_reducedRHS2.print(cout);*/
  }
  //solving reduced system
  void solveReducedSystem()
  {
     if (_procID != _nprocs-1)
        denseLUsolve(_reducedA1, _pp1, _reducedRHS1); //solution z1 is stored in reducedRHS1
     if ( _procID != 0 )
        denseLUsolve(_reducedA2, _pp2, _reducedRHS2);//solution z2 is stored in reducedRHS2
  }
  //dens LU solver
  void denseLUsolve( const Array2D<Diag>& A, const Array<int>& pp, Array<X>& rhs )
  {
	  /*rhs.print(cout);*/
     const int n = A.getRow();
     //forward solve
     for ( int i = 0; i < n-1; i++ ){
        //swap components of y: y(i) <---> y(pp(i))
        if ( pp[i] > i ){
	    X tmp = rhs[i];
	    rhs[i] = rhs[pp[i]];
	    rhs[pp[i]]=tmp;
        }
	for ( int j = i+1; j < n; j++ ){
	   rhs[j] -= A(j,i) * rhs[i];
	}
     }

     //back solve (later =/ might be useful to define)
     if ( n >= 1 ) //prevent accessing arrays if n=0 size
        rhs[n-1] = rhs[n-1] /  A(n-1,n-1);

     for ( int i = n-2; i >=0; i-- ){
        for ( int j = i+1; j < n; j++ ){
	   rhs[i] -= A(i,j) * rhs[j];
        }
	rhs[i] = rhs[i] /  A(i,i);
     }
     /*cout << endl;*/
     /*rhs.print(cout);*/
     /*cout << endl;*/
  }

//exchange reduced system solution between neighbourhodds
  void exchange_reducedSol()
  {
     _JokerZ1.zero();
     _JokerZ2.zero();
     #ifdef FVM_PARALLEL
     //send-recv single call since each process is involving send and recv
     MPI::Status status;  
     if (_procID != _nprocs-1)
     MPI::COMM_WORLD.Sendrecv(_reducedRHS1.getData(), _reducedRHS1.getDataSize(), MPI::BYTE, _procID+1, 3199,
                              _JokerZ2.getData(), _JokerZ2.getDataSize(), MPI::BYTE, _procID+1, 4199, status );
     if ( _procID != 0 )
     MPI::COMM_WORLD.Sendrecv(_reducedRHS2.getData(), _reducedRHS2.getDataSize(), MPI::BYTE, _procID-1, 4199,
                 	      _JokerZ1.getData(), _JokerZ1.getDataSize(), MPI::BYTE, _procID-1, 3199, status );
     #endif
     /*cout << endl;*/
     /*_JokerZ1.print(cout);*/
     /*cout << endl;*/
     /*_JokerZ2.print(cout);*/
     /*cout << endl;*/
	
  }
//setting rhs for  final 
  void setRHS(const Array<X>& f,  bool negate_rhs=false)
  {
      _RHS.zero();
      if ( negate_rhs == true ){
         for ( int i = 0; i < _ncells; i++ )
	     _RHS[i] = -f[i];
      } else {
         for ( int i = 0; i < _ncells; i++ )
	     _RHS[i] = f[i];
      }
      //top part
      for ( int i = 0; i < _bandwidth; i++ ){
         X dot_product = NumTypeTraits<X>::getZero();
         for ( int j = 0; j < _bandwidth; j++ ){
		 dot_product +=   _L(i,j) * _JokerZ1[j];
         }
         _RHS[i] -=  dot_product;
      }
      //bottom part
      for ( int i = 0; i < _bandwidth; i++ ){
         X dot_product = NumTypeTraits<X>::getZero();
         const int indx = _ncells-_bandwidth+i;
         for ( int j = 0; j < _bandwidth; j++ ){
		 dot_product +=  _R(i,j) * _JokerZ2[j];
         }
         _RHS[indx] -=  dot_product;
      }
      /*_RHS.print(cout);*/
      /*cout << endl;*/
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
  Array2D<Diag>  _reducedA1;    // bxb matrix (I-V * W )z1 = g1 - V g2 (system above processor boundary line)
  Array2D<Diag>  _reducedA2;    // bxb matrix (I-W * V )z2 = g2 - V g1 (system below processor boundayr line)
  Array<X>       _reducedRHS1;  // g1 - V g2
  Array<X>       _reducedRHS2;  // g2 - W g1
  Array<X>       _JokerZ1;      //store z1 solution from other process(rank-1)
  Array<X>       _JokerZ2;      //store z2 solution from other process(rank+1)
  Array<X>       _RHS;          //final step f1 - C2*z1 - B2*z4
  Array<X>       _g;       //g matrix LU g = f
  Array<X>       _gB;    //bottom part of g, gB = g(ncells-b,ncells-1)
  Array<X>       _gT;    //top part of g, gT = g(0:b-1)
  Array<X>       _JokergB;  //from top process 
  Array<X>       _JokergT;  //from bottom process
  Array2D<Diag>  _yL; 
  Array2D<Diag>  _yR;      //joker vector for intermiediate steps, 
                           //LU x = f,  Ux = y first solve L y = f and then solve U x = y 
  Array<X>        _y;      //joker vector			   
  Array<int>      _pp1;
  Array<int>      _pp2;   //permutation vectors
};


#endif
