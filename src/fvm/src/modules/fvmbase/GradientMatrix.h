// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _GRADIENTMATRIX_H_
#define _GRADIENTMATRIX_H_

#include <set>

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

#include "Mesh.h"
#include "CRConnectivity.h"
#include "Array.h"
#include "StorageSite.h"
#include "Gradient.h"

class GradientMatrixBase
{
 public:
  typedef pair<const StorageSite*, const StorageSite*> EntryIndex;
  typedef map<EntryIndex, shared_ptr<ArrayBase> > GhostArrayMap;
  GradientMatrixBase() {}
  virtual void syncLocal() {};
  virtual ~GradientMatrixBase() {}
};

template<class T_Scalar>
class GradientMatrix : public GradientMatrixBase
{
public:
  typedef Vector<T_Scalar,3> Coord;
  
  GradientMatrix(const Mesh& mesh) :
    _mesh(mesh),
#ifdef FVM_PARALLEL
    _conn(mesh.getCellCellsGhostExt()),
    _row(_conn.getRow() ),
    _col(_conn.getCol() ),
#else
    _conn(mesh.getCellCells()),
    _row(_conn.getRow()),
    _col(_conn.getCol()),
#endif
    _coeffs(_col.getLength()),
    _pairWiseAssemblers()
  {
  }

  virtual ~GradientMatrix()
  {}
  
  template<class X>
  shared_ptr<Array<Gradient<X> > >
  getGradient(const Array<X>& x) const
  {
    typedef Gradient<X> GradType;
    typedef Array<GradType> GradArray;

    shared_ptr<GradArray> gradXPtr(new GradArray(x.getLength()));
    GradArray& gradX = *gradXPtr;
    
    const int nRows = getConnectivity().getRowSite().getSelfCount();
    for(int nr=0; nr<nRows; nr++)
    {
        gradX[nr].zero();
        for (int nb = _row[nr]; nb<_row[nr+1]; nb++)
        {
            const int j = _col[nb];
            gradX[nr].accumulate(_coeffs[nb],x[j]-x[nr]);
        }
    }
    return gradXPtr;
  }
  
  template<class X>
  void
  computeGradient(Gradient<X>& g, const Array<X>& x, int i) const
  {
    g.zero();
   
    // for boundaries use the adjacent cell
    if (_row[i+1] - _row[i]  == 1)
      i = _col[_row[i]];
    
    for (int nb = _row[i]; nb<_row[i+1]; nb++)
    {
        const int j = _col[nb];
        g.accumulate(_coeffs[nb],x[j]-x[i]);
    }
  }

  template<class X>
    X
    computeR(const Gradient<X>& g, const Array<X>& x, const Coord dist, int i, int j) const
    {//Darwish and Moukalled, Int. J. H. M. T., 46 (2003) 599-611

      X den=x[j]-x[i];
      X num=2.*(g*dist);
      if (den!=0.)
	return num/den-1.;
      return 0;
    }

  template<class X>
    void
    computeFaceGradient(Gradient<X>& g, const Array<X>& x, int i) const
    {
      g.zero();
   
    // for boundaries use the adjacent cell
    if (_row[i+1] - _row[i]  == 1)
      i = _col[_row[i]];
    
    for (int nb = _row[i]; nb<_row[i+1]; nb++)
    {
        const int j = _col[nb];
        g.accumulate(_coeffs[nb],x[j]);
    }
  }
  
  const CRConnectivity& getConnectivity() const {return _conn;}
  
  Array<Coord>& getCoeffs() {return _coeffs;}
  const Array<Coord>& getCoeffs() const {return _coeffs;}

  Coord& getCoeff(const int i,  const int j)
  {
    for (int nnb = _row[i]; nnb<_row[i+1]; nnb++)
    {
        if (_col[nnb] == j)
          return _coeffs[nnb];
    }
   
    throw CException("invalid indices");
  }
  

  const Coord& getCoeff(const int i,  const int j) const
  {
    for (int nnb = _row[i]; nnb<_row[i+1]; nnb++)
    {
        if (_col[nnb] == j)
          return _coeffs[nnb];
    }
    throw CException("invalid indices");
  }




  class PairWiseAssembler
  {
  public:
    PairWiseAssembler(Array<Coord>& coeffs,
                      const Array<Vector<int,2> >& pairToCol) :
      _coeffs(coeffs),
      _pairToCol(pairToCol)
    {}

    Coord& getCoeff01(const int np)
    {
      return _coeffs[_pairToCol[np][0]];
    }
    
    Coord& getCoeff10(const int np)
    {
      return _coeffs[_pairToCol[np][1]];
    }
  private:
  

    Array<Coord>& _coeffs;
    const Array<Vector<int,2> >& _pairToCol;
  };
  
  PairWiseAssembler& getPairWiseAssembler(const CRConnectivity& pairs)
  {
    if (_pairWiseAssemblers.find(&pairs) == _pairWiseAssemblers.end())
    {
        _pairWiseAssemblers[&pairs] =
          new PairWiseAssembler(_coeffs,
                                getConnectivity().getPairToColMapping(pairs));
    }
    return *_pairWiseAssemblers[&pairs];
  }
  

    //fill scatterArray (both mesh and partition) and only gatherArray for mesh
    void    createScatterGatherValuesBuffer()
    {
#ifdef FVM_PARALLEL
       //SENDING allocation and filling
       const StorageSite& site     = _mesh.getCells();
       const CRConnectivity& cellCells = _mesh.getCellCells(); 
       const StorageSite::ScatterMap& scatterMap = site.getScatterMap();
       foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
           const StorageSite&  oSite = *mpos.first;
           const Array<int>& scatterArray = *mpos.second;
           //loop over surround indices and itself for sizing
           EntryIndex e(&site,&oSite);
           int sendSize = 0;
           for ( int i = 0; i < scatterArray.getLength(); i++ ){
              sendSize += cellCells.getCount( scatterArray[i] );
           }
           //allocate array
           _sendValues[e]  = shared_ptr< Array<Coord> > ( new Array<Coord> (sendSize) );
           //fill send array
           Array<Coord>& valueArray = dynamic_cast< Array<Coord>& > ( *_sendValues[e]  );
           int indx = 0;
           for( int i = 0; i < scatterArray.getLength(); i++ ){
              const int ii = scatterArray[i];
              for ( int j = 0; j < cellCells.getCount(ii); j++ ){
                 const int jj = cellCells(ii,j);
                 valueArray[indx] = getCoeff(ii,jj);
                 indx++;
              }
           }
       }
       //RECIEVING allocation (filling will be done by MPI Communication)
       const StorageSite::GatherMap& gatherMap = site.getGatherMap();
       foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap){
          const StorageSite&  oSite = *mpos.first;
          EntryIndex e(&oSite,&site);
          //get recvSize
          const Array<int>& recvCounts   =  dynamic_cast< const Array<int>& > (_mesh.getRecvCounts(e));
          int recvSize = 0;
          for ( int i = 0; i < recvCounts.getLength(); i++ ){
             recvSize += recvCounts[i];
           }
          //allocate array
          _recvValues [e] = shared_ptr< Array<Coord> > ( new Array<Coord> (recvSize) );
       }
#endif
    }
    
       //fill scatterArray (both mesh and partition) and only gatherArray for mesh
    void    recvScatterGatherValuesBufferLocal()
    {
#ifdef FVM_PARALLEL
       //RECIEVING allocation (filling will be done by MPI Communication)
       const StorageSite& site     = _mesh.getCells();
       const StorageSite::GatherMap& gatherMap = site.getGatherMap();
       foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap){
          const StorageSite&  oSite = *mpos.first;
          EntryIndex e(&oSite,&site);
          //get recvSize
          //mesh interface can be done know
          if ( oSite.getGatherProcID() == - 1) {
             *_recvValues[e] = *_sendValues [e];
          } 
       }
#endif
    }

    //sending values
    void syncValues()
    {
#ifdef FVM_PARALLEL
          //SENDING
          const int  request_size = get_request_size();
          MPI::Request   request_send[ request_size ];
          MPI::Request   request_recv[ request_size ];
          int indxSend = 0;
          int indxRecv = 0;
          const StorageSite&    site      = _mesh.getCells();
          const StorageSite::ScatterMap& scatterMap = site.getScatterMap();
          foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
             const StorageSite&  oSite = *mpos.first;
             EntryIndex e(&site,&oSite);
             //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
             ArrayBase& sendArray = *_sendValues[e];
             //loop over surround indices and itself
             int to_where  = oSite.getGatherProcID();
             if ( to_where != -1 ){
                int mpi_tag = oSite.getTag();
                request_send[indxSend++] =  
                    MPI::COMM_WORLD.Isend( sendArray.getData(), sendArray.getDataSize(), MPI::BYTE, to_where, mpi_tag );
             }
          }
          //RECIEVING
          //getting values from other meshes to fill g
          const StorageSite::GatherMap& gatherMap = site.getGatherMap();
          foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap){
             const StorageSite&  oSite = *mpos.first;
             //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
             EntryIndex e(&oSite,&site);
             ArrayBase& recvArray = *_recvValues[e];
             int from_where       = oSite.getGatherProcID();
             if ( from_where != -1 ){
                int mpi_tag = oSite.getTag();
                request_recv[indxRecv++] =  
                      MPI::COMM_WORLD.Irecv( recvArray.getData(), recvArray.getDataSize(), MPI::BYTE, from_where, mpi_tag );
             }
          }

          int count  = get_request_size();
          MPI::Request::Waitall( count, request_recv );
          MPI::Request::Waitall( count, request_send );

//           const StorageSite&    site      = _mesh.getCells();
//           const StorageSite::GatherMap& gatherMap = site.getGatherMap();

	  //replacing values
          const map<int,int>&   globalToLocal = _mesh.getGlobalToLocal();
          foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap){
             const StorageSite&  oSite = *mpos.first;
             const Array<int>& gatherArray = dynamic_cast< const Array<int>& > (*mpos.second);
             //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
             int from_where       = oSite.getGatherProcID();
             if ( from_where != -1 ){
                EntryIndex e(&oSite,&site);
                const Array<int>  & recv_counts  = dynamic_cast< const Array<int>&   > (_mesh.getRecvCounts(e) );
                const Array<int>  & recv_indices = dynamic_cast< const Array<int>&   > (_mesh.getRecvIndices(e));
                const Array<Coord>& recv_values  = dynamic_cast< const Array<Coord>& > (*_recvValues [e]);

                //loop over gatherArray
                int indx = 0;
                for ( int i = 0; i < gatherArray.getLength(); i++ ){
                   const int ii  = gatherArray[i];
                   const int nnb = recv_counts[i]; //give getCount() 
                   for ( int nb = 0; nb < nnb; nb++ ){
                      const int jj = globalToLocal.find( recv_indices[indx] )->second;
                       Coord& matrix_entry = getCoeff(ii,jj);
                       matrix_entry = recv_values[indx];
                       indx++;
                   }
                }
             }
          }

#endif
    }





private:

  int get_request_size()
  {
      int indx =  0;
      const StorageSite& site = _mesh.getCells();
      const StorageSite::ScatterMap& scatterMap = site.getScatterMap();
      foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
         const StorageSite&  oSite = *mpos.first;
         //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
         if ( oSite.getGatherProcID() != -1 )
            indx++;
      }
      return indx;
  }
 
 
  virtual void printRow(const int i) const
  {
    
    for (int nb = _row[i]; nb<_row[i+1]; nb++)
    {
        const int j = _col[nb];
        cout << "coeff (" <<  i << "," << j << ") = " << _coeffs[nb] << endl;
    }
  }

  const Mesh&   _mesh;
  const CRConnectivity& _conn;
  const Array<int>& _row;
  const Array<int>& _col;

  
  Array<Coord> _coeffs;
  mutable map<const CRConnectivity*,PairWiseAssembler*> _pairWiseAssemblers;

#ifdef FVM_PARALLEL
  GhostArrayMap   _sendValues;
  GhostArrayMap   _recvValues;
#endif
};


#endif
