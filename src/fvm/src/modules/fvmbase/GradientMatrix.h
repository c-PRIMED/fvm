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

  typedef pair<const StorageSite*, const StorageSite*> EntryIndex;
  typedef map<EntryIndex, shared_ptr<ArrayBase> > GhostArrayMap;

public:
  GradientMatrixBase() {}
  virtual ~GradientMatrixBase() {}
};

template<class T_Scalar>
class GradientMatrix : public GradientMatrixBase
{
public:
  typedef Vector<T_Scalar,3> Coord;
  
  GradientMatrix(const Mesh& mesh) :
    _mesh(mesh),
    _row(0),
    _col(0),
    _coeffs(_col.getLength()),
    _pairWiseAssemblers()
  {
     init();
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
    
    const int nRows = _conn->getRowSite().getSelfCount();
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
  

  const CRConnectivity& getConnectivity() const {return *_conn;}

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


  void sync()
  {
       createScatterGatherValuesBuffer();
       syncValues();
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
                                _conn->getPairToColMapping(pairs));
    }
    return *_pairWiseAssemblers[&pairs];
  }
  
private:


   void init()
   {
       //create recvCounts buffer
       createScatterGatherCountsBuffer();
       //syn counts
       syncCounts();
       //create recvindices buffer
       createScatterGatherIndicesBuffer();
      //sync indices 
       syncIndices();
      //create CRConnectivity (cellCells) new one take account into ghostcell connecitivty
       CRConn(); 
       //resize class members arrays, _row, _col, _colls
       fillClassMemberArrays();
   }

    void CRConn()
    {
       createRowColSiteCRConn();
       countCRConn();
       addCRConn();
       //CRConnectivityPrint(*_conn, 2, "cellCells");
 
    }

    void createRowColSiteCRConn()
    { 
        //counting interface counts
        //counting interfaces
        set<int> interfaceCells;
        const Array<int>& localToGlobal = _mesh.getLocalToGlobal();
        const map<int,int>& globalToLocal = _mesh.getGlobalToLocal();
        const StorageSite& site = _mesh.getCells();
        const int originalCount = site.getCount();
        const StorageSite::GatherMap& gatherMap = site.getGatherMap();
        int countLevel0 = 0;
        foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap){
           const StorageSite&  oSite = *mpos.first;
           const Array<int>& gatherArray = dynamic_cast< const Array<int>& > (*mpos.second);
           countLevel0 += gatherArray.getLength();
           EntryIndex e(&oSite,&site);
           const Array<int>  & recv_indices = dynamic_cast< const Array<int>& > (*_recvIndices[e]);
           const Array<int>  & recv_counts  = dynamic_cast< const Array<int>& > (*_recvCounts [e]);
           //loop over gatherArray
           int indx = 0;
           for ( int i = 0; i < gatherArray.getLength(); i++ ){
              const int nnb = recv_counts[i]; //give getCount() 
              for ( int nb = 0; nb < nnb; nb++ ){
                 const int localID = globalToLocal.find( recv_indices[indx] )->second;
                 if ( localID >= originalCount ){
                    interfaceCells.insert( recv_indices[indx] );
                 }
                 indx++;
              }
           }
        }
        const int selfCount   = site.getSelfCount();
        const int countLevel1 = int(interfaceCells.size());
        //ghost cells = sum of boundary and interfaces
        const int nghost = getNumBounCells() +  countLevel0 + countLevel1;
        _newSite = shared_ptr<StorageSite> ( new StorageSite(selfCount, nghost) );
         //constructing new CRConnecitvity;
        _conn = shared_ptr< CRConnectivity> ( new CRConnectivity( *_newSite, *_newSite) );

    }
 

    void countCRConn()
    {
       CRConnectivity& conn = *_conn;
       conn.initCount();
       const StorageSite& site = _mesh.getCells();
       const CRConnectivity& cellCells = _mesh.getCellCells();
       int ncount = site.getSelfCount() + getNumBounCells();
       //loop over old connectivity (inner + boundary) 
       for ( int i = 0; i < ncount; i++ ){
          conn.addCount(i, cellCells.getCount(i) );
       }
       // now interfaces
       const StorageSite::GatherMap& gatherMap = site.getGatherMap();
       foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap){
          const StorageSite&  oSite = *mpos.first;
          const Array<int>& gatherArray = dynamic_cast< const Array<int>& > (*mpos.second);
          EntryIndex e(&oSite,&site);
          const Array<int>& recv_counts = dynamic_cast< const Array<int>& > (*_recvCounts [e]);
          //loop over gatherArray
          for ( int i = 0; i < gatherArray.getLength(); i++ ){
             conn.addCount(gatherArray[i], recv_counts[i] );
          }
       }
       //finishCount
       conn.finishCount();
    }

    void addCRConn()
    {
       CRConnectivity& conn = *_conn;
       const StorageSite& site = _mesh.getCells();
       const CRConnectivity& cellCells = _mesh.getCellCells();
       int ncount = site.getSelfCount() + getNumBounCells();
       //loop over olde connectivity (inner + boundary) 
       for ( int i = 0; i < ncount; i++ ){
          for ( int j = 0; j < cellCells.getCount(i); j++ ){
             conn.add(i, cellCells(i,j));
          }
       }
       // now interfaces
       const map<int,int>& globalToLocal = _mesh.getGlobalToLocal();
       const StorageSite::GatherMap& gatherMap = site.getGatherMap();
       foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap){
          const StorageSite&  oSite = *mpos.first;
          const Array<int>& gatherArray = dynamic_cast< const Array<int>& > (*mpos.second);
          EntryIndex e(&oSite,&site);
          const Array<int>& recv_counts  = dynamic_cast< const Array<int>& > (*_recvCounts [e]);
          const Array<int>& recv_indices = dynamic_cast< const Array<int>& > (*_recvIndices[e]);
          //loop over gatherArray
          int indx = 0;
          for ( int i = 0; i < gatherArray.getLength(); i++ ){
             const int ncount = recv_counts[i];
             for ( int j = 0; j < ncount; j++ ){
                const int addCell = globalToLocal.find( recv_indices[indx] )->second;
                conn.add(gatherArray[i], addCell);
                indx++;
             }
          }
       }
       //finish add
       conn.finishAdd();

    }

    int getNumBounCells()
    {
       //boundary information has been stored
       const FaceGroupList&  boundaryFaceGroups = _mesh.getBoundaryFaceGroups();
       int nBounElm = 0;
       for ( int bounID = 0; bounID < int(boundaryFaceGroups.size()); bounID++){
          nBounElm += boundaryFaceGroups.at(bounID)->site.getCount();
       }
       return nBounElm; 
    }


    //fill countArray (both mesh and partition) and only gatherArray for mesh
    void    createScatterGatherCountsBuffer()
    {
       //SENDING allocation and filling
       const StorageSite& site     = _mesh.getCells();
       const CRConnectivity& cellCells = _mesh.getCellCells();
       const StorageSite::ScatterMap& scatterMap = site.getScatterMap();
       foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
           const StorageSite&  oSite = *mpos.first;
           //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
           const Array<int>& scatterArray = *mpos.second;
           cout << endl;
           //key site
           EntryIndex e(&site,&oSite);
           //allocate array
           _sendCounts[e] = shared_ptr< Array<int>   > ( new Array<int> (scatterArray.getLength()) ); 
           //fill send array
           Array<int>& sendArray = dynamic_cast< Array<int>& > ( *_sendCounts[e] );
           for( int i = 0; i < scatterArray.getLength(); i++ ){
              const int cellID = scatterArray[i];
              sendArray[i] = cellCells.getCount(cellID);
           }
       }
       //RECIEVING allocation (filling will be done by MPI Communication)
       const StorageSite::GatherMap& gatherMap = site.getGatherMap();
       foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap){
          const StorageSite&  oSite = *mpos.first;
          const Array<int>& gatherArray = *mpos.second;
          //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
          EntryIndex e(&oSite,&site);
          //allocate array
          _recvCounts[e] = shared_ptr< Array<int> > ( new Array<int> (gatherArray.getLength()) ); 
          //mesh interface can be done know
          if ( oSite.getGatherProcID() == - 1) {
             *_recvCounts[e] = *_sendCounts[e];
          } 
       }
    }

 void syncCounts()
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
           ArrayBase& sendArray = *_sendCounts[e];

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
          ArrayBase& recvArray = *_recvCounts[e];
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
#endif

    }



    //fill scatterArray (both mesh and partition) and only gatherArray for mesh
    void    createScatterGatherIndicesBuffer()
    {
       //SENDING allocation and filling
       const StorageSite& site     = _mesh.getCells();
       const CRConnectivity& cellCells = _mesh.getCellCells();
       const Array<int>&   localToGlobal = _mesh.getLocalToGlobal();
       const StorageSite::ScatterMap& scatterMap = site.getScatterMap();
       foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
           const StorageSite&  oSite = *mpos.first;
           const Array<int>& scatterArray = *mpos.second;
           //loop over surround indices and itself for sizing
           EntryIndex e(&site,&oSite);
           //allocate array
           int sendSize = 0;
           for ( int i = 0; i < scatterArray.getLength(); i++ ){
              sendSize += cellCells.getCount( scatterArray[i] );
           }
           _sendIndices[e] = shared_ptr< Array<int>   > ( new Array<int> (sendSize) ); 
           //fill send array
           Array<int>& sendArray = dynamic_cast< Array<int>&   > ( *_sendIndices[e] );
           int indx = 0;
           for( int i = 0; i < scatterArray.getLength(); i++ ){
              const int cellID = scatterArray[i];
              for ( int j = 0; j < cellCells.getCount(cellID); j++ ){
                 sendArray[indx] = localToGlobal[ cellCells(cellID,j) ];
                 indx++;
              }
           }
       }
       //RECIEVING allocation (filling will be done by MPI Communication)
       const StorageSite::GatherMap& gatherMap = site.getGatherMap();
       foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap){
          const StorageSite&  oSite = *mpos.first;
          EntryIndex e(&oSite,&site);
          const Array<int>& recvCounts   =  dynamic_cast< const Array<int>& > (*_recvCounts[e]);
          int recvSize = 0;
          for ( int i = 0; i < recvCounts.getLength(); i++ ){
             recvSize += recvCounts[i];
          }
          //allocate array
          _recvIndices[e] = shared_ptr< Array<int> > ( new Array<int>     (recvSize) ); 
          //mesh interface can be done know
          if ( oSite.getGatherProcID() == - 1) {
             *_recvIndices[e] = *_sendIndices[e];
          } 
       }

    }


  void syncIndices()
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
           ArrayBase& sendArray = *_sendIndices[e];
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
          EntryIndex e(&oSite,&site);
          ArrayBase& recvArray = *_recvIndices[e];
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
#endif

    }


    int  get_request_size()
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

    //fill scatterArray (both mesh and partition) and only gatherArray for mesh
    void    createScatterGatherValuesBuffer()
    {
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
          const Array<int>& recvCounts   =  dynamic_cast< const Array<int>& > (*_recvCounts[e]);
          int recvSize = 0;
          for ( int i = 0; i < recvCounts.getLength(); i++ ){
             recvSize += recvCounts[i];
           }
          //allocate array
          _recvValues [e] = shared_ptr< Array<Coord> > ( new Array<Coord> (recvSize) );
          //mesh interface can be done know
          if ( oSite.getGatherProcID() == - 1) {
             *_recvValues[e] = *_sendValues [e];
          } 
       }

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
#endif


#ifndef FVM_PARALLEL
          const StorageSite&    site      = _mesh.getCells();
          const StorageSite::GatherMap& gatherMap = site.getGatherMap();
#endif
          //replacing values
          const map<int,int>&   globalToLocal = _mesh.getGlobalToLocal();
          foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap){
             const StorageSite&  oSite = *mpos.first;
             const Array<int>& gatherArray = dynamic_cast< const Array<int>& > (*mpos.second);
             //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
             int from_where       = oSite.getGatherProcID();
             if ( from_where != -1 ){
                EntryIndex e(&oSite,&site);
                const Array<int>  & recv_counts  = dynamic_cast< const Array<int>&   > (*_recvCounts [e]);
                const Array<int>  & recv_indices = dynamic_cast< const Array<int>&   > (*_recvIndices[e]);
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


    }

    //resize class member arrays
    void fillClassMemberArrays()
    { 
       //get new row and col array from new connectivity
       const Array<int>&  row = _conn->getRow();
       const Array<int>&  col = _conn->getCol();
       //resize member _row and _col and coeff arrays
       _row.resize   ( row.getLength() );
       _col.resize   ( col.getLength() );
       _coeffs.resize( col.getLength() );
       //assign row and col
       _row = row;
       _col = col;
    }
 

#ifdef FVM_PARALLEL
    void CRConnectivityPrint( const CRConnectivity& conn, int procID, const string& name )
    {
       if ( MPI::COMM_WORLD.Get_rank() == procID ){
          cout <<  name << " :" << endl;
          const Array<int>& row = conn.getRow();
          const Array<int>& col = conn.getCol();
          for ( int i = 0; i < row.getLength()-1; i++ ){
             cout << " i = " << i << ",    ";
             for ( int j = row[i]; j < row[i+1]; j++ )
                cout << col[j] << "  ";
             cout << endl;
             }
       }

     }
#endif

  virtual void printRow(const int i) const
  {
    
    for (int nb = _row[i]; nb<_row[i+1]; nb++)
    {
        const int j = _col[nb];
        cout << "coeff (" <<  i << "," << j << ") = " << _coeffs[nb] << endl;
    }
  }

  const Mesh&   _mesh;
  shared_ptr< CRConnectivity> _conn;
  Array<int> _row;
  Array<int> _col;
  Array<Coord> _coeffs;
  shared_ptr<StorageSite> _newSite;
  mutable map<const CRConnectivity*,PairWiseAssembler*> _pairWiseAssemblers;
  GhostArrayMap   _sendCounts;
  GhostArrayMap   _recvCounts;
  GhostArrayMap   _sendIndices;
  GhostArrayMap   _recvIndices;
  GhostArrayMap   _sendValues;
  GhostArrayMap   _recvValues;

};


#endif
