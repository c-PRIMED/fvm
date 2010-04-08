#ifndef _CRMATRIX_H_
#define _CRMATRIX_H_

#include "Matrix.h"
#include "CRConnectivity.h"
#include "Array.h"
#include "StorageSite.h"
#include "LinearSystemMerger.h"

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

#include <set>

/**
 * Sparse matrix stored using a compressed row format. The sparsity
 * pattern is provided by a CRConnectivity object that is required at
 * construction time. Note that we assume that diagonal is implicitly
 * always present (i.e., col[row[i]] through col[row[i+1]] do not
 * contain i) and store diagonal and off diagonal entries in separate
 * arrays.
 *
 * The class is templated with the types of diagonal, off-diagonal and
 * the corresponding x arrays as template parameters.
 * 
 */

template<class T_Diag, class T_OffDiag, class X>
class CRMatrix : public Matrix
{
public:

   //friend class LinearSystemMerger;

  typedef T_Diag Diag;
  typedef T_OffDiag OffDiag;
  typedef Array<Diag> DiagArray;
  typedef Array<int> IntArray;
  
  typedef Array<OffDiag> OffDiagArray;
  typedef Array<X> XArray;
  typedef shared_ptr< Array<int> >  IntArrayPtr;
  typedef shared_ptr< Array<Diag> > DiagArrayPtr;
  
  typedef shared_ptr< Array<double> >  ArrayDblePtr;
  /**
   * Embedded class used for easy (ie. no search) access to matrix
   * entries for the special case of face based finite volume
   * discretizations. Works in conjunction with a pairwise access
   * mapping provided by CRConnectivity
   * 
   */

  class PairWiseAssembler
  {
  public:
    PairWiseAssembler(Array<OffDiag>& coeffs,
                      const Array<Vector<int,2> >& pairToCol) :
      _coeffs(coeffs),
      _pairToCol(pairToCol)
    {}

    OffDiag& getCoeff01(const int np)
    {
      return _coeffs[_pairToCol[np][0]];
    }
    
    OffDiag& getCoeff10(const int np)
    {
      return _coeffs[_pairToCol[np][1]];
    }
    
    void addCoeffsSymmetric(const int np, const OffDiag& c)
    {
      _coeffs[_pairToCol[np][0]] += c;
      _coeffs[_pairToCol[np][1]] += c;
    }

    void addCoeffs(const int np, const OffDiag& c01, const OffDiag& c10)
    {
      _coeffs[_pairToCol[np][0]] += c01;
      _coeffs[_pairToCol[np][1]] += c10;
    }

    void addCoeff01(const int np, const OffDiag& c01)
    {
      _coeffs[_pairToCol[np][0]] += c01;
    }
    
    void addCoeff10(const int np, const OffDiag& c10)
    {
      _coeffs[_pairToCol[np][1]] += c10;
    }
  private:
    Array<OffDiag>& _coeffs;
    const Array<Vector<int,2> >& _pairToCol;
  };

  typedef map<const CRConnectivity*,PairWiseAssembler*> PairWiseAssemblerMap;
  

  CRMatrix(const CRConnectivity& conn) :
    Matrix(),
    _conn(conn),
    _row(_conn.getRow()),
    _col(_conn.getCol()),
    _diag(_conn.getRowDim()),
    _offDiag(_col.getLength()),
    _pairWiseAssemblers(),
    _isBoundary(_conn.getRowDim())
  {
    logCtor();
    _isBoundary = false;
  }

  
  DEFINE_TYPENAME("CRMatrix<"
                  +NumTypeTraits<Diag>::getTypeName()+","
                  +NumTypeTraits<OffDiag>::getTypeName()+","
                  +NumTypeTraits<X>::getTypeName()
             +">");

  virtual void initAssembly()
  {
    _diag.zero();
    _offDiag.zero();
    _isBoundary = false;
  }

  /**
   * y = this * x
   * 
   */

  virtual void multiply(IContainer& yB, const IContainer& xB) const
  {
    XArray& y = dynamic_cast<XArray&>(yB);
    const XArray& x = dynamic_cast<const XArray&>(xB);
    
    const int nRows = _conn.getRowSite().getSelfCount();
    for(int nr=0; nr<nRows; nr++)
    {
        y[nr] = _diag[nr]*x[nr];
        for (int nb = _row[nr]; nb<_row[nr+1]; nb++)
        {
            const int j = _col[nb];
            y[nr] += _offDiag[nb]*x[j];
        }
    }
  }

  /**
   * y += this * x
   * 
   */

  virtual void multiplyAndAdd(IContainer& yB, const IContainer& xB) const
  {
    XArray& y = dynamic_cast<XArray&>(yB);
    const XArray& x = dynamic_cast<const XArray&>(xB);
    
    const int nRows = _conn.getRowSite().getSelfCount();
    for(int nr=0; nr<nRows; nr++)
    {
        y[nr] += _diag[nr]*x[nr];
        for (int nb = _row[nr]; nb<_row[nr+1]; nb++)
        {
            const int j = _col[nb];
            y[nr] += _offDiag[nb]*x[j];
        }
    }
  }

  /**
   * forward GaussSeidel update for this * x  = b
   * 
   */

  virtual void forwardGS(IContainer& xB, IContainer& bB, IContainer&) const
  {
    XArray& x = dynamic_cast<XArray&>(xB);
    const XArray& b = dynamic_cast<const XArray&>(bB);
    
    const int nRows = _conn.getRowSite().getSelfCount();

    X sum;
    for(int nr=0; nr<nRows; nr++)
    {
        sum = b[nr];
        for (int nb = _row[nr]; nb<_row[nr+1]; nb++)
        {
            const int j = _col[nb];
            sum += _offDiag[nb]*x[j];
        }
        x[nr] = -sum/_diag[nr];
    }
  }
  
  /**
   * reverse GaussSeidel update for this * x  = b
   * 
   */

  virtual void reverseGS(IContainer& xB, IContainer& bB, IContainer&) const
  {
    XArray& x = dynamic_cast<XArray&>(xB);
    const XArray& b = dynamic_cast<const XArray&>(bB);
    
    const int nRows = _conn.getRowSite().getSelfCount();

    X sum;
    for(int nr=nRows-1; nr>=0; nr--)
    {
        sum = b[nr];
        for (int nb = _row[nr]; nb<_row[nr+1]; nb++)
        {
            const int j = _col[nb];
            sum += _offDiag[nb]*x[j];
        }
        x[nr] = -sum/_diag[nr];
    }
  }

  /**
   * Jacobi update for this * x  = b
   * 
   */

  virtual void Jacobi(IContainer& xnewB, const IContainer& xoldB,
                         const IContainer& bB) const
  {
    XArray& xnew = dynamic_cast<XArray&>(xnewB);
    const XArray& xold = dynamic_cast<const XArray&>(xoldB);
    const XArray& b = dynamic_cast<const XArray&>(bB);
    
    const int nRows = _conn.getRowSite().getSelfCount();

    X sum;
    for(int nr=0; nr<nRows; nr++)
    {
        sum = b[nr];
        for (int nb = _row[nr]; nb<_row[nr+1]; nb++)
        {
            const int j = _col[nb];
            sum += _offDiag[nb]*xold[j];
        }
        xnew[nr] = -sum/_diag[nr];
    }
  }

  virtual void iluSolve(IContainer& xB, const IContainer& bB, const IContainer&) const
  {
    XArray& x = dynamic_cast<XArray&>(xB);
    shared_ptr<XArray> y = dynamic_pointer_cast<XArray>(x.newClone());
    const XArray& b = dynamic_cast<const XArray&>(bB);

    if (!_iluConnPtr)
      compute_ILU0();

    lowerSolve(*y,b);
    upperSolve(x,*y);
  }
  
  /**
   * r = b + this*x
   * 
   */

  virtual void computeResidual(const IContainer& xB, const IContainer& bB,
                               IContainer& rB) const
  {
    const XArray& x = dynamic_cast<const XArray&>(xB);
    const XArray& b = dynamic_cast<const XArray&>(bB);
    XArray& r = dynamic_cast<XArray&>(rB);
    
    const int nRows = _conn.getRowSite().getSelfCount();

    for(int nr=0; nr<nRows; nr++)
    {
        r[nr] = b[nr] + _diag[nr]*x[nr];
        for (int nb = _row[nr]; nb<_row[nr+1]; nb++)
        {
            const int j = _col[nb];
            r[nr] += _offDiag[nb]*x[j];
        }
    }
  }

  /**
   * solve rows marked as boundary
   * 
   */

  virtual void solveBoundary(IContainer& xB, IContainer& bB, IContainer&) const
  {
    XArray& x = dynamic_cast<XArray&>(xB);
    const XArray& b = dynamic_cast<const XArray&>(bB);
    
    const int nRowsInterior = _conn.getRowSite().getSelfCount();
    const int nRowsExtra = _conn.getRowSite().getCount();

    X sum;
    for(int nr=nRowsInterior; nr<nRowsExtra; nr++)
      if (_isBoundary[nr])
      {
          sum = b[nr];
          for (int nb = _row[nr]; nb<_row[nr+1]; nb++)
          {
              const int j = _col[nb];
              sum += _offDiag[nb]*x[j];
          }
          x[nr] = -sum/_diag[nr];
      }
  }
  
  /**
   * create a matrix coarsening by agglomerating rows with
   * neighbours that have the largest coefficients. groupSize is the
   * approximate number of fine rows grouped in each coarse row and
   * weightRatioThreshold is used to dtermine if the coefficient is
   * large enough.
   *
   * The coarsening is defined by the coarseIndex array that stores
   * for each fine row the index of the coarse row into which it has
   * been grouped. The function returns the size of the coarse level.
   */

  virtual int createCoarsening(IContainer& gCoarseIndex, const int groupSize,
                               const double weightRatioThreshold)
  {
    Array<int>&  coarseIndex = dynamic_cast<Array<int>& >(gCoarseIndex);

    const int nRows = _conn.getRowSite().getSelfCount();

    coarseIndex = -1;

    // the number of rows in the coarse level matrix
    int nCoarseRows=0;

    // temp storage to keep track of how many in each coarse group
    Array<int> coarseCount(nRows);

    coarseCount = 0;

    for(int nr=0; nr<nRows; nr++)
      if (coarseIndex[nr] == -1 && !_isBoundary[nr])
      {
          // the row that we are going to go through looking for
          // neighbours to group with; it starts off being the current
          // row but for group sizes greater than 2 it will be reset
          // to the ones that we group this one with

          int current = nr;

          int colMaxGrouped=-1; // columnn with largest coeff among grouped ones
          int colMaxUngrouped=-1; // columnn with largest coeff among ungrouped ones

          int nGrouped;

          // assign current row to a new coarse row

          coarseIndex[current] = nCoarseRows;

          for(nGrouped=1; nGrouped<groupSize; nGrouped++)
          {
              // find largest coeff among both grouped and ungrouped columns

              double maxWeightUngrouped = 0;
              double maxWeightGrouped = 0;
              colMaxGrouped = -1;
              colMaxUngrouped = -1;

              for(int nb=_row[current]; nb<_row[current+1]; nb++)
              {
                  const int nc = _col[nb];

                  //skip ghost and boundaries
                  if (nc < nRows  && !_isBoundary[nc])
                  {
                      double diagMeasure0 =
                        NumTypeTraits<Diag>::doubleMeasure(_diag[nr]);
                      double diagMeasure1 =
                        NumTypeTraits<Diag>::doubleMeasure(_diag[nc]);
                      double coeffMeasure =
                        NumTypeTraits<OffDiag>::doubleMeasure(_offDiag[nb]);
                      
                      double thisWeight = fabs(coeffMeasure /
                                               max(diagMeasure0,diagMeasure1));

                      if (coarseIndex[nc] == -1)
                      {
                          if (colMaxUngrouped == -1 || (thisWeight > maxWeightUngrouped))
                          {
                              colMaxUngrouped = nc;
                              maxWeightUngrouped = thisWeight;
                          }
                      }
                      else if (coarseIndex[nc] != coarseIndex[nr])
                      {
                          if (colMaxGrouped == -1 || (thisWeight > maxWeightGrouped))
                          {
                              colMaxGrouped = nc;
                              maxWeightGrouped = thisWeight;
                          }
                      }
                  }
              }

              // if we found at least one ungrouped row, group it
              // with the current one as long as it is large
              // enough compared to already grouped ones
              
              if ( (colMaxUngrouped != -1) &&
                   (colMaxGrouped == -1 ||
                    (maxWeightUngrouped > weightRatioThreshold*maxWeightGrouped)))
              {
                  coarseIndex[colMaxUngrouped] = coarseIndex[current];
                  coarseCount[coarseIndex[current]]++;
                  // next we will look at neighbours of this row
                  current = colMaxUngrouped;
              }
              else
                // no point trying to find more for the current coarse group
                break;
          }
          
          // if we managed to find at least one ungrouped
          // neighbour to group the current row with or if there
          // is no already grouped neighbour that is too crowded
          // then accept the new group
          if (nGrouped > 1 || colMaxGrouped == -1 ||
              coarseCount[coarseIndex[colMaxGrouped]] > groupSize+2)
          {
              coarseCount[coarseIndex[nr]]++;
              nCoarseRows++;
          }
          else
          {
              // put the current row in an existing group
              coarseIndex[nr] = coarseIndex[colMaxGrouped];
              coarseCount[coarseIndex[colMaxGrouped]]++;
          }
      }

    return nCoarseRows;
  }


  /**
   * Create the connectivity for a coarse level matrix given the
   * coarsening (ie. the fine to coarse index mapping) and its
   * transpose (ie coarse to fine mapping, provided as a
   * CRConnectivity object)
   * 
   */

  shared_ptr<CRConnectivity>
  createCoarseConnectivity(const IContainer& gCoarseIndex,
                           const CRConnectivity& coarseToFine,
                           const StorageSite& coarseRowSite,
                           const StorageSite& coarseColSite)
  {
    const Array<int>&  coarseIndex = dynamic_cast<const Array<int>& >(gCoarseIndex);

    const int nCoarseRows = coarseRowSite.getCount();

    shared_ptr<CRConnectivity> coarseCR(new CRConnectivity(coarseRowSite,coarseColSite));


    coarseCR->initCount();

    Array<bool> coarseCounted(nCoarseRows);

    coarseCounted = false;
    
    for(int nrCoarse=0; nrCoarse<nCoarseRows; nrCoarse++)
    {
        const int nFine = coarseToFine.getCount(nrCoarse); 
        for(int nfr=0; nfr<nFine; nfr++)
        {
            const int nrFine = coarseToFine(nrCoarse,nfr);
            for (int nb = _row[nrFine]; nb<_row[nrFine+1]; nb++)
            {
                const int nc = _col[nb];
                const int ncCoarse = coarseIndex[nc];


                if (ncCoarse>=0 && nrCoarse!=ncCoarse && !coarseCounted[ncCoarse])
                {
                    coarseCounted[ncCoarse] = true;
                    coarseCR->addCount(nrCoarse,1);
                }
            }
        }

        // reset counted 
        for(int nfr=0; nfr<nFine; nfr++)
        {
            const int nrFine = coarseToFine(nrCoarse,nfr);
            
            for (int nb = _row[nrFine]; nb<_row[nrFine+1]; nb++)
            {
                const int nc = _col[nb];
                const int ncCoarse = coarseIndex[nc];
                if (ncCoarse>=0)
                  coarseCounted[ncCoarse] = false;
            }
        }
    }

    coarseCR->finishCount();

    for(int nrCoarse=0; nrCoarse<nCoarseRows; nrCoarse++)
    {
        const int nFine = coarseToFine.getCount(nrCoarse);
        for(int nfr=0; nfr<nFine; nfr++)
        {
            const int nrFine = coarseToFine(nrCoarse,nfr);

            for (int nb = _row[nrFine]; nb<_row[nrFine+1]; nb++)
            {
                const int nc = _col[nb];
                const int ncCoarse = coarseIndex[nc];
                if (ncCoarse>=0 && nrCoarse!=ncCoarse && !coarseCounted[ncCoarse])
                {
                    coarseCounted[ncCoarse] = true;
                    coarseCR->add(nrCoarse,ncCoarse);
                }
            }
        }

        // reset counted 
        for(int nfr=0; nfr<nFine; nfr++)
        {
            const int nrFine = coarseToFine(nrCoarse,nfr);
            
            for (int nb = _row[nrFine]; nb<_row[nrFine+1]; nb++)
            {
                const int nc = _col[nb];
                const int ncCoarse = coarseIndex[nc];
                if (ncCoarse>=0)
                coarseCounted[ncCoarse] = false;
            }
        }
    }


    coarseCR->finishAdd();
 
    return coarseCR;
  }

  /**
   * create the coarse matrix given the coarsening and the coarse
   * level connectivity.
   * 
   */

  shared_ptr<Matrix>
  createCoarseMatrix(const IContainer& gCoarseIndex,
                     const CRConnectivity& coarseToFine,
                     const CRConnectivity& coarseConnectivity)
  {
    const Array<int>&  coarseIndex = dynamic_cast<const Array<int>& >(gCoarseIndex);
    const int nCoarseRows = coarseConnectivity.getRowDim();

    shared_ptr<CRMatrix> coarseMatrix(new CRMatrix(coarseConnectivity));

    Array<Diag>& coarseDiag = coarseMatrix->getDiag();
    Array<OffDiag>& coarseOffDiag = coarseMatrix->getOffDiag();

    const Array<int>& coarseConnRow = coarseConnectivity.getRow();
    const Array<int>& coarseConnCol = coarseConnectivity.getCol();

    coarseDiag.zero();
    coarseOffDiag.zero();

    //used to avoid searches when inserting coeffs
    Array<int> coarseCoeffPos(nCoarseRows);


    for(int nrCoarse=0; nrCoarse<nCoarseRows; nrCoarse++)
    {
        // for easy indexing when inserting coefficients set col
        // positions into the coarse connectivity
        for(int nb=coarseConnRow[nrCoarse]; nb<coarseConnRow[nrCoarse+1]; nb++)
          coarseCoeffPos[coarseConnCol[nb]] = nb;

        // loop over the fine rows that make up this coarse row
        const int nFine = coarseToFine.getCount(nrCoarse);
        for(int nfr=0; nfr<nFine; nfr++)
        {
            const int nrFine = coarseToFine(nrCoarse,nfr);

            coarseDiag[nrCoarse] += _diag[nrFine];
            
            for (int nb = _row[nrFine]; nb<_row[nrFine+1]; nb++)
            {
                const int nc = _col[nb];
                const int ncCoarse = coarseIndex[nc];

                if (ncCoarse<0) continue;
                
                if (nrCoarse!=ncCoarse)
                {
                    const int pos = coarseCoeffPos[ncCoarse];
                    coarseOffDiag[pos] += _offDiag[nb];
                }
                else 
                {
                    coarseDiag[nrCoarse] += _offDiag[nb];
                }
            }
        }
    }

    return coarseMatrix;
  }

#ifdef FVM_PARALLEL

shared_ptr<Matrix>
createMergeMatrix( const LinearSystemMerger& mergeLS )
{
   const CRConnectivity&  conn = mergeLS.getConnectivity();

   const map<int,IntArrayPtr>& localToGlobalMap = mergeLS.getLocalToGlobal();

   const Array<int>& globalToProc  = mergeLS.getGlobalToProc();
   const Array<int>& globalToLocal = mergeLS.getGlobalToLocal();
   const Array<int>& selfCounts    = mergeLS.getSelfCounts();

   const map< int, IntArrayPtr >&  rowMap = mergeLS.getLocalConnRow();
   const map< int, IntArrayPtr >&  colMap = mergeLS.getLocalConnCol();
   const vector< map<int,int> >&   gatherIDsLocalToGlobal = mergeLS.getGatherIDsLocalToGlobal();  //[proc][localid] = globalID

   const map<int, ArrayDblePtr> &  diagMap   =  mergeLS.getDiag();
   const map<int, ArrayDblePtr> & offDiagMap =  mergeLS.getOffDiag();

   const MPI::Intracomm& comm = mergeLS.getComm();
   const set<int>& group = mergeLS.getGroup();

   shared_ptr<CRMatrix> mergeMatrix( new CRMatrix(conn) );
   mergeMatrix->initAssembly();

   Array<Diag>&    mergeDiag    = mergeMatrix->getDiag();
   Array<OffDiag>& mergeOffDiag = mergeMatrix->getOffDiag();
   const Array<int>& mergeRow = conn.getRow();
   const Array<int>& mergeCol = conn.getCol();
   foreach ( const set<int>::value_type proc_id, group ){
      const Array<double>& diag    = *diagMap.find(proc_id)->second;
      const Array<double>& offDiag = *offDiagMap.find(proc_id)->second;

      const Array<int>& localToGlobal = *localToGlobalMap.find(proc_id)->second;
      const map<int,int>& gatherIDsMap = gatherIDsLocalToGlobal[proc_id];

      const Array<int>& row = *rowMap.find(proc_id)->second;
      const Array<int>& col = *colMap.find(proc_id)->second;

      for ( int i = 0; i < selfCounts[proc_id]; i++ ){
          int glbl_indx = localToGlobal[i];
          mergeDiag[glbl_indx] = diag[i];
      }

      for ( int i = 0; i < selfCounts[proc_id]; i++ ){
           int global_indx = localToGlobal[i];
           for ( int np = row[i]; np < row[i+1]; np++ ){ 
                int local_j = col[np];

                //if -1 means boundary skip it
                if ( local_j >= 0 ){
                    int global_j = -1;
                    if ( local_j < selfCounts[proc_id] )
                        global_j = localToGlobal[ local_j ];
                    else  //ghost 
                        global_j = gatherIDsMap.find(local_j)->second;


                     for ( int npg = mergeRow[global_indx]; npg < mergeRow[global_indx+1]; npg++ ){
                         if ( mergeCol[npg] == global_j ){ 
                             mergeOffDiag[npg] = offDiag[np];
                         }
                     }
                }
           }	
      }

   }


    return mergeMatrix;

}

#endif
  
  virtual const CRConnectivity& getConnectivity() const {return _conn;}

  Array<Diag>& getDiag() {return _diag;}
  Array<OffDiag>& getOffDiag() {return _offDiag;}

  const Array<Diag>& getDiag() const {return _diag;}
  const Array<OffDiag>& getOffDiag() const {return _offDiag;}
 
  void *getDiagData() const { return _diag.getData(); }
  void *getOffDiagData() const { return _offDiag.getData(); }
  int   getDiagDataSize() const { return _diag.getDataSize(); }
  int   getOffDiagDataSize() const { return _offDiag.getDataSize(); }
  
  PairWiseAssembler& getPairWiseAssembler(const CRConnectivity& pairs)
  {
    if (_pairWiseAssemblers.find(&pairs) == _pairWiseAssemblers.end())
    {
        _pairWiseAssemblers[&pairs] =
          new PairWiseAssembler(_offDiag,
                                _conn.getPairToColMapping(pairs));
    }
    return *_pairWiseAssemblers[&pairs];
  }

  virtual ~CRMatrix()
  {
    for (typename PairWiseAssemblerMap::iterator pos = _pairWiseAssemblers.begin();
         pos != _pairWiseAssemblers.end();
         ++pos)
      delete pos->second;
    logDtor();
  }

  void setDirichlet(const int nr)
  {
    for (int nb = _row[nr]; nb<_row[nr+1]; nb++)
      _offDiag[nb] = NumTypeTraits<OffDiag>::getZero();
    _diag[nr] = NumTypeTraits<Diag>::getNegativeUnity();
    _isBoundary[nr] = true;
  }
  
  void setBoundary(const int nr)
  {
    _isBoundary[nr] = true;
  }


private:

  // computes the level 0 (i.e. with no additional non zero entries)
  // incomplete LU factorization
  void compute_ILU0() const
  {
    const int nRows = _conn.getRowSite().getSelfCount();

    // create connectivity matrix for the ILU factored matrix; for ILU
    // 0 the only difference from this CRMatrix's connectivity is that
    // the diagonal is explicitly included
    //
    
    _iluConnPtr = shared_ptr<CRConnectivity>
      (new CRConnectivity(_conn.getRowSite(),_conn.getColSite()));
    CRConnectivity& iluConn = *_iluConnPtr;

    iluConn.initCount();

    for(int nr=0; nr<nRows; nr++)
    {
        // add one for diagonal
        iluConn.addCount(nr,1);
        for(int nb=_row[nr]; nb<_row[nr+1]; nb++)
        {
            // only include neighbours in this matrix
            if (_col[nb] < nRows)
              iluConn.addCount(nr,1);
        }
    }

    iluConn.finishCount();

    const int nnz = iluConn.getCol().getLength();
    _iluCoeffsPtr = DiagArrayPtr(new DiagArray(nnz));
    _iluDiagIndexPtr = IntArrayPtr(new IntArray(nRows));

    DiagArray& iluCoeffs = *_iluCoeffsPtr;
    IntArray& iluDiagIndex = *_iluDiagIndexPtr;
    
    // now we can fill the iluCoeffs with this matrix's coeffs
    for(int nr=0; nr<nRows; nr++)
    {
        // lower coeffs first
        for(int nb=_row[nr]; nb<_row[nr+1]; nb++)
        {
            const int j = _col[nb];
            if (j < nRows && j < nr)
            {
                const int pos = iluConn.add(nr,j);
                iluCoeffs[pos] = _offDiag[nb];
            }

        }

        // add diagonal next
        const int pos = iluConn.add(nr,nr);
        iluCoeffs[pos] = _diag[nr];
        iluDiagIndex[nr] = pos;

        // finally the upper coeffs 
        for(int nb=_row[nr]; nb<_row[nr+1]; nb++)
        {
            const int j = _col[nb];
            if (j < nRows && j > nr)
            {
                const int pos = iluConn.add(nr,j);
                iluCoeffs[pos] = _offDiag[nb];
            }

        }
    }

    // finalize the iluConnectivity
    iluConn.finishAdd();

    // ilu decomposition

    const IntArray& iluRow = iluConn.getRow();
    const IntArray& iluCol = iluConn.getCol();
    
    
    //work array iw(n) and diagonal pointer uptr
    Array<int> iw(nRows);
    Array<int> uptr(nRows);
    iw = 0;
    uptr = 0;

    const Diag unity(NumTypeTraits<Diag>::getUnity());
    
    //main loop
    for (int k = 0; k < nRows; k++ )
    {
        const int j1 = iluRow[k];
        const int j2 = iluRow[k+1];
        for (int j = j1; j <j2; j++ )
          iw[ iluCol[j] ] = j;
	  
        int j = j1;
        do
        {
            int jrow = iluCol[j];
            if ( jrow < k )
            {
                const Diag t1 = iluCoeffs[j] * iluCoeffs[ uptr[jrow] ];
                iluCoeffs[j] = t1;
                //perform linear combination
                for (int jj = uptr[jrow]+1; jj < iluRow[jrow+1]; jj++ )
                {
                    int jw = iw[  iluCol[jj] ];
                    if ( jw != 0 )
		         iluCoeffs[jw] -= t1 * iluCoeffs[jj];
                }
                j++;
            }
            else
            {
	        uptr[k] = j;
		break;
            }
        }
        while( j < j2 );

        //inverting diagonal
        iluCoeffs[j]  = unity / iluCoeffs[j];	     
        //zeroing iw
        for ( int i = j1; i < j2; i++ )
          iw[ iluCol[i] ] = 0;
    }  	    
  }
  
  
  //Solve L y = b. This assumes that the factorization has already
  //been performed
  void lowerSolve( XArray& y,  const XArray& b) const
  {
    const IntArray& iluRow = _iluConnPtr->getRow();
    const IntArray& iluCol = _iluConnPtr->getCol();
    const IntArray& iluDiagIndex = *_iluDiagIndexPtr;
    const DiagArray& iluCoeffs = *_iluCoeffsPtr;

    const int nRows = _conn.getRowSite().getSelfCount();

    for (int j = 0; j < nRows; j++)
    {
        X yj = -b[j];
        for ( int k = iluRow[j]; k < iluDiagIndex[j]; k++ )
        {
            yj -= iluCoeffs[k] * y[ iluCol[k] ];
        }
        y[j] = yj;
    }
  }
  
  //solve U x = y
  void upperSolve( XArray& x,  const XArray& y ) const
  {
    const IntArray& iluRow = _iluConnPtr->getRow();
    const IntArray& iluCol = _iluConnPtr->getCol();
    const IntArray& iluDiagIndex = *_iluDiagIndexPtr;
    const DiagArray& iluCoeffs = *_iluCoeffsPtr;
    
    const int nRows = _conn.getRowSite().getSelfCount();

    for ( int j = nRows-1; j >= 0; j-- )
    {
        X xj = y[j];
        for ( int k =  iluDiagIndex[j]+1; k <iluRow[j+1]; k++ )
        {
            xj -= iluCoeffs[k] * x[ iluCol[k] ];
        }
        x[j] = iluCoeffs[ iluDiagIndex[j] ]*xj;
    }
  }
  

  const CRConnectivity& _conn;
  const Array<int>& _row;
  const Array<int>& _col;
  Array<Diag> _diag;
  Array<OffDiag> _offDiag;
  PairWiseAssemblerMap _pairWiseAssemblers;
  Array<bool> _isBoundary;

  // connectivity for the ILU factorization. This connectivity
  // explicitly includes diagonal and is also ordered such that for a
  // particular row the lower triangle entries are first, followed by
  // the diagonal and then the upper triangle entries. The _iluDiagIndexPtr
  // array is used to keep the index of the diagonal location
  mutable shared_ptr<CRConnectivity> _iluConnPtr;

  mutable IntArrayPtr _iluDiagIndexPtr;

  // ilu coeffs
  mutable DiagArrayPtr _iluCoeffsPtr;

 
};


#include "Vector.h"
#include "DiagonalTensor.h"

template<class T>
struct MatrixTraitHelper
{
  typedef CRMatrix<T,T,T> T_Matrix;
};

template<class T, int N>
struct MatrixTraitHelper<Vector<T,N> >
{
  typedef Vector<T,N> T_Vector;
  typedef DiagonalTensor<T,N> T_DiagTensor;
  typedef CRMatrix<T_DiagTensor,T_DiagTensor,T_Vector> T_Matrix;
};

#endif
