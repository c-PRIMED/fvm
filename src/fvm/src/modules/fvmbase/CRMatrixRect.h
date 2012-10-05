// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _CRMATRIXRECT_H_
#define _CRMATRIXRECT_H_

#include "Matrix.h"
#include "CRConnectivity.h"
#include "Array.h"
#include "StorageSite.h"


/**
 * Sparse matrix stored using a compressed row format. The sparsity
 * pattern is provided by a CRConnectivity object that is required at
 * construction time. Note that we assume that diagonal is implicitly
 * always present (i.e., col[row[i]] through col[row[i+1]] do not
 * contain i) and store diagonal and off diagonal entries in separate
 * arrays.
 *
 * The class is templated with the types of diagonal, off-diagonal and
 * the corresponsing x arrays as template parameters.
 * 
 */

template<class T_Coeff, class X, class B>
class CRMatrixRect : public Matrix
{
public:
  typedef T_Coeff Coeff;
  typedef T_Coeff Diag;
  typedef T_Coeff OffDiag;
  typedef Array<Coeff> CoeffArray;
  typedef Array<Coeff> DiagArray;
  typedef Array<X> XArray;
  typedef Array<B> BArray;

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
  

  CRMatrixRect(const CRConnectivity& conn) :
    Matrix(),
    _conn(conn),
    _row(_conn.getRow()),
    _col(_conn.getCol()),
    _diag(_conn.getRowDim()),
    _offDiag(_col.getLength()),
    _pairWiseAssemblers()
  {
    logCtor();
  }

  
  DEFINE_TYPENAME("CRMatrixRect<"
             +NumTypeTraits<Coeff>::getTypeName()+","
             +NumTypeTraits<X>::getTypeName()+","
             +NumTypeTraits<B>::getTypeName()
             +">");

  virtual void initAssembly()
  {
    _diag.zero();
    _offDiag.zero();
  }

  /**
   * y = this * x
   * 
   */

  virtual void multiply(IContainer& yB, const IContainer& xB) const
  {
    BArray& y = dynamic_cast<BArray&>(yB);
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
    BArray& y = dynamic_cast<BArray&>(yB);
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

  
  const CRConnectivity& getConnectivity() const {return _conn;}

  Array<Coeff>& getDiag() {return _diag;}
  Array<Coeff>& getOffDiag() {return _offDiag;}

  const Array<Coeff>& getDiag() const {return _diag;}
  const Array<Coeff>& getOffDiag() const {return _offDiag;}

  virtual ~CRMatrixRect()
  {
    logDtor();
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
    const Array<int>&  coarseIndex =
      dynamic_cast<const Array<int>& >(gCoarseIndex);

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

  virtual
  shared_ptr<Matrix>
  createCoarseMatrix(const IContainer& gCoarseIndex,
                     const CRConnectivity& coarseToFine,
                     const CRConnectivity& coarseConnectivity)
  {
    const Array<int>&  coarseIndex =
      dynamic_cast<const Array<int>& >(gCoarseIndex);
    
    const int nCoarseRows = coarseConnectivity.getRowDim();

    shared_ptr<CRMatrixRect> coarseMatrix(new CRMatrixRect(coarseConnectivity));

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


private:
  const CRConnectivity& _conn;
  const Array<int>& _row;
  const Array<int>& _col;
  Array<Coeff> _diag;
  Array<Coeff> _offDiag;
  PairWiseAssemblerMap _pairWiseAssemblers;
};


#endif
