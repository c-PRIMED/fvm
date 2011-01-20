#ifndef _GRADIENTMATRIX_H_
#define _GRADIENTMATRIX_H_

#include "CRConnectivity.h"
#include "Array.h"
#include "StorageSite.h"
#include "Gradient.h"

class GradientMatrixBase
{
public:
  GradientMatrixBase() {}
  virtual ~GradientMatrixBase() {}
};

template<class T_Scalar>
class GradientMatrix : public GradientMatrixBase
{
public:
  typedef Vector<T_Scalar,3> Coord;
  
  GradientMatrix(const CRConnectivity& conn) :
    _conn(conn),
    _row(_conn.getRow()),
    _col(_conn.getCol()),
    _coeffs(_col.getLength()),
    _pairWiseAssemblers()
  {}

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
    
    const int nRows = _conn.getRowSite().getSelfCount();
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
                                _conn.getPairToColMapping(pairs));
    }
    return *_pairWiseAssemblers[&pairs];
  }
  
private:
  const CRConnectivity& _conn;
  const Array<int>& _row;
  const Array<int>& _col;
  Array<Coord> _coeffs;
  mutable map<const CRConnectivity*,PairWiseAssembler*> _pairWiseAssemblers;
};


#endif
