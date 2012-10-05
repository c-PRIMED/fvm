// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _SQUARETENSOR_H_
#define _SQUARETENSOR_H_

#include "NumType.h"
#include "Vector.h"
#include <sstream>
#include "CRConnectivity.h"


template <class T, int N>
class SquareTensor
{
public:

  enum { NSQR = N*N };
  
  typedef SquareTensor<T,N> This_T;
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef typename NumTypeTraits<T>::T_BuiltIn T_BuiltIn;
  
  SquareTensor()
  {}

  SquareTensor(const SquareTensor& o)
  {
    for(int i=0; i<NSQR; i++)
      _data[i] = o._data[i];
  }
  
  SquareTensor(const T& s)
  {
    for(int i=0; i<N; i++)
      for(int j=0; j<N; j++)
        (*this)(i,j) = (i==j) ? s : 0;
  }
  
  
  static string getTypeName()
  {
    return "SquareTensor<" + NumTypeTraits<T>::getTypeName() +
      "," + intAsString(N) +
      ">";
  }
  
  static int getDimension() {return 1;}
  
  static void getShape(int *shp) { *shp = NSQR;}
  static int getDataSize()
  {
    return NSQR*NumTypeTraits<T>::getDataSize();
  }
  
  //  T& operator[](int n) {return _data[n];}
  T& operator()(int i, int j) {return _data[i*N+j];}
  
  //const T& operator[](int n) const {return _data[n];}
  
  const T& operator()(int i, int j) const {return _data[i*N+j];}
    
  void printFromC(ostream &os) const
  {
    os << "[ " ;
    for(int i=0;i<NSQR;i++)
      os << _data[i] << " " ;
    os << "]";
  }

  static void write(FILE* fp, const SquareTensor& x)
  {
    for(int i=0; i<NSQR; i++)
    {
        NumTypeTraits<T>::write(fp,x[i]);
        fprintf(fp, " ");
    }
  }
  
  SquareTensor& operator=(const T& s)
  {
    for(int i=0; i<N; i++)
      for(int j=0; j<N; j++)
        (*this)(i,j) = (i==j) ? s : T(0.);
    return *this;
  }
  
  SquareTensor& operator=(const SquareTensor& o)
  {
    for(int i=0;i<NSQR;i++)
      _data[i] = o._data[i];
    return *this;
  }

  SquareTensor operator-()
  {
    SquareTensor r;
    for(int i=0;i<NSQR;i++)
      r._data[i]=-_data[i];
    return r;
  }

  SquareTensor& operator+=(const SquareTensor& o)
  {
    for(int i=0;i<NSQR;i++)
      _data[i] += o._data[i];
    return *this;
  }

  SquareTensor& operator+=(const T s)
  {
    for(int i=0;i<N;i++)
      (*this)(i,i) += s;
    return *this;
  }

  SquareTensor& operator-=(const SquareTensor& o)
  {
    for(int i=0;i<NSQR;i++)
      _data[i] -= o._data[i];
    return *this;
  }

  SquareTensor& operator-=(const T s)
  {
    for(int i=0;i<N;i++)
      (*this)(i,i) -= s;
    return *this;
  }

  SquareTensor& operator/=(const T s)
  {
    for(int i=0;i<NSQR;i++)
      _data[i] /= s;
    return *this;
  }

  SquareTensor& operator/=(const SquareTensor& o)
  {
    *this *= inverse(o);
    return *this;
  }

  SquareTensor& operator*=(const T s)
  {
    for(int i=0;i<NSQR;i++)
      _data[i] *= s;
    return *this;
  }

  SquareTensor& operator*=(const SquareTensor& o)
  {
    SquareTensor p;
    for(int i=0;i<N;i++)
      for(int j=0;j<N;j++)
      {
          p(i,j) = 0;
          for(int k=0;k<N;k++)
            p(i,j) += (*this)(i,k) * o(k,j);
      }
    
    *this = p;
    return *this;
  }

  void zero()
  {
    for(int i=0;i<NSQR;i++) _data[i] = NumTypeTraits<T>::getZero();
  }
  
  T mag2() const
  {
    T r(NumTypeTraits<T>::getZero());
    for(int i=0; i<N; i++)
      r += (*this)(i,i) * (*this)(i,i);
    return r;
  }

  bool operator<(const double tolerance) const
  {
    return mag2() < tolerance*tolerance;
  }

  static SquareTensor getZero()
  {
    SquareTensor z;
    z.zero();
    return z;
  }

  static double doubleMeasure(const SquareTensor& x)
  {
    double m=0;
    for (int i=0; i<N; i++)
      m += NumTypeTraits<T>::doubleMeasure(x(i,i));
    return m;
  }

  static SquareTensor getNegativeUnity()
  {
    SquareTensor n(getZero());
    for(int i=0; i<N; i++)
      n(i,i) = NumTypeTraits<T>::getNegativeUnity();
    
    return n;
  }

  static SquareTensor getUnity()
  {
    SquareTensor n(getZero());
    for(int i=0; i<N; i++)
      n(i,i) = NumTypeTraits<T>::getUnity();
    
    return n;
  }

  static void accumulateOneNorm(SquareTensor& sum, const SquareTensor& v)
  {
    for(int i=0; i<NSQR; i++)
      NumTypeTraits<T>::accumulateOneNorm(sum._data[i],v._data[i]);
  }

  static void accumulateDotProduct(SquareTensor& sum, const SquareTensor& v0,
                                   const SquareTensor& v1)
  {
    for(int i=0; i<NSQR; i++)
      NumTypeTraits<T>::accumulateDotProduct(sum._data[i],v0._data[i],v1._data[i]);
  }

  static void reduceSum(T_Scalar& sum, const This_T& x)
  {
    for(int i=0; i<NSQR; i++)
      NumTypeTraits<T>::reduceSum(sum,x._data[i]);
  }
  
  static void safeDivide(SquareTensor& x, const SquareTensor& y)
  {
    for(int i=0; i<NSQR; i++)
      NumTypeTraits<T>::safeDivide(x._data[i],y._data[i]);
  }

  static void normalize(SquareTensor& x, const SquareTensor& y)
  {
    for(int i=0; i<NSQR; i++)
      NumTypeTraits<T>::normalize(x._data[i],y._data[i]);
  }

  static void setMax(SquareTensor& x, const SquareTensor& y)
  {
    for(int i=0; i<NSQR; i++)
      NumTypeTraits<T>::setMax(x._data[i],y._data[i]);
  }
  
  private:
  T _data[NSQR];
  };

template<class T, int N>
inline ostream& operator<<(ostream &os,
                           const SquareTensor<T,N> &v)
{
  v.printFromC(os);
  return os;
}


template<class T, int N>
SquareTensor<T,N>
operator+(const SquareTensor<T,N>& a, const SquareTensor<T,N>& b)
{
  return SquareTensor<T,N>(a) += b;
}
  
  template<class T, int N>
SquareTensor<T,N>
operator-(const SquareTensor<T,N>& a, const SquareTensor<T,N>& b)
{
  return SquareTensor<T,N>(a) -= b;
}

template<class T, int N>
SquareTensor<T,N>
operator-(const SquareTensor<T,N>& a)
{
  return -SquareTensor<T,N>(a);
}

template<class T, int N>
SquareTensor<T,N>
operator*(const SquareTensor<T,N>& a, const SquareTensor<T,N>& b)
{
  return SquareTensor<T,N>(a) *= b;
}

template<class T, int N>
SquareTensor<T,N>
operator*(const T s, const SquareTensor<T,N>& a)
{
  return SquareTensor<T,N>(a) *= s;
}

template<class T, int N>
SquareTensor<T,N>
operator*(const SquareTensor<T,N>& a, const T s)
{
  return SquareTensor<T,N>(a) *= s;
}

template<class T, int N>
Vector<T,N>
operator*(const SquareTensor<T,N>& a, const Vector<T,N>& b)
{
  Vector<T,N> r;
  for(int i=0; i<N; i++)
  {
      r[i] = 0;
      for(int j=0; j<N; j++)
        r[i] += a(i,j)*b[j];
  }
  return r;
}

template<class T, int N>
SquareTensor<T,N>
operator/(const SquareTensor<T,N>& a, const T s)
{
  return SquareTensor<T,N>(a) /= s;
}

template<class T, int N>
Vector<T,N>
operator/(const Vector<T,N>& a, const SquareTensor<T,N>& b) 
{
  return inverse(b)*a;
}

template<class T, int N>
SquareTensor<T,N>
operator/(const SquareTensor<T,N>& a, const SquareTensor<T,N>& b) 
{

  return inverse(b)*a;
}

template<class T, int N>
SquareTensor<T,N>
operator/(const T s, const SquareTensor<T,N>& a)
{
  SquareTensor<T,N> r(0);
  for(int i=0; i<N; i++) r(i,i) = s;
  return inverse(a)*r;
}


template<class T>
SquareTensor<T,2>
inverse(const SquareTensor<T,2>& a)
{
  SquareTensor<T,2> inv;
  T det = a(0,0)*a(1,1)-a(0,1)*a(1,0);
  inv(0,0) = a(1,1) / det;
  inv(0,1) = -a(0,1) / det;
  inv(1,0) = -a(1,0) / det;
  inv(1,1) = a(0,0) / det;

  return inv;
}


template<class T>
SquareTensor<T,3>
inverse(const SquareTensor<T,3>& a)
{
  SquareTensor<T,3> inv;
  T det = a(0,0)*(a(1,1)*a(2,2)-a(1,2)*a(2,1))
    -a(0,1)*(a(1,0)*a(2,2) - a(1,2)*a(2,0))
    +a(0,2)*(a(1,0)*a(2,1) - a(1,1)*a(2,0));
  
  inv(0,0) =  (a(1,1)*a(2,2) - a(1,2)*a(2,1)) / det;
  inv(0,1) =  (a(0,2)*a(2,1) - a(0,1)*a(2,2)) / det;
  inv(0,2) =  (a(0,1)*a(1,2) - a(0,2)*a(1,1)) / det;
  inv(1,0) =  (a(1,2)*a(2,0) - a(1,0)*a(2,2)) / det;
  inv(1,1) =  (a(0,0)*a(2,2) - a(0,2)*a(2,0)) / det;
  inv(1,2) =  (a(0,2)*a(1,0) - a(0,0)*a(1,2)) / det;
  inv(2,0) =  (a(1,0)*a(2,1) - a(1,1)*a(2,0)) / det;
  inv(2,1) =  (a(0,1)*a(2,0) - a(0,0)*a(2,1)) / det;
  inv(2,2) =  (a(0,0)*a(1,1) - a(0,1)*a(1,0)) / det;

  return inv;
  
}


template<int N>
void setFlatCoeffs(Array<double>& flatCoeffs,
                   const CRConnectivity& flatConnectivity,
                   const Array<SquareTensor<double,N> >& diag,
                   const Array<SquareTensor<double,N> >& offDiag,
                   const CRConnectivity& connectivity)
{
    const Array<int>& myRow = connectivity.getRow();
    const Array<int>& myCol = connectivity.getCol();
    const int rowDim = connectivity.getRowDim();
    
    for(int i=0; i<rowDim; i++)
    {
        for(int ndr=0; ndr<N; ndr++)
          for(int ndc=0; ndc<N; ndc++)
          {
              const int nfr = i*N + ndr;
              const int nfc = i*N + ndc;

              const int fp = flatConnectivity.getCoeffPosition(nfc,nfr);
              flatCoeffs[fp] = diag[i](ndr,ndc);
          }
      
        for(int jp=myRow[i]; jp<myRow[i+1]; jp++)
        {
            const int j = myCol[jp];
            
            for(int ndr=0; ndr<N; ndr++)
              for(int ndc=0; ndc<N; ndc++)
              {
                  const int nfr = i*N + ndr;
                  const int nfc = j*N + ndc;
                  const int fp = flatConnectivity.getCoeffPosition(nfc,nfr);
                  flatCoeffs[fp] = offDiag[jp](ndr,ndc);
            }
      }
  }
}



/*
template<class T>
SquareTensor<T,10>
  ludcmp(const SquareTensor<T,10>a, int n, Vector<T,10> indx, double d){
  const int NMAX(500);
  const double TINY(1.0e-20);
  int i,imax,j,k;
  double aamax,dum,sum,;
  Array<T> vv;  // vv stores the implicit scaling of each row.
  d=1.0;
  for( i=0;i<n;i++){ 
    aamax=0.0;
    for (j=0;j<n;j++){    
      if (abs(a(i,j)) > aamax) {aamax=abs(a(i,j)); }
    }
    if (aamax == 0.) {break; }  //’singular matrix in ludcmp’ 
    vv(i)=1./aamax ;
    
  }
  
  for( j=0;j<n;j++)  {
    for( i=0;i<j;i++){
      sum=a(i,j);
      for(k=0;k<i;k++){
	sum=sum-a(i,k)*a(k,j);
      }
      a(i,j)=sum;
    }
    aamax=0.; 
    for(i=j;i<n;i++){ 
      sum=a(i,j); 
      for(k=0;k<j;k++){       
	sum=sum-a(i,k)*a(k,j);
      }
      a(i,j)=sum;
      dum=vv(i)*abs(sum); 
      if (dum >= aamax){imax=i;aamax=dum;}
    }
    if (j ~= imax){ //Do we need to interchange rows?    
      for(k=0;k<n;k++){ //Yes, do so...
	dum=a(imax,k);
	a(imax,k)=a(j,k);
	a(j,k)=dum;
      }
      d=-d; //...and change the parity of d.
      vv(imax)=vv(j); //! Also interchange the scale factor.
    }
    indx(j)=imax;
    if(a(j,j) == 0.){a(j,j)=TINY;}
    //If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
    // For some applications on singular matrices, it is desirable to substitute TINY for zero. 
    if(j.~= n){ //Now, finally, divide by the pivot element. 
      dum=1./a(j,j); 
      for( i=j+1;i<n;i++) {a(i,j)=a(i,j)*dum;}
    }
  } //Go back for the next column in the reduction. 
  
  
}


void lubksb(const SquareTensor<T,10>a, int n, Array<T,10> indx, Array<T,10> b)
{
  Array<Int>& indx;
  REAL a(n,n),b(n);
  int i,ii,j,ll;
  double sum;
  ii=0;
  for(i=0;i<n;i++){
    ll=indx(i);
    sum=b(ll);
    b(ll)=b(i);
    if (ii.ne.0){
      for(j=ii;j<i;j++){
	sum=sum-a(i,j)*b(j);
      }
    }
    else if (sum ~= 0.){ 
      ii=i;    //A nonzero element was encountered, so from now on we will
    }          //to do the sums in the loop above.
    b(i)=sum;
  }
  for( i=n-1;i>=0;i--)
    {//backsubstitution, equation (2.3.7).
    sum=b(i);
    for(j=i+1;j<n;j++){
      sum=sum-a(i,j)*b(j);
    }
    b(i)=sum/a(i,i); //Store a component of the solution vector X.
  }
}
*/
  
#endif

