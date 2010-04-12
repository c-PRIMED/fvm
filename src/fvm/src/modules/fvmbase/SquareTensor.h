#ifndef _SQUARETENSOR_H_
#define _SQUARETENSOR_H_

#include "NumType.h"
#include "Vector.h"
#include <sstream>


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
        (*this)(i,i) = (i==j) ? s : 0;
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
  T& operator()(int i, int j) {return _data[j*N+i];}
  
  //const T& operator[](int n) const {return _data[n];}
  
  const T& operator()(int i, int j) const {return _data[j*N+i];}
    
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
        (*this)(i,i) = (i==j) ? s : 0;
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
    SquareTensor n(0);
    for(int i=0; i<N; i++)
      n(i,i) = NumTypeTraits<T>::getNegativeUnity();
    
    return n;
  }

  static SquareTensor getUnity()
  {
    SquareTensor n(0);
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
    -a(0,1)*(a(1,0)*a(2,2) - a(1,2)*a(2,1))
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
  
#endif

