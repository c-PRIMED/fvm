// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _ELECOFFDIAGONALTENSOR_H_
#define _ELECOFFDIAGONALTENSOR_H_

#include "NumType.h"
#include "Vector.h"
#include <sstream>
#include "ElecDiagonalTensor.h"
//this structure is used to store the off dianogal tensor in dielectric charging model
//the full diag tensor looks like
/*
  | 0  0  0 ... 0 |
  | 0  0  0 ... 0 |
  | 0  0  0 ... 0 |
  | 0  0  0 ... Nc |
*/
// only one element is non-zero. it connects to the diagonal tensor via drift model
// thus, a scalar is stored to represent the whole tensor


template <class T, int N>
class ElecOffDiagonalTensor
{
public:
  enum { TN = N+1 } ;
  typedef ElecOffDiagonalTensor<T,N> This_T;
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef typename NumTypeTraits<T>::T_BuiltIn T_BuiltIn;
  
  ElecOffDiagonalTensor()
  {}

  ElecOffDiagonalTensor(const ElecOffDiagonalTensor& o)
  {
    _data = o._data;
  }

  ElecOffDiagonalTensor(const T& o)
  {
    _data = o;
  }

  static string getTypeName()
  {
    return "ElecOffDiagonalTensor<" + NumTypeTraits<T>::getTypeName() +
      "," + intAsString(N) +
      ">";
  }

  static int getDimension() {return 1;}

  static void getShape(int *shp) { *shp = TN;}

  static int getDataSize()
  {
    return NumTypeTraits<T>::getDataSize();
  }

  T& operator[](int n) {return _data;}

  const T& operator[](int n) const {return _data;}

  void printFromC(ostream &os) const
  {
    os << "[ "  << _data << " " << "]";
  }

  static void write(FILE* fp, const ElecOffDiagonalTensor& x)
  {
    NumTypeTraits<T>::write(fp,x._data);
  }
  
  ElecOffDiagonalTensor& operator=(const T& o)
  {
    _data = o;
    return *this;
  }

  ElecOffDiagonalTensor& operator=(const ElecOffDiagonalTensor& o)
  {
    _data = o._data;
    return *this;
  }

  ElecOffDiagonalTensor operator-()
  {
    ElecOffDiagonalTensor r;
    r._data=-_data;
    return r;
  }

  ElecOffDiagonalTensor& operator+=(const ElecOffDiagonalTensor& o)
  {
    _data += o._data;
    return *this;
  }

  ElecOffDiagonalTensor& operator+=(const T s)
  {
    _data += s;
    return *this;
  }

  ElecOffDiagonalTensor& operator-=(const ElecOffDiagonalTensor& o)
  {
    _data -= o._data;
    return *this;
  }

  ElecOffDiagonalTensor& operator-=(const T s)
  {
    _data -= s;
    return *this;
  }

  ElecOffDiagonalTensor& operator/=(const T s)
  {
    throw CException("operator not defined for elec offdiag /= s");
    _data /= s;
    return *this;
  }

 
  ElecOffDiagonalTensor& operator/=(const ElecOffDiagonalTensor& o)
  {
    throw CException("no operator defined for /= elec offdiag");
  }
  

  ElecOffDiagonalTensor& operator*=(const T s)
  {
    throw CException("operator not defined for elec offdiag *= s");
    _data *= s;
    return *this;
  }

  
  ElecOffDiagonalTensor& operator*=(const ElecOffDiagonalTensor& o)
  {
    throw CException("no operator defined for *= elec offdiag");
  }
  
  void zero()
  {
    _data = NumTypeTraits<T>::getZero();
  }

  T mag2() const
  {
    T r(NumTypeTraits<T>::getZero());
    r+=_data * _data;
    return r;
  }

  bool operator<(const double tolerance) const
  {
    return mag2() < tolerance*tolerance;
  }

  static ElecOffDiagonalTensor getZero()
  {
    ElecOffDiagonalTensor z;
    z.zero();
    return z;
  }

  static double doubleMeasure(const ElecOffDiagonalTensor& x)
  {
     return 0.0;  
  }

  static ElecOffDiagonalTensor getNegativeUnity()
  {
    ElecOffDiagonalTensor n;
    n._data = NumTypeTraits<T>::getNegativeUnity();
    
    return n;
  }

  static ElecOffDiagonalTensor getUnity()
  {
    ElecOffDiagonalTensor n;
    n._data = NumTypeTraits<T>::getUnity();
    
    return n;
  }

  static void accumulateOneNorm(ElecOffDiagonalTensor& sum, const ElecOffDiagonalTensor& v)
  {
    throw CException("elec offdiag accumlateOneNorm is not defined!");
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::accumulateOneNorm(sum[i],v[i]);
  }

  static void accumulateDotProduct(ElecOffDiagonalTensor& sum, const ElecOffDiagonalTensor& v0,
                                const ElecOffDiagonalTensor& v1)
  {
    throw CException("elec offdiag accumlateDotProduct is not defined!");
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::accumulateDotProduct(sum[i],v0[i],v1[i]);
  }

  static void reduceSum(T_Scalar& sum, const This_T& x)
  {
    throw CException(" elec offdiag reduceSume is not defined!");
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::reduceSum(sum,x[i]);
  }
  
  static void safeDivide(ElecOffDiagonalTensor& x, const ElecOffDiagonalTensor& y)
  {
    throw CException(" elec offdiag safeDivide is not defined!");
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::safeDivide(x[i],y[i]);
  }

  static void normalize(ElecOffDiagonalTensor& x, const ElecOffDiagonalTensor& y)
  {
    throw CException(" elec offdiag normalize is not defined!");
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::normalize(x[i],y[i]);
  }

  static void setMax(ElecOffDiagonalTensor& x, const ElecOffDiagonalTensor& y)
  {
    throw CException(" elec offdiag setMax is not defined!");
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::setMax(x[i],y[i]);
  }



 private:
  T _data;
  
};



template<class T, int N>
inline ostream& operator<<(ostream &os,
                           const ElecOffDiagonalTensor<T,N> &v)
{
  v.printFromC(os);
  return os;
}


template<class T, int N>
ElecOffDiagonalTensor<T,N>
operator+(const ElecOffDiagonalTensor<T,N>& a, const ElecOffDiagonalTensor<T,N>& b)
{
  return ElecOffDiagonalTensor<T,N>(a) += b;
}

template<class T, int N>
ElecOffDiagonalTensor<T,N>
operator-(const ElecOffDiagonalTensor<T,N>& a, const ElecOffDiagonalTensor<T,N>& b)
{
  return ElecOffDiagonalTensor<T,N>(a) -= b;
}

template<class T, int N>
ElecOffDiagonalTensor<T,N>
operator-(const ElecOffDiagonalTensor<T,N>& a)
{
  return -ElecOffDiagonalTensor<T,N>(a);
}

template<class T, int N>
ElecOffDiagonalTensor<T,N>
operator*(const ElecOffDiagonalTensor<T,N>& a, const ElecOffDiagonalTensor<T,N>& b)
{
  throw CException("no operator defined for diag * diag");
}

template<class T, int N>
ElecOffDiagonalTensor<T,N>
operator*(const T s, const ElecOffDiagonalTensor<T,N>& a)
{
  return ElecOffDiagonalTensor<T,N>(a) *= s;
}

template<class T, int N>
ElecOffDiagonalTensor<T,N>
operator*(const ElecOffDiagonalTensor<T,N>& a, const T s)
{
  return ElecOffDiagonalTensor<T,N>(a) *= s;
}

template<class T, int N>
Vector<T,N+1>
operator*(const ElecOffDiagonalTensor<T,N>& a, const Vector<T,N+1>& b)
{
  Vector<T,N+1> r;
  r = 0;
  r[N] += a[N] * b[N]; 
  
  return r;
}

template<class T, int N>
ElecOffDiagonalTensor<T,N>
operator/(const ElecOffDiagonalTensor<T,N>& a, const T s)
{
  return ElecOffDiagonalTensor<T,N>(a) /= s;
}

template<class T, int N>
Vector<T,N+1>
operator/(const Vector<T,N+1>& a, const ElecOffDiagonalTensor<T,N>& b) 
{
  Vector<T,N+1> x;
  x = 0;
  x[N] = a[N] / b[N];
  return x;
}


template<class T, int N>
ElecOffDiagonalTensor<T,N>
operator/(const ElecOffDiagonalTensor<T,N>& a, const ElecOffDiagonalTensor<T,N>& b) 
{
  throw CException("operator not defined for offdiag/offdiag"); 
}



template<class T, int N>
ElecOffDiagonalTensor<T,N>
operator/(const T s, const ElecOffDiagonalTensor<T,N>& a)
{
  throw CException("operator not defined for s/diag");
}
  
template<class T, int N>
T DiagToOffDiag(const ElecOffDiagonalTensor<T,N>& x)
{
  throw CException("diag to offdiag not defined");
}



#endif








