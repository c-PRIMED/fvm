// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _VECTORTRANSPOSE_H_
#define _VECTORTRANSPOSE_H_

#include "NumType.h"

#include "Vector.h"
#include "SquareTensor.h"

template<class T, int N>
class VectorTranspose
{
public:
  typedef VectorTranspose<T,N> This_T;
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef typename NumTypeTraits<T>::T_BuiltIn T_BuiltIn;
  typedef T T_NumType;
  enum {Length = N};

  static int NumPy_TypeNum;

  VectorTranspose()
  {}
  
  VectorTranspose(const VectorTranspose& x) :
    _v(x._v)
  {}
  
  VectorTranspose(const Vector<T,N>& v) :
    _v(v)
  {}

  static string getTypeName()
  {
    return "VectorTranspose<" + NumTypeTraits<T>::getTypeName() +
      "," + intAsString(N) +
      ">";
  }
  static int getDimension() {return NumTypeTraits<T>::getDimension()+1;}
  
  static void getShape(int *shp) { *shp = N; NumTypeTraits<T>::getShape(shp+1);}

  static int getDataSize()
  {
    return N*NumTypeTraits<T>::getDataSize();
  }

  T& operator[](int n) {return _v[n];}
  const T& operator[](int n) const {return _v[n];}

  VectorTranspose& operator=(const T& o)
  {
    _v = o;
    return *this;
  }

  VectorTranspose& operator=(const int o)
  {
    _v = o;
    return *this;
  }

  
  VectorTranspose& operator=(const VectorTranspose& o)
  {
    _v = o._v;
    return *this;
  }

  VectorTranspose operator-()
  {
    return VectorTranspose(-_v);
  }

  VectorTranspose& operator+=(const VectorTranspose& o)
  {
    _v += o._v;
    return *this;
  }

  VectorTranspose& operator-=(const VectorTranspose& o)
  {
    _v -= o._v;
    return *this;
  }

  VectorTranspose& operator/=(const T s)
  {
    _v /= s;
    return *this;
  }

  // elementwise operation
  VectorTranspose& operator/=(const VectorTranspose& o)
  {
    _v /= o._v;
    return *this;
  }

  VectorTranspose& operator*=(const T s)
  {
    _v *= s;
    return *this;
  }

  T operator*=(const Vector<T,N>& o)
  {
    return dot(_v,o);
  }


  // elementwise operation
  VectorTranspose& operator*=(const VectorTranspose& o)
  {
    _v *= o._v;
    return *this;
  }

  SquareTensor<T,N> getTensor(const Vector<T,N>& o)
  {
      SquareTensor<T,N> r;
      r.zero();
      for(int i=0;i<N;i++)
	for(int j=0;j<N;j++)
	  r(i,j)=o[i]*_v[j];
      return r;
  }

  void zero()
  {
    _v.zero();
  }
  
  static VectorTranspose getZero()
  {
    VectorTranspose z;
    z.zero();
    return z;
  }

  bool operator<(const double tolerance) const
  {
    return _v < tolerance;
  }

  void printFromC(ostream &os) const
  {
    _v.printFromC(os);
    os << "^T";
  }

  static void accumulateOneNorm(VectorTranspose& sum, const VectorTranspose& v)
  {
    Vector<T,N>::accumulateOneNorm(sum._v,v._v);
  }

  static void accumulateDotProduct(VectorTranspose& sum,
                                   const VectorTranspose& v0,
                                   const VectorTranspose& v1)
  {
    Vector<T,N>::accumulateDotProduct(sum._v, v0._v, v1._v);
  }

  static void safeDivide(VectorTranspose& x, const VectorTranspose& y)
  {
    Vector<T,N>::safeDivide(x._v,y._v);
  }

  static void normalize(VectorTranspose& x, const VectorTranspose& y)
  {
    Vector<T,N>::normalize(x._v,y._v);
  }

  static void setMax(VectorTranspose& x, const VectorTranspose& y)
  {
    Vector<T,N>::setMax(x._v,y._v);
  }

  static void reduceSum(T_Scalar& sum, const This_T& x)
  {
    Vector<T,N>::reduceSum(sum,x._v);
  }

private:
  Vector<T,N> _v;
};

template<class T, int N>
int VectorTranspose<T,N>::NumPy_TypeNum = -1;

template<class T, int N>
VectorTranspose<T,N>
operator+(const VectorTranspose<T,N>& a, const VectorTranspose<T,N>& b)
{
  return VectorTranspose<T,N>(a) += b;
}

template<class T, int N>
VectorTranspose<T,N>
operator-(const VectorTranspose<T,N>& a, const VectorTranspose<T,N>& b)
{
  return VectorTranspose<T,N>(a) -= b;
}

template<class T, int N>
VectorTranspose<T,N>
operator-(const VectorTranspose<T,N>& a)
{
  return -VectorTranspose<T,N>(a);
}

template<class T, int N>
VectorTranspose<T,N>
operator*(const T s, const VectorTranspose<T,N>& a)
{
  return VectorTranspose<T,N>(a) *= s;
}

template<class T, int N>
VectorTranspose<T,N>
operator*(const VectorTranspose<T,N>& a, const T s)
{
  return VectorTranspose<T,N>(a) *= s;
}

template<class T, int N>
T
operator*(const VectorTranspose<T,N>& a, const Vector<T,N>& b)
{
  return VectorTranspose<T,N>(a) *= b;
}

// does elemenwise operation
template<class T, int N>
VectorTranspose<T,N>
operator*(const VectorTranspose<T,N>& a, const VectorTranspose<T,N>& b)
{
  return VectorTranspose<T,N>(a) *= b;
}

template<class T, int N>
VectorTranspose<T,N>
operator/(const VectorTranspose<T,N>& a, const T s)
{
  return VectorTranspose<T,N>(a) /= s;
}

// does elemenwise operation
template<class T, int N>
VectorTranspose<T,N>
operator/(const VectorTranspose<T,N>& a, const VectorTranspose<T,N>& b)
{
  return VectorTranspose<T,N>(a) /= b;
}

template<class T, int N>
inline ostream& operator<<(ostream &os,
                           const VectorTranspose<T,N> &v)
{
  v.printFromC(os);
  return os;
}

#endif
