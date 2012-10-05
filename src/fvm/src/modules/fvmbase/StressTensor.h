// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _STRESSTENSOR_H_
#define _STRESSTENSOR_H_

#include "NumType.h"
#include "Vector.h"
#include <sstream>


/**
 * Symmetric 3x3 tensor to store stresses. Components are stored in the following order:
 * xx, yy, zz, xy, yz, zx
 *
 * Only used for exporting the stresses to MPM so most of the algebra
 * operations aren't currently implemented.
 * 
 */

template <class T>
class StressTensor
{
public:
  typedef StressTensor<T> This_T;
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef typename NumTypeTraits<T>::T_BuiltIn T_BuiltIn;

  static string getTypeName()
  {
    return "StressTensor<" + NumTypeTraits<T>::getTypeName() +
      ">";
  }

  // treated as one d for numpy export purposes
  static int getDimension() {return 1;}
  
  static void getShape(int *shp) { *shp = 6;}

  static int getDataSize()
  {
    return 6*NumTypeTraits<T>::getDataSize();
  }

  T& operator[](int n) {return _data[n];}
  const T& operator[](int n) const {return _data[n];}

  StressTensor& operator=(const T& o)
  {
    for(int i=0;i<6;i++)
      _data[i] = o;
    return *this;
  }

  StressTensor& operator=(const StressTensor& o)
  {
    for(int i=0;i<6;i++)
      _data[i] = o[i];
    return *this;
  }

  StressTensor operator-()
  {
    StressTensor r;
    for(int i=0;i<6;i++)
      r[i]=-_data[i];
    return r;
  }

  StressTensor& operator+=(const StressTensor& o)
  {
    for(int i=0;i<6;i++)
      _data[i] += o[i];
    return *this;
  }

  StressTensor& operator+=(const T s)
  {
    for(int i=0;i<6;i++)
      _data[i] += s;
    return *this;
  }

  StressTensor& operator-=(const StressTensor& o)
  {
    for(int i=0;i<6;i++)
      _data[i] -= o[i];
    return *this;
  }

  StressTensor& operator-=(const T s)
  {
    for(int i=0;i<6;i++)
      _data[i] -= s;
    return *this;
  }

  StressTensor& operator/=(const T s)
  {
    for(int i=0;i<6;i++)
      _data[i] /= s;
    return *this;
  }

  StressTensor& operator/=(const StressTensor& o)
  {
    throw;
  }

  StressTensor& operator*=(const T s)
  {
    for(int i=0;i<6;i++)
      _data[i] *= s;
    return *this;
  }

  StressTensor& operator*=(const StressTensor& o)
  {
    throw;
  }

  void zero()
  {
    for(int i=0;i<6;i++) _data[i] = NumTypeTraits<T>::getZero();
  }
  
  T mag2() const
  {
    throw;
  }

  bool operator<(const double tolerance) const
  {
    return mag2() < tolerance*tolerance;
  }

  static StressTensor getZero()
  {
    StressTensor z;
    z.zero();
    return z;
  }

  static double doubleMeasure(const StressTensor& x)
  {
    throw;
  }

  static StressTensor getNegativeUnity()
  {
    throw;
  }

  static StressTensor getUnity()
  {
    throw;
  }

  static void accumulateOneNorm(StressTensor& sum, const StressTensor& v)
  {
    throw;
  }

  static void accumulateDotProduct(StressTensor& sum, const StressTensor& v0,
                                const StressTensor& v1)
  {
    throw;
  }

  static void reduceSum(T_Scalar& sum, const This_T& x)
  {
    throw;
  }
  
  static void safeDivide(StressTensor& x, const StressTensor& y)
  {
    throw;
  }

  static void normalize(StressTensor& x, const StressTensor& y)
  {
    throw;
  }

  static void setMax(StressTensor& x, const StressTensor& y)
  {
    throw;
  }

private:
  T _data[6];
};

template<class T>
inline ostream& operator<<(ostream &os,
                           const StressTensor<T> &v)
{
  throw;
}

template<class T>
StressTensor<T>
operator+(const StressTensor<T>& a, const StressTensor<T>& b)
{
  return StressTensor<T>(a) += b;
}

template<class T>
StressTensor<T>
operator-(const StressTensor<T>& a, const StressTensor<T>& b)
{
  return StressTensor<T>(a) -= b;
}

template<class T>
StressTensor<T>
operator-(const StressTensor<T>& a)
{
  return -StressTensor<T>(a);
}

template<class T>
StressTensor<T>
operator*(const StressTensor<T>& a, const StressTensor<T>& b)
{
  throw;
}

template<class T>
StressTensor<T>
operator*(const T s, const StressTensor<T>& a)
{
  throw;
}

template<class T>
StressTensor<T>
operator*(const StressTensor<T>& a, const T s)
{
  throw;
}

template<class T>
StressTensor<T>
operator/(const StressTensor<T>& a, const T s)
{
  throw;
}


template<class T>
StressTensor<T>
operator/(const T s, const StressTensor<T>& a)
{
  throw;
}


#endif
