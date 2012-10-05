// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _GRADIENT_H_
#define _GRADIENT_H_


#include "Vector.h"

template <class T>
class Gradient
{
public:
  typedef Gradient<T> This_T;
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef typename NumTypeTraits<T>::T_BuiltIn T_BuiltIn;
  
  typedef Vector<T_Scalar,3> Coord;

  static string getTypeName()
  {
    return "Gradient<" + NumTypeTraits<T>::getTypeName() +      ">";
  }
  
  static int getDimension() {return NumTypeTraits<T>::getDimension()+1;}
  
  static void getShape(int *shp)
  {
    *shp = 3;
    NumTypeTraits<T>::getShape(shp+1);
  }

  static int getDataSize()
  {
    return 3*NumTypeTraits<T>::getDataSize();
  }

  T& operator[](int n) {return _data[n];}
  const T& operator[](int n) const {return _data[n];}


  Gradient& operator=(const T_Scalar& o)
  {
    for(int i=0;i<3;i++)
      _data[i] = o;
    return *this;
  }

  Gradient& operator=(const Gradient& o)
  {
    for(int i=0;i<3;i++)
      _data[i] = o._data[i];
    return *this;
  }

  void accumulate(const Coord& wt, const T& v)
  {
    for(int i=0;i<3;i++)
      _data[i] += wt[i]*v;
  }

  Gradient& operator+=(const Gradient& o)
  {
    for(int i=0;i<3;i++)
      _data[i] += o._data[i];
    return *this;
  }

  Gradient& operator-=(const Gradient& o)
  {
    for(int i=0;i<3;i++)
      _data[i] -= o._data[i];
    return *this;
  }

  Gradient operator-() const
  {
    Gradient r;
    for(int i=0;i<3;i++)
      r._data[i] = -_data[i];
    return r;
  }

  Gradient& operator*=(const T_Scalar& s)
  {
    for(int i=0;i<3;i++)
      _data[i] *= s;
    return *this;
  }

  Gradient& operator*=(const Gradient& o)
  {
    for(int i=0;i<3;i++)
      _data[i] *= o._data[i];
    return *this;
  }

  T operator*=(const Coord& v)
  {
    T r(NumTypeTraits<T>::getZero());
    for(int i=0;i<3;i++)
      r += _data[i]*v[i];
    return r;
  }

  Gradient& operator/=(const T_Scalar& s)
  {
    for(int i=0;i<3;i++)
      _data[i] /= s;
    return *this;
  }

  // elementwise
  Gradient& operator/=(const Gradient& o)
  {
    for(int i=0;i<3;i++)
      _data[i] /= o._data[i];
    return *this;
  }

  void zero()
  {
    for(int i=0;i<3;i++) _data[i] = NumTypeTraits<T>::getZero();
  }
  
  static void accumulateOneNorm(Gradient& sum, const Gradient& v)
  {
    for(int i=0; i<3; i++)
      NumTypeTraits<T>::accumulateOneNorm(sum._data[i],v[i]);
  }

  static void accumulateDotProduct(Gradient& sum, const Gradient& v0, const Gradient& v1)
  {
    for(int i=0; i<3; i++)
      NumTypeTraits<T>::accumulateDotProduct(sum._data[i],v0._data[i],v1._data[i]);
  }

  static void reduceSum(T_Scalar& sum, const This_T& x)
  {
    for(int i=0; i<3; i++)
      NumTypeTraits<T>::reduceSum(sum,x[i]);
  }
  
  static void safeDivide(Gradient& x, const Gradient& y)
  {
    for(int i=0; i<3; i++)
      NumTypeTraits<T>::safeDivide(x._data[i],y._data[i]);
  }

  static void normalize(Gradient& x, const Gradient& y)
  {
    for(int i=0; i<3; i++)
      NumTypeTraits<T>::normalize(x._data[i],y._data[i]);
  }

  static void setMax(Gradient& x, const Gradient& y)
  {
    for(int i=0; i<3; i++)
      NumTypeTraits<T>::setMax(x._data[i],y._data[i]);
  }
  static Gradient getZero()
  {
    Gradient z;
    z.zero();
    return z;
  }

  void print(ostream &os) const
  {
    os << "[ " ;
    for(int i=0;i<3;i++)
      os << _data[i] << " " ;
    os << "]";
  }

  static void write(FILE* fp, const Gradient& x)
  {
    for(int i=0; i<3; i++)
    {
        NumTypeTraits<T>::write(fp,x._data[i]);
        fprintf(fp, " ");
    }
  }

  T mag2() const
  {
    T r(NumTypeTraits<T>::getZero());
    for(int i=0; i<3; i++)
      r+=_data[i]*_data[i];
    return r;
  }

  bool operator<(const double tolerance) const
  {
    return mag2() < tolerance*tolerance;
  }
private:
  T _data[3];
};

template<class T>
inline ostream& operator<<(ostream &os,
                           const Gradient<T> &v)
{
  v.print(os);
  return os;
}

template<class T>
Gradient<T>
operator+(const Gradient<T>& a, const Gradient<T>& b)
{
  return Gradient<T>(a) += b;
}

template<class T>
Gradient<T>
operator*(const Gradient<T>& a, const typename NumTypeTraits<T>::T_Scalar& s)
{
  return Gradient<T>(a) *= s;
}

template<class T>
T
operator*(const Gradient<T>& a, const typename Gradient<T>::Coord& v)
{
  return Gradient<T>(a) *= v;
}

// does elemenwise operation
template<class T>
Gradient<T>
operator*(const Gradient<T>& a, const Gradient<T>& b)
{
  return Gradient<T>(a) *= b;
}

template<class T>
Gradient<T>
operator/(const Gradient<T>& a, const typename NumTypeTraits<T>::T_Scalar& s)
{
  return Gradient<T>(a) /= s;
}

// for statements like x /= 2;
template<class T>
Gradient<T>
operator/(const Gradient<T>& a, const int& s)
{
  return Gradient<T>(a) /= typename NumTypeTraits<T>::T_Scalar(s);
}

template<class T>
Gradient<T>
operator-(const Gradient<T>& a)
{
  return -Gradient<T>(a);
}


#endif

