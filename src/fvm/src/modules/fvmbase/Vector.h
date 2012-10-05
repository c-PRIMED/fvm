// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _VECTOR_H_
#define _VECTOR_H_

#include "NumType.h"

#include <sstream>
inline string intAsString(const int i)
{
  ostringstream osBuffer;
  osBuffer << i;
  return osBuffer.str();
}

template <class T, int N>
class Vector
{
public:
  typedef Vector<T,N> This_T;
  
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef typename NumTypeTraits<T>::T_BuiltIn T_BuiltIn;
  typedef Vector<T_BuiltIn,N> VBuiltIn;

  typedef T T_NumType;
  
  enum {Length = N};
    
  // not initialized
  Vector() {}
  
  Vector(const Vector& o)
  {
    for(int i=0;i<N;i++)
      _data[i] = o[i];
  }

  static string getTypeName()
  {
    return "Vector<" + NumTypeTraits<T>::getTypeName() +
      "," + intAsString(N) +
      ">";
  }
  static int getDimension() {return NumTypeTraits<T>::getDimension()+1;}
  
  static void getShape(int *shp) { *shp = N; NumTypeTraits<T>::getShape(shp+1);}

  void* getData() 
  {
     return _data;
  }
  static int getDataSize()
  {
    return N*NumTypeTraits<T>::getDataSize();
  }

  T& operator[](int n) {return _data[n];}
  const T& operator[](int n) const {return _data[n];}

  void printFromC(ostream &os) const
  {
    os << "[ " ;
    os << scientific;
    for(int i=0;i<N;i++)
      os << _data[i] << " " ;
    os << "]";
  }

  static void write(FILE* fp, const Vector& x)
  {
    for(int i=0; i<N; i++)
    {
        NumTypeTraits<T>::write(fp,x[i]);
        fprintf(fp, " ");
    }
  }
  
  Vector& operator=(const T& o)
  {
    for(int i=0;i<N;i++)
      _data[i] = o;
    return *this;
  }

#if 0
  Vector& operator=(const int o)
  {
    for(int i=0;i<N;i++)
      _data[i] = o;
    return *this;
  }
#endif
  
  Vector& operator=(const Vector& o)
  {
    for(int i=0;i<N;i++)
      _data[i] = o[i];
    return *this;
  }

  Vector operator-()
  {
    Vector r;
    for(int i=0;i<N;i++)
      r[i]=-_data[i];
    return r;
  }

  Vector& operator+=(const Vector& o)
  {
    for(int i=0;i<N;i++)
      _data[i] += o[i];
    return *this;
  }

  Vector& operator-=(const Vector& o)
  {
    for(int i=0;i<N;i++)
      _data[i] -= o[i];
    return *this;
  }

  Vector& operator/=(const T s)
  {
    for(int i=0;i<N;i++)
      _data[i] /= s;
    return *this;
  }

  // elementwise operation
  Vector& operator/=(const Vector& o)
  {
    for(int i=0;i<N;i++)
      _data[i] /= o[i];
    return *this;
  }

  Vector& operator*=(const T s)
  {
    for(int i=0;i<N;i++)
      _data[i] *= s;
    return *this;
  }

  // elementwise operation
  Vector& operator*=(const Vector& o)
  {
    for(int i=0;i<N;i++)
      _data[i] *= o[i];
    return *this;
  }

  void zero()
  {
    for(int i=0;i<N;i++) _data[i] = NumTypeTraits<T>::getZero();
  }
  
  T mag2() const
  {
    T r(NumTypeTraits<T>::getZero());
    for(int i=0; i<N; i++)
      r+=_data[i]*_data[i];
    return r;
  }

  bool operator<(const double tolerance) const
  {
    return mag2() < tolerance*tolerance;
  }

  bool greaterThanElementwise(const Vector& y) const
  {
    for(int i=0; i<N; i++)
      if (_data[i] <= y[i])
        return false;
    return true;
  }
  
  static Vector getZero()
  {
    Vector z;
    z.zero();
    return z;
  }

  static void accumulateOneNorm(Vector& sum, const Vector& v)
  {
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::accumulateOneNorm(sum[i],v[i]);
  }

  static void accumulateDotProduct(Vector& sum, const Vector& v0, const Vector& v1)
  {
    for(int i=0; i<N; i++)
      sum[i]+=v0[i]*v1[i];
  }

  static void reduceSum(T_Scalar& sum, const This_T& x)
  {
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::reduceSum(sum,x[i]);
  }

  static void safeDivide(Vector& x, const Vector& y)
  {
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::safeDivide(x[i],y[i]);
  }

  static void normalize(Vector& x, const Vector& y)
  {
    const T ymag = sqrt(y.mag2());
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::safeDivide(x[i],ymag);
  }
 

  static void setMax(Vector& x, const Vector& y)
  {
    if (x.mag2() < y.mag2())
      for(int i=0; i<N; i++)
        x[i]=y[i];
  }
 
private:
  T _data[N];
};

template<class T, int N>
inline ostream& operator<<(ostream &os,
                           const Vector<T,N> &v)
{
  v.printFromC(os);
  return os;
}

template<class T>
Vector<T,3>
cross(const Vector<T,3>& a, const Vector<T,3>& b)
{
  Vector<T,3> c;
  c[0] = a[1]*b[2]-a[2]*b[1];
  c[1] = a[2]*b[0]-a[0]*b[2];
  c[2] = a[0]*b[1]-a[1]*b[0];
  return c;
}

template<class T>
T
dot(const Vector<T,3>& a, const Vector<T,3>& b)
{
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

template<class T>
T
mag(const Vector<T,3>& a)
{
  return sqrt(dot(a,a));
}

template<class T>
T
mag2(const Vector<T,3>& a)
{
  return dot(a,a);
}

template<class T, int N>
Vector<T,N>
operator+(const Vector<T,N>& a, const Vector<T,N>& b)
{
  return Vector<T,N>(a) += b;
}

template<class T, int N>
Vector<T,N>
operator-(const Vector<T,N>& a, const Vector<T,N>& b)
{
  return Vector<T,N>(a) -= b;
}

template<class T, int N>
Vector<T,N>
operator-(const Vector<T,N>& a)
{
  return -Vector<T,N>(a);
}

template<class T, int N>
Vector<T,N>
operator*(const T s, const Vector<T,N>& a)
{
  return Vector<T,N>(a) *= s;
}

template<class T, int N>
Vector<T,N>
operator*(const Vector<T,N>& a, const T s)
{
  return Vector<T,N>(a) *= s;
}

// does elemenwise operation
template<class T, int N>
Vector<T,N>
operator*(const Vector<T,N>& a, const Vector<T,N>& b)
{
  return Vector<T,N>(a) *= b;
}

template<class T, int N>
Vector<T,N>
operator/(const Vector<T,N>& a, const T s)
{
  return Vector<T,N>(a) /= s;
}

// does elementwise operation
template<class T, int N>
Vector<T,N>
operator/(const Vector<T,N>& a, const Vector<T,N>& b)
{
  return Vector<T,N>(a) /= b;
}

template<class T, int N>
struct ArrayScalarTraits<Vector<T,N> >
{
  static void limit(Vector<T,N>& val, const double min, const double max)
  {
    for(int i=0; i<N; i++)
    {
        if (val[i] < min)
          val[i] = min;
        else if (val[i] > max)
          val[i] = max;
    }
  }
};

#endif
