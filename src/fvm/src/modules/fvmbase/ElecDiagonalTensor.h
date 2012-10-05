// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _ELECDIAGONALTENSOR_H_
#define _ELECDIAGONALTENSOR_H_

#include "NumType.h"
#include "Vector.h"
#include <sstream>
#include "ElecOffDiagonalTensor.h"

//this structure is used to store the dianogal tensor in dielectric charging model
//the full diag tensor looks like
/*
  | Nt0                        Nt0c   |
  |      nt1                   Nt1c   |
  |           ...              ...    |
  |                Nt(n-1)    Nt(n-1)c |
  | Nct0 Nct1 ...  Nct(n-1)   Nc     |
*/
// here Nt is charges in trap and Nc is charges in conduction band
// if the level of trap depth is n, then there are 3n+1 non zeros
// it is stored in a 1D array data[3n+1]
// here is how the indices work
// (0, n-1)  refers to diagonal elements Nt0, ... Nt(n-1)
// (n, 2n-1) refers to upper diag elements Nt0c, .... Nt(n-1)c
// (2n, 3n-1) refers to lower diag elements
// (3n)  refers to Nc

// the offdiag would be just a scalar, since the only connection between
// neighbors is due to nc drift

template <class T, int N>
class ElecDiagonalTensor
{
public:
  enum { TN = 3*N+1 } ;
  typedef ElecDiagonalTensor<T,N> This_T;
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef typename NumTypeTraits<T>::T_BuiltIn T_BuiltIn;
  
  ElecDiagonalTensor()
  {}

  ElecDiagonalTensor(const ElecDiagonalTensor& o)
  {
    for(int i=0; i<TN; i++)
      _data[i] = o._data[i];
  }

  ElecDiagonalTensor(const T& o)
  {
    for(int i=0; i<N; i++)
      _data[i] = o;
    for(int i=N; i<TN; i++)
      _data[i] = 0;
  }

  static string getTypeName()
  {
    return "ElecDiagonalTensor<" + NumTypeTraits<T>::getTypeName() +
      "," + intAsString(N) +
      ">";
  }

  static int getDimension() {return 1;}

  static void getShape(int *shp) { *shp = TN;}

  static int getDataSize()
  {
    return TN*NumTypeTraits<T>::getDataSize();
  }

  T& operator[](int n) {return _data[n];}

  const T& operator[](int n) const {return _data[n];}

  void printFromC(ostream &os) const
  {
    os << "[ " ;
    for(int i=0;i<TN;i++)
      os << _data[i] << " " ;
    os << "]";
  }

  static void write(FILE* fp, const ElecDiagonalTensor& x)
  {
    for(int i=0; i<TN; i++)
    {
        NumTypeTraits<T>::write(fp,x[i]);
        fprintf(fp, " ");
    }
  }
  
  ElecDiagonalTensor& operator=(const T& o)
  {
    for(int i=0;i<N;i++)
      _data[i] = o;
   for(int i=N; i<TN; i++)
      _data[i] = 0; 
    return *this;
  }

  ElecDiagonalTensor& operator=(const ElecDiagonalTensor& o)
  {
    for(int i=0;i<TN;i++)
      _data[i] = o[i];
    return *this;
  }
  
  ElecDiagonalTensor& operator=(const ElecOffDiagonalTensor<T,N>& o)
  {
    for(int i=0;i<TN;i++)
      _data[i] = 0;
    _data[3*N] = o[N];
    return *this;
  }
  
  ElecDiagonalTensor operator-()
  {
    ElecDiagonalTensor r;
    for(int i=0;i<TN;i++)
      r[i]=-_data[i];
    return r;
  }

  ElecDiagonalTensor& operator+=(const ElecDiagonalTensor& o)
  {
    for(int i=0;i<TN;i++)
      _data[i] += o[i];
    return *this;
  }

  
  ElecDiagonalTensor& operator+=(const ElecOffDiagonalTensor<T,N> & o)
  {
    _data[3*N] += o[N];
    return *this;
  }
  ElecDiagonalTensor& operator+=(const T s)
  {
    for(int i=0;i<N;i++)
      _data[i] += s;
    _data[3*N] += s;
    return *this;
  }

  ElecDiagonalTensor& operator-=(const ElecDiagonalTensor& o)
  {
    for(int i=0;i<TN;i++)
      _data[i] -= o[i];
    return *this;
  }

  ElecDiagonalTensor& operator-=(const T s)
  {
    for(int i=0;i<N;i++)
      _data[i] -= s;
    _data[3*N] -= s;
    return *this;
  }

  ElecDiagonalTensor& operator/=(const T s)
  {
    throw CException("operator not defined for diag /= s");
    for(int i=0;i<TN;i++)
      _data[i] /= s;
    return *this;
  }

 
  ElecDiagonalTensor& operator/=(const ElecDiagonalTensor& o)
  {
    throw CException("no operator defined for /= diag");
  }
  

  ElecDiagonalTensor& operator*=(const T s)
  {
    throw CException("operator not defined for diag *= s");
    for(int i=0;i<TN;i++)
      _data[i] *= s;
    return *this;
  }

  
  ElecDiagonalTensor& operator*=(const ElecDiagonalTensor& o)
  {
    throw CException("no operator defined for *= diag");
  }
  
  void zero()
  {
    for(int i=0;i<TN;i++) _data[i] = NumTypeTraits<T>::getZero();
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

  static ElecDiagonalTensor getZero()
  {
    ElecDiagonalTensor z;
    z.zero();
    return z;
  }

  static double doubleMeasure(const ElecDiagonalTensor& x)
  {
      return 0.0;
  }

  static ElecDiagonalTensor getNegativeUnity()
  {
    ElecDiagonalTensor n;
    for(int i=0; i<N; i++)
      n._data[i] = NumTypeTraits<T>::getNegativeUnity();
    
    return n;
  }

  static ElecDiagonalTensor getUnity()
  {
    ElecDiagonalTensor n;
    for(int i=0; i<N; i++)
      n._data[i] = NumTypeTraits<T>::getUnity();
    
    return n;
  }

  static void accumulateOneNorm(ElecDiagonalTensor& sum, const ElecDiagonalTensor& v)
  {
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::accumulateOneNorm(sum[i],v[i]);
  }

  static void accumulateDotProduct(ElecDiagonalTensor& sum, const ElecDiagonalTensor& v0,
                                const ElecDiagonalTensor& v1)
  {
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::accumulateDotProduct(sum[i],v0[i],v1[i]);
  }

  static void reduceSum(T_Scalar& sum, const This_T& x)
  {
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::reduceSum(sum,x[i]);
  }
  
  static void safeDivide(ElecDiagonalTensor& x, const ElecDiagonalTensor& y)
  {
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::safeDivide(x[i],y[i]);
  }

  static void normalize(ElecDiagonalTensor& x, const ElecDiagonalTensor& y)
  {
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::normalize(x[i],y[i]);
  }

  static void setMax(ElecDiagonalTensor& x, const ElecDiagonalTensor& y)
  {
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::setMax(x[i],y[i]);
  }



 private:
  T _data[3*N+1];
  
};



template<class T, int N>
inline ostream& operator<<(ostream &os,
                           const ElecDiagonalTensor<T,N> &v)
{
  v.printFromC(os);
  return os;
}


template<class T, int N>
ElecDiagonalTensor<T,N>
operator+(const ElecDiagonalTensor<T,N>& a, const ElecDiagonalTensor<T,N>& b)
{
  return ElecDiagonalTensor<T,N>(a) += b;
}

template<class T, int N>
ElecDiagonalTensor<T,N>
operator-(const ElecDiagonalTensor<T,N>& a, const ElecDiagonalTensor<T,N>& b)
{
  return ElecDiagonalTensor<T,N>(a) -= b;
}

template<class T, int N>
ElecDiagonalTensor<T,N>
operator-(const ElecDiagonalTensor<T,N>& a)
{
  return -ElecDiagonalTensor<T,N>(a);
}

template<class T, int N>
ElecDiagonalTensor<T,N>
operator*(const ElecDiagonalTensor<T,N>& a, const ElecDiagonalTensor<T,N>& b)
{
  throw CException("no operator defined for diag * diag");
}

template<class T, int N>
ElecDiagonalTensor<T,N>
operator*(const T s, const ElecDiagonalTensor<T,N>& a)
{
  return ElecDiagonalTensor<T,N>(a) *= s;
}

template<class T, int N>
ElecDiagonalTensor<T,N>
operator*(const ElecDiagonalTensor<T,N>& a, const T s)
{
  return ElecDiagonalTensor<T,N>(a) *= s;
}

template<class T, int N>
Vector<T,N+1>
operator*(const ElecDiagonalTensor<T,N>& a, const Vector<T,N+1>& b)
{
  Vector<T,N+1> r;
  r = 0;
  for(int i=0; i<N; i++){
    r[i] = a[i]*b[i]+a[N+i]*b[N];
  }
  for(int i=0; i<=N; i++)
    r[N] += a[2*N+i]*b[i]; 
  
  return r;
}

template<class T, int N>
ElecDiagonalTensor<T,N>
operator/(const ElecDiagonalTensor<T,N>& a, const T s)
{
  return ElecDiagonalTensor<T,N>(a) /= s;
}

template<class T, int N>
Vector<T,N+1>
operator/(const Vector<T,N+1>& a, const ElecDiagonalTensor<T,N>& b) 
{
  //solve bx = a and return x
  //rename Mx = v where M=b v=a
  Vector<T,N+1> x;
  Vector<T,N+1> v;
  ElecDiagonalTensor<T,N> M;
  x = 0;
  M = b;
  v = a;
  
  //first step: eleminate M to upper diagonal  
  for(int i=0; i<N;i++){
    if(M[2*N+i] !=0.0 && M[i] != 0.0){
      M[3*N] -= M[N+i] * (M[2*N+i]/M[i]);
      v[N] -= v[i] * (M[2*N+i]/M[i]);
    }
  }
  
  //second step: go backward to solve x
  if (M[3*N] != 0.0)
    x[N] = v[N] / M[3*N];
  for(int i=0; i<N; i++){
    if (M[i] != 0.0)
      x[i] = (v[i] - M[N+i]*x[N]) / M[i];
  }
    
  return x;
}


template<class T, int N>
ElecDiagonalTensor<T,N>
operator/(const ElecDiagonalTensor<T,N>& a, const ElecDiagonalTensor<T,N>& b) 
{
  throw CException("operator not defined for diag/diag"); 
}
  
template<class T, int N>
ElecDiagonalTensor<T,N>
operator/(const T s, const ElecDiagonalTensor<T,N>& a)
{
  throw CException("operator not defined for s/diag");
}
  
template<class T, int N>
T DiagToOffDiag(const ElecDiagonalTensor<T,N>& x)
{
  throw CException("diag to offdiag not defined");
}


#endif








