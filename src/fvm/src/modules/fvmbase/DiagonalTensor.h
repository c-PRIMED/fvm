#ifndef _DIAGONALTENSOR_H_
#define _DIAGONALTENSOR_H_

#include "NumType.h"
#include "Vector.h"
#include <sstream>


template <class T, int N>
class DiagonalTensor
{
public:
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef typename NumTypeTraits<T>::T_BuiltIn T_BuiltIn;

  static string getTypeName()
  {
    return "DiagonalTensor<" + NumTypeTraits<T>::getTypeName() +
      "," + intAsString(N) +
      ">";
  }
  static int getDimension() {return 1;}
  
  static void getShape(int *shp) { *shp = N;}
  static int getDataSize()
  {
    return N*NumTypeTraits<T>::getDataSize();
  }

  T& operator[](int n) {return _data[n];}
  const T& operator[](int n) const {return _data[n];}

  void printFromC(ostream &os) const
  {
    os << "[ " ;
    for(int i=0;i<N;i++)
      os << _data[i] << " " ;
    os << "]";
  }

  static void write(FILE* fp, const DiagonalTensor& x)
  {
    for(int i=0; i<N; i++)
    {
        NumTypeTraits<T>::write(fp,x[i]);
        fprintf(fp, " ");
    }
  }
  
  DiagonalTensor& operator=(const T& o)
  {
    for(int i=0;i<N;i++)
      _data[i] = o;
    return *this;
  }

  DiagonalTensor& operator=(const DiagonalTensor& o)
  {
    for(int i=0;i<N;i++)
      _data[i] = o[i];
    return *this;
  }

  DiagonalTensor operator-()
  {
    DiagonalTensor r;
    for(int i=0;i<N;i++)
      r[i]=-_data[i];
    return r;
  }

  DiagonalTensor& operator+=(const DiagonalTensor& o)
  {
    for(int i=0;i<N;i++)
      _data[i] += o[i];
    return *this;
  }

  DiagonalTensor& operator+=(const T s)
  {
    for(int i=0;i<N;i++)
      _data[i] += s;
    return *this;
  }

  DiagonalTensor& operator-=(const DiagonalTensor& o)
  {
    for(int i=0;i<N;i++)
      _data[i] -= o[i];
    return *this;
  }

  DiagonalTensor& operator-=(const T s)
  {
    for(int i=0;i<N;i++)
      _data[i] -= s;
    return *this;
  }

  DiagonalTensor& operator/=(const T s)
  {
    for(int i=0;i<N;i++)
      _data[i] /= s;
    return *this;
  }

  DiagonalTensor& operator/=(const DiagonalTensor& o)
  {
    for(int i=0;i<N;i++)
      _data[i] /= o[i];
    return *this;
  }

  DiagonalTensor& operator*=(const T s)
  {
    for(int i=0;i<N;i++)
      _data[i] *= s;
    return *this;
  }

  DiagonalTensor& operator*=(const DiagonalTensor& o)
  {
    for(int i=0;i<N;i++)
      _data[i] *= o[i];
    return *this;
  }

  void zero()
  {
    for(int i=0;i<N;i++) _data[i] = NumTypeTraits<T>::getZero();
  }
  
  static DiagonalTensor getZero()
  {
    DiagonalTensor z;
    z.zero();
    return z;
  }

  static double doubleMeasure(const DiagonalTensor& x)
  {
    double m=0;
    for (int i=0; i<N; i++)
      m += NumTypeTraits<T>::doubleMeasure(x[i]);
    return m;
  }

  static DiagonalTensor getNegativeUnity()
  {
    DiagonalTensor n;
    for(int i=0; i<N; i++)
      n._data[i] = NumTypeTraits<T>::getNegativeUnity();
    
    return n;
  }

  static DiagonalTensor getUnity()
  {
    DiagonalTensor n;
    for(int i=0; i<N; i++)
      n._data[i] = NumTypeTraits<T>::getUnity();
    
    return n;
  }

  static void accumulateOneNorm(DiagonalTensor& sum, const DiagonalTensor& v)
  {
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::accumulateOneNorm(sum[i],v[i]);
  }

  static void accumulateDotProduct(DiagonalTensor& sum, const DiagonalTensor& v0,
                                const DiagonalTensor& v1)
  {
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::accumulateDotProduct(sum[i],v0[i],v1[i]);
  }

private:
  T _data[N];
};

template<class T, int N>
inline ostream& operator<<(ostream &os,
                           const DiagonalTensor<T,N> &v)
{
  v.printFromC(os);
  return os;
}


template<class T, int N>
DiagonalTensor<T,N>
operator+(const DiagonalTensor<T,N>& a, const DiagonalTensor<T,N>& b)
{
  return DiagonalTensor<T,N>(a) += b;
}

template<class T, int N>
DiagonalTensor<T,N>
operator-(const DiagonalTensor<T,N>& a, const DiagonalTensor<T,N>& b)
{
  return DiagonalTensor<T,N>(a) -= b;
}

template<class T, int N>
DiagonalTensor<T,N>
operator-(const DiagonalTensor<T,N>& a)
{
  return -DiagonalTensor<T,N>(a);
}

template<class T, int N>
DiagonalTensor<T,N>
operator*(const DiagonalTensor<T,N>& a, const DiagonalTensor<T,N>& b)
{
  return DiagonalTensor<T,N>(a) *= b;
}

template<class T, int N>
DiagonalTensor<T,N>
operator*(const T s, const DiagonalTensor<T,N>& a)
{
  return DiagonalTensor<T,N>(a) *= s;
}

template<class T, int N>
DiagonalTensor<T,N>
operator*(const DiagonalTensor<T,N>& a, const T s)
{
  return DiagonalTensor<T,N>(a) *= s;
}

template<class T, int N>
Vector<T,N>
operator*(const DiagonalTensor<T,N>& a, const Vector<T,N>& b)
{
  Vector<T,N> r;
  for(int i=0; i<N; i++) r[i] = a[i]*b[i];
  return r;
}

template<class T, int N>
DiagonalTensor<T,N>
operator/(const DiagonalTensor<T,N>& a, const T s)
{
  return DiagonalTensor<T,N>(a) /= s;
}

template<class T, int N>
Vector<T,N>
operator/(const Vector<T,N>& a, const DiagonalTensor<T,N>& b) 
{
  Vector<T,N> r;
  for(int i=0; i<N; i++) r[i] = a[i]/b[i];
  return r;
}

template<class T, int N>
DiagonalTensor<T,N>
operator/(const T s, const DiagonalTensor<T,N>& a)
{
  DiagonalTensor<T,N> r;
  for(int i=0; i<N; i++) r[i] = s/a[i];
  return r;
}


#endif

