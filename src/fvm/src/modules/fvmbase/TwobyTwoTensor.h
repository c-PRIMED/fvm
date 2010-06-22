#ifndef _TWOBYTWOTENSOR_H_
#define _TWOBYTWOTENSOR_H_

#include "NumType.h"
#include "Vector.h"
#include <sstream>

const int  N = 4;

template <class T>
class TwobyTwoTensor
{
 public:
  typedef TwobyTwoTensor<T> This_T;
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef typename NumTypeTraits<T>::T_BuiltIn T_BuiltIn;
  //static const int  N = 4;

  TwobyTwoTensor() {}

  TwobyTwoTensor(const TwobyTwoTensor& o){
    for (int i=0; i<N; i++)
       _data[i] = o._data[i];
  }

  TwobyTwoTensor(const T& s){
    _data[0] = s;
    _data[3] = s;
  }

  static string getTypeName()
  {
    return "TwobyTwoTensor<" + NumTypeTraits<T>::getTypeName() +
      "," + intAsString(N) +
      ">";
  }
  
  //what is dimension?
  static int getDimension() {return 1;}

  //what is shape?
  static void getShape(int *shp) { *shp = 4;}

  static int getDataSize()
  {
    return N*NumTypeTraits<T>::getDataSize();
  }

  T& operator[] (int n)  { return _data[n]; }

  const T& operator[](int n) const {return _data[n];}

  void printFromC(ostream &os) const
  {
    os << "[ " << _data[0] << ", " << _data[1] << endl;
    os << "  " << _data[2] << ", " << _data[3] << " ]" << endl;    
  }

  static void write(FILE* fp, const TwobyTwoTensor& x)
  {
    for(int i=0; i<N; i++)
    {
        NumTypeTraits<T>::write(fp,x[i]);
        fprintf(fp, " ");
    }
  }

  TwobyTwoTensor& operator=(const T& o)
  {
    //purpose here is to setup the coeffcient of diag
    _data[0] = o;
    _data[3] = o;
    return *this;
  }
  
  TwobyTwoTensor& operator=(const TwobyTwoTensor& o)
  {
    for(int i=0;i<N;i++)
      _data[i] = o[i];
    return *this;
  }

  TwobyTwoTensor operator-()
  {
    TwobyTwoTensor r;
    for(int i=0;i<N;i++)
      r[i]=-_data[i];
    return r;
  }

  TwobyTwoTensor& operator+=(const TwobyTwoTensor& o)
  {
    for(int i=0;i<N;i++)
      _data[i] += o[i];
    return *this;
  }

  TwobyTwoTensor& operator+=(const T o)
  {
    _data[0] += o;
    _data[3] += o;
    return *this;
  }
  
  TwobyTwoTensor& operator-=(const TwobyTwoTensor& o)
  {
    for(int i=0;i<N;i++)
      _data[i] -= o[i];
    return *this;
  }

  TwobyTwoTensor& operator-=(const T o)
  {
    //for(int i=0;i<N;i++)
    //  _data[i] -= o;
    _data[0] -= o;
    _data[3] -= o;
    return *this;
  }
  
  TwobyTwoTensor& operator*=(const T o)
  {
    for(int i=0;i<N;i++)
      _data[i] *= o;
    return *this;
  }

  TwobyTwoTensor& operator*=(const TwobyTwoTensor& o)
  {
    TwobyTwoTensor r;
    r[0] = _data[0] * o[0] + _data[1] * o[2];
    r[1] = _data[0] * o[1] + _data[1] * o[3];
    r[2] = _data[2] * o[0] + _data[3] * o[2];
    r[3] = _data[2] * o[1] + _data[3] * o[3];
    return *this = r;
  }
  
  TwobyTwoTensor& inverse()
  {
    TwobyTwoTensor& o = *this;
    TwobyTwoTensor r;
    T det = o[0]*o[3] - o[1]*o[2];
    T detInv = T(1.0) / det;
    r[0] = o[3] * detInv;
    r[1] = -o[1] * detInv;
    r[2] = -o[2] * detInv;
    r[3] = o[0] * detInv;
    return *this = r;
  }

  const TwobyTwoTensor inverse(const TwobyTwoTensor& o)
  {
    TwobyTwoTensor r;
    T det = o[0]*o[3] - o[1]*o[2];
    T detInv = T(1.0) / det;
    r[0] = o[3] * detInv;
    r[1] = -o[1] * detInv;
    r[2] = -o[2] * detInv;
    r[3] = o[0] * detInv;
    return r;
  }

  TwobyTwoTensor& operator/=(const T o)
  {
    for(int i=0;i<N;i++)
      _data[i] /= o;
    return *this;
  }

  TwobyTwoTensor& operator/=(const TwobyTwoTensor& o)
  {
    return (*this) *= inverse(o); 
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

  static TwobyTwoTensor getZero()
  {
    TwobyTwoTensor z;
    z.zero();
    return z;
  }

  static double doubleMeasure(const TwobyTwoTensor& x)
  {
    double m=0;
    for (int i=0; i<N; i++)
      m += NumTypeTraits<T>::doubleMeasure(x[i]);
    return m;
  }

  static TwobyTwoTensor getNegativeUnity()
  {
    TwobyTwoTensor x;
    for(int i=0; i<N; i++)
      x._data[i] = NumTypeTraits<T>::getNegativeUnity();    
    return x;
  }

  static TwobyTwoTensor getUnity()
  {
    TwobyTwoTensor x;
    for(int i=0; i<N; i++)
      x._data[i] = NumTypeTraits<T>::getUnity();
    return x;
  }  

  static void accumulateOneNorm(TwobyTwoTensor& sum, const TwobyTwoTensor& v)
  {
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::accumulateOneNorm(sum[i],v[i]);
  }

  static void accumulateDotProduct(TwobyTwoTensor& sum, const TwobyTwoTensor& v0,
                                const TwobyTwoTensor& v1)
  {
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::accumulateDotProduct(sum[i],v0[i],v1[i]);
  }

  static void reduceSum(T_Scalar& sum, const This_T& x)
  {
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::reduceSum(sum,x[i]);
  }
  
  static void safeDivide(TwobyTwoTensor& x, const TwobyTwoTensor& y)
  {
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::safeDivide(x[i],y[i]);
  }

  static void normalize(TwobyTwoTensor& x, const TwobyTwoTensor& y)
  {
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::normalize(x[i],y[i]);
  }

  static void setMax(TwobyTwoTensor& x, const TwobyTwoTensor& y)
  {
    for(int i=0; i<N; i++)
      NumTypeTraits<T>::setMax(x[i],y[i]);
  }


 private:
  T  _data[4]; 

};

template<class T>
inline ostream& operator<<(ostream &os,
                           const TwobyTwoTensor<T> &v)
{
  v.printFromC(os);
  return os;
}


template<class T>
TwobyTwoTensor<T>
operator+(const TwobyTwoTensor<T>& a, const TwobyTwoTensor<T>& b)
{
  return TwobyTwoTensor<T>(a) += b;
}

template<class T>
TwobyTwoTensor<T>
operator-(const TwobyTwoTensor<T>& a, const TwobyTwoTensor<T>& b)
{
  return TwobyTwoTensor<T>(a) -= b;
}

template<class T>
TwobyTwoTensor<T>
operator-(const TwobyTwoTensor<T>& a)
{
  return -TwobyTwoTensor<T>(a);
}

template<class T>
TwobyTwoTensor<T>
operator*(const TwobyTwoTensor<T>& a, const TwobyTwoTensor<T>& b)
{
  return TwobyTwoTensor<T>(a) *= b;
}

template<class T>
TwobyTwoTensor<T>
operator*(const T s, const TwobyTwoTensor<T>& a)
{
  return TwobyTwoTensor<T>(a) *= s;
}

template<class T>
TwobyTwoTensor<T>
operator*(const TwobyTwoTensor<T>& a, const T s)
{
  return TwobyTwoTensor<T>(a) *= s;
}

template<class T>
Vector<T, 2>
operator*(const TwobyTwoTensor<T>& a, const Vector<T, 2>& b)
{
  Vector<T, 2> r;
  r[0] = a[0]*b[0] + a[1]*b[1];
  r[1] = a[2]*b[0] + a[3]*b[1];
  return r;
}

template<class T>
Vector<T, 2>
operator*(const Vector<T, 2>& a, const TwobyTwoTensor<T>& b )
{
  Vector<T, 2> r;
  r[0] = a[0]*b[0] + a[1]*b[2];
  r[1] = a[0]*b[1] + a[1]*b[3];
  return r;
}


template<class T>
TwobyTwoTensor<T>
operator/(const TwobyTwoTensor<T>& a, const T s)
{
  return TwobyTwoTensor<T>(a) /= s;
}

template<class T>
TwobyTwoTensor<T>
operator/(const TwobyTwoTensor<T>& a, const TwobyTwoTensor<T>& s)
{
  return TwobyTwoTensor<T>(a) /= s;
}

template<class T>
Vector<T,2>
operator/(const Vector<T,2>& a, const TwobyTwoTensor<T>& b) 
{
  Vector<T,2> r;
  TwobyTwoTensor<T> bInv = b;
  r = a * bInv.inverse();
  return r;
}


template<class T>
TwobyTwoTensor<T>
operator/(const T s, const TwobyTwoTensor<T>& a)
{
  TwobyTwoTensor<T> r;
  for(int i=0; i<N; i++) r[i] = s/a[i];
  return r;
}


#endif
