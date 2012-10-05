// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _PC_H_
#define _PC_H_


#include "NumType.h"

#include "misc.h"
#include "PCSet.h"

// helper class for convenient initialization of PolynomialChaos type
// variables. This class, along with the overloaded comma (,) operator
// in PC allows us to write statements like
// a = 1, 2, 3;

template<class X, class T_iterator>
class Initializer
{
public:

  Initializer(T_iterator iter):
    _iter(iter)
  {}
  
  Initializer operator,(const X& x)
  {
    *_iter = x;
    return Initializer(_iter+1);
  }

protected:
  T_iterator _iter;
};


// used for compile time determination of the total size of the PC
// class's data member

template <int N>
struct Factorial 
{
  enum { value = N * Factorial<N - 1>::value };
};
 
template <>
struct Factorial<0> 
{
  enum { value = 1 };
};
 

class PCSet;

PCSet* createPCSet(const int order, const int dim);


template<int ORDER, int DIM>
class PC 
{
public:

  enum {N = Factorial<ORDER+DIM>::value/(Factorial<ORDER>::value*Factorial<DIM>::value) };
  
  typedef PC This_T;
  typedef PC T_Scalar;
  typedef NumTypeTraits<double>::T_BuiltIn T_BuiltIn;

  typedef Initializer<double,double*> PCInitializer;
  
  static string getTypeName()
  {
    return "PC";
  }
  
  PC()
  {
    for (int i=0;i<N;i++)
      _data[i]=0;
  }

  
  explicit PC(const double v)
  {
    _data[0] = v;
    for (int i=1;i<N;i++)
      _data[i]=0;
  }


  PC(const PC& o)
  {
    for (int i=0;i<N;i++)
      _data[i]=o._data[i];
  }
  
  ~PC()
  {}

  double& operator[](int n) {return _data[n];}
  const double& operator[](int n) const {return _data[n];}

  PCInitializer operator,(const double& x)
  {
    _data[1] = x;
    return PCInitializer(_data+2);
  }


  PC& operator=(const PC& o)
  {
    if (this == &o)
      return *this;
    for (int i=0;i<N;i++)
      _data[i]=o._data[i];
    return *this;
  }

  PC& operator=(const double& f)
  {
    _data[0] = f;
    for (int i=1;i<N;i++)
      _data[i]=0;
    return *this;
  }

  PC& operator+=(const PC& o)
  {
    for (int i=0;i<N;i++)
      _data[i] += o._data[i];
    return *this;
  }

  PC& operator-=(const PC& o)
  {
    for (int i=0;i<N;i++)
      _data[i] -= o._data[i];
    return *this;
  }

  PC& operator*=(const PC& o)
  {
    double pdata[N];
    _pcSet->Prod(_data,o._data,pdata);
    for(int i=0; i<N; i++)
      _data[i] = pdata[i];
    return *this;
  }


  PC& operator/=(const PC& o)
  {
    double pdata[N];
    _pcSet->Div(_data,o._data,pdata);
    for(int i=0; i<N; i++)
      _data[i] = pdata[i];
    return *this;
  }




  PC& operator+=(const double& o)
  {
    _data[0] += o;
    return *this;
  }

  PC& operator-=(const double& o)
  {
    _data[0] -= o;
    return *this;
  }

  PC& operator*=(const double& o)
  {
    for(int i=0; i<N; i++)
      _data[i] *= o;
    return *this;
  }

  PC& operator/=(const double& o)
  {
    for(int i=0; i<N; i++)
      _data[i] /= o;
    return *this;
  }

#define PC_RELATIONAL_OPS(opname,_op_)          \
  bool opname(const PC& o) const                \
  {                                             \
    return (_data[0] _op_ o._data[0]);          \
  }                                             \
    bool opname(const double& o) const          \
  {                                             \
    return (_data[0] _op_ o);                   \
  }                                             

  PC_RELATIONAL_OPS(operator>,>);
  PC_RELATIONAL_OPS(operator>=,>=);
  PC_RELATIONAL_OPS(operator<,<);
  PC_RELATIONAL_OPS(operator<=,<=);
  PC_RELATIONAL_OPS(operator==,==);
  PC_RELATIONAL_OPS(operator!=,!=);
  
#undef PC_RELATIONAL_OPS
  
  void print(ostream &os) const
  {
    os << "[ ";
    for(int i=0; i<N; i++)
      os << _data[i] << " ";
    os << "]";
  }

  static PC getZero()
  {
    return PC();
  }

  static PC getUnity()
  {
    PC one;
    one._data[0] = 1.0;
    return one;

  }

  static PC getNegativeUnity()
  {
    PC m_one;
    m_one._data[0] = -1.0;
    return m_one;
  }

  static double doubleMeasure(const PC& x)
  {
    return NumTypeTraits<double>::doubleMeasure(x._data[0]);
  }

  static void setFloat(PC& t, const int i, const double& val)
  {
    t._data[i]=val;
  }
  static double getFloat(const PC& t, const int i)
  {
    return t._data[i];
  }

  double stdDev() const
  {
    return _pcSet->StDv(_data);
  }

  // only printing the value and not the derivative
  static void write(FILE* fp, const PC& x) {throw;}

  static int getDimension() {return NumTypeTraits<double>::getDimension()+1;}
  
  static void getShape(int *shp) { *shp = N; NumTypeTraits<double>::getShape(shp+1);}

  static int getDataSize()
  {
    return N*NumTypeTraits<double>::getDataSize();
  }

  PC fabs() const
  {
    PC x(*this);
    x._data[0] = ::fabs(_data[0]);
    return x;
  }
    
  
  static void accumulateOneNorm(PC& sum, const PC& v) { sum += v.fabs();}
  
  static void accumulateDotProduct(PC& sum, const PC& v0, const PC& v1)
  {
    sum += v0*v1;
  }

  static void reduceSum(PC& sum, const PC& x) {sum+=x;}

  static void safeDivide(PC& x, const PC& y) {if (y._data[0]!=0) x/=y;}
  static void normalize(PC& x, const PC& y) {if (y._data[0]!=0) x/=y._data[0];}

  static void setMax(PC& x, const PC& y) {if (y._data[0]>x._data[0]) x._data[0]=y._data[0];}


  // the coefficients of expansion
  double _data[N];

  // the PCSet object from the UQ toolkit that implements the
  // algebraic and other operations for this particular version of PC
  // class
  
  static PCSet *_pcSet;
};


#define PC_BINARY_OP(opname,_op_)                                       \
  template<int ORDER, int DIM>                                          \
  PC<ORDER,DIM> opname(const PC<ORDER,DIM>& a, const PC<ORDER,DIM>& b)  \
  {                                                                     \
    return PC<ORDER,DIM>(a) _op_ b;                                     \
  }                                                                     \
  template<int ORDER, int DIM>                                          \
  PC<ORDER,DIM> opname(const PC<ORDER,DIM>& a, const double& b)         \
  {                                                                     \
    return PC<ORDER,DIM>(a) _op_ b;                                     \
  }                                                                     \
  template<int ORDER, int DIM>                                          \
  PC<ORDER,DIM> opname(const double& a, const PC<ORDER,DIM>& b)         \
  {                                                                     \
    return PC<ORDER,DIM>(a) _op_ b;                                     \
  }                                                                     \
 
PC_BINARY_OP(operator+,+=);
PC_BINARY_OP(operator-,-=);
PC_BINARY_OP(operator*,*=);
PC_BINARY_OP(operator/,/=);

template<int ORDER, int DIM>               
PC<ORDER,DIM> operator+(const PC<ORDER,DIM>& a)
{
  return a;
}

template<int ORDER, int DIM>               
PC<ORDER,DIM> operator-(const PC<ORDER,DIM>& a)
{
  PC<ORDER,DIM> z;
  return z-a;
}

template<int ORDER, int DIM>               
PC<ORDER,DIM> fabs(const PC<ORDER,DIM>& a)
{
  PC<ORDER,DIM> x(a);
  x._data[0]=fabs(a._data[0]);
  return x;
}

template<int ORDER, int DIM>               
inline ostream &operator<<(ostream &os, const PC<ORDER,DIM> &a)
{
  a.print(os);
  return os;
}

template<int ORDER, int DIM>
PC<ORDER,DIM>
sqrt(const PC<ORDER,DIM>& a)
{
  PC<ORDER,DIM> r;
  PC<ORDER,DIM>::_pcSet->RPow(a._data,r._data,0.5);
  return r;
}

template<int ORDER, int DIM>               
double stdDev(const PC<ORDER,DIM>& a)
{
  return PC<ORDER,DIM>::_pcSet->StDv(a._data);
}



template<int ORDER, int DIM>
PCSet* PC<ORDER,DIM>::_pcSet = createPCSet(ORDER,DIM);

template<>
template<int ORDER, int DIM>
struct ArrayScalarTraits<PC< ORDER, DIM > >
{
  static void limit(PC<ORDER,DIM>& val, const double min, const double max)
  {
    if (val._data[0] < min)
      val._data[0] = min;
    else if (val._data[0] > max)
      val._data[0] = max;
  }
};

#endif

