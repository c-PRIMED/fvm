#ifndef _TANGENT_H_
#define _TANGENT_H_

#include "NumType.h"

#include "misc.h"

class Tangent 
{
public:

  typedef Tangent This_T;
  typedef Tangent T_Scalar;
  typedef NumTypeTraits<double>::T_BuiltIn T_BuiltIn;

  static string getTypeName()
  {
    return "Tangent";
  }
  
  // no initialization for empty ctor - consistent with builtins
  Tangent()
  {}

  
  
  Tangent(const double v, const double dv) :
    _v(v),
    _dv(dv)
  {}

  explicit Tangent(const double v) :
    _v(v),
    _dv(NumTypeTraits<double>::getZero())
  {}


  //to allow statements like T x(0)
  explicit Tangent(const int v) :
    _v(double(v)),
    _dv(NumTypeTraits<double>::getZero())
  {}

  Tangent(const Tangent& o) :
    _v(o._v),
    _dv(o._dv)
  {}
  
  ~Tangent()
  {}

  
  Tangent& operator=(const Tangent& o)
  {
    if (this == &o)
      return *this;
    _v = o._v;
    _dv = o._dv;
    return *this;
  }

  Tangent& operator=(const double& f)
  {
    _v = f;
    _dv = NumTypeTraits<double>::getZero();
    return *this;
  }

  // to allow statements like x = 0;
  Tangent& operator=(const int& i)
  {
    _v = double(i);
    _dv = NumTypeTraits<double>::getZero();
    return *this;
  }

  Tangent& operator+=(const Tangent& o)
  {
    _v += o._v;
    _dv += o._dv;
    return *this;
  }

  Tangent& operator-=(const Tangent& o)
  {
    _v -= o._v;
    _dv -= o._dv;
    return *this;
  }

  Tangent& operator*=(const Tangent& o)
  {
    _dv = _dv*o._v + _v*o._dv;
    _v *= o._v;
    return *this;
  }

  Tangent& operator/=(const Tangent& o)
  {
    _dv = (_dv*o._v - _v*o._dv)/(o._v*o._v);
    _v /= o._v;
    return *this;
  }

  Tangent& operator+=(const double& o)
  {
    _v += o;
    return *this;
  }

  Tangent& operator-=(const double& o)
  {
    _v -= o;
    return *this;
  }

  Tangent& operator*=(const double& o)
  {
    _v *= o;
    _dv *= o;
    return *this;
  }

  Tangent& operator/=(const double& o)
  {
    _v /= o;
    _dv /= o;
    return *this;
  }

  Tangent& operator+=(const int& i)
  {
    _v += double(i);
    return *this;
  }
  
  Tangent& operator-=(const int& i)
  {
    _v -= double(i);
    return *this;
  }
  
  Tangent& operator*=(const int& i)
  {
    _v *= double(i);
    _dv *= double(i);
    return *this;
  }
  
  Tangent& operator/=(const int& i)
  {
    _v /= double(i);
    _dv /= double(i);
    return *this;
  }

#define TANGENT_RELATIONAL_OPS(opname,_op_)     \
  bool opname(const Tangent& o) const           \
  {                                             \
    return (_v _op_ o._v);                      \
  }                                             \
    bool opname(const double& o) const          \
  {                                             \
    return (_v _op_ o);                         \
  }                                             \
    bool opname(const int& o) const             \
  {                                             \
    return (_v _op_ double(o));                 \
  }

  TANGENT_RELATIONAL_OPS(operator>,>);
  TANGENT_RELATIONAL_OPS(operator>=,>=);
  TANGENT_RELATIONAL_OPS(operator<,<);
  TANGENT_RELATIONAL_OPS(operator<=,<=);
  TANGENT_RELATIONAL_OPS(operator==,==);
  TANGENT_RELATIONAL_OPS(operator!=,!=);
  
#undef TANGENT_RELATIONAL_OPS
  
  void print(ostream &os) const
  {
    os << "< "  << _v << " , " << _dv << ">";
  }

  static Tangent getZero()
  {
    double zero = NumTypeTraits<double>::getZero();
    return Tangent(zero,zero);
  }

  static Tangent getUnity()
  {
    return Tangent(NumTypeTraits<double>::getUnity(),
                   NumTypeTraits<double>::getZero());
  }

  static Tangent getNegativeUnity()
  {
    return Tangent(NumTypeTraits<double>::getNegativeUnity(),
                   NumTypeTraits<double>::getZero());
  }

  static double doubleMeasure(const Tangent& x)
  {
    return NumTypeTraits<double>::doubleMeasure(x._v);
  }

  static void setFloat(Tangent& t, const int i, const double& val)
  {
    if (i == 0)
      t._v = val;
    else
      t._dv = val;
  }
  static double getFloat(const Tangent& t, const int i)
  {
    if (i==0)
      return t._v;
    else
      return t._dv;
  }


  // only printing the value and not the derivative
  static void write(FILE* fp, const Tangent& x) {fprintf(fp,"%f",x._v);}

  static int getDimension() {return NumTypeTraits<double>::getDimension()+1;}
  
  static void getShape(int *shp) { *shp = 2; NumTypeTraits<double>::getShape(shp+1);}

  static int getDataSize()
  {
    return 2*NumTypeTraits<double>::getDataSize();
  }

  Tangent fabs() const
  {
    return Tangent(::fabs(_v), (_v > 0 ? _dv : -_dv));
  }
    
  
  static void accumulateOneNorm(Tangent& sum, const Tangent& v) { sum += v.fabs();}
  
  static void accumulateDotProduct(Tangent& sum, const Tangent& v0, const Tangent& v1)
  {
    sum._v += v0._v*v1._v;
    sum._dv += v0._v*v1._dv + v0._dv*v1._v;
  }

  
  double _v;
  double _dv;
};


#define TANGENT_BINARY_OP(opname,_op_)                  \
  Tangent opname(const Tangent& a, const Tangent& b)    \
  {                                                     \
    return Tangent(a) _op_ b;                            \
  }                                                      \
  Tangent opname(const Tangent& a, const double& b)      \
  {                                                      \
    return Tangent(a) _op_ b;                            \
  }                                                      \
  Tangent opname(const double& a, const Tangent& b)      \
  {                                                      \
    return Tangent(b) _op_ a;                            \
  }                                                      \
  Tangent opname(const Tangent& a, const int& b)         \
  {                                                      \
    return Tangent(a) _op_ b;                            \
  }                                                      \
  Tangent opname(const int& a, const Tangent& b)         \
  {                                                      \
    return Tangent(b) _op_ a;                            \
  }                                                      \
                                                         
TANGENT_BINARY_OP(operator+,+=);
TANGENT_BINARY_OP(operator-,-=);
TANGENT_BINARY_OP(operator*,*=);
TANGENT_BINARY_OP(operator/,/=);

Tangent operator+(const Tangent& a)
{
  return a;
}

Tangent operator-(const Tangent& a)
{
  return Tangent(-a._v,-a._dv);
}

inline ostream &operator<<(ostream &os, const Tangent &a)
{
  a.print(os);
  return os;
}

Tangent sin(const Tangent& a)
{
  return Tangent(sin(a._v), a._dv*cos(a._v));
}

Tangent fabs(const Tangent& a)
{
  return Tangent(fabs(a._v), (a._v > 0 ? a._dv : -a._dv));
}

Tangent sqrt(const Tangent& a)
{
  double sqv = sqrt(a._v);
  return Tangent(sqv, sqv==0.0 ? 0 : 0.5*a._dv/sqv);
}


#endif
