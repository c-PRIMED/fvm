%{
#include "Tangent.h"
  %}

#if 0
%typemap(in) Tangent
{
  if (PyFloat_Check($input))
  {
      $1 = Tangent(PyFloat_AsDouble($input));
  }
  else if (PyInt_Check($input))
  {
      $1 = Tangent(double(PyInt_AsLong($input)));
  }
  else if (PyTuple_Check($input) && (2==PyTuple_Size($input)))
  {
      $1 = Tangent(PyFloat_AsDouble(PyTuple_GetItem($input,0)),
                   PyFloat_AsDouble(PyTuple_GetItem($input,1)));
  }
  else
    throw CException("invalid Tangent input");
}
#endif

class Tangent
{
public:
  Tangent();
  Tangent(double a, double b);
  double _v;
  double _dv;
};

typedef Tangent ATYPE;


//%include "atype.swg"

#define USING_ATYPE_TANGENT 

#define ATYPE_STR Tangent

