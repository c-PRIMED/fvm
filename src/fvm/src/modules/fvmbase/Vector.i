%{
#include "Vector.h"
%}

#ifdef USING_ATYPE_TANGENT

%include "atype.i"

#endif

template <class T, int N>
class Vector
{
public:
  %extend
  {
    T __getitem__(int i) {return (*$self)[i];}

    void __setitem__(int i, double x) {(*$self)[i]=x;}

  };
};
