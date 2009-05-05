%{
#include "Vector.h"
%}

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
