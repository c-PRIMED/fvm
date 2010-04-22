%{
#include "TwobyTwoTensor.h"
%}

#ifdef USING_ATYPE_TANGENT
%include "atype.i"
#endif

template <class T>
class TwobyTwoTensor
{
 public:
  static string getTypeName();
  static int getDimension();
  static void getShape(int *shp);
  static int getDataSize();
  void zero();
  static TwobyTwoTensor getZero();
  static double doubleMeasure(const TwobyTwoTensor& x);
  static TwobyTwoTensor getNegativeUnity();
  static TwobyTwoTensor getUnity();
  TwobyTwoTensor operator-();

  %extend{
     T __getitem__(int i)  {return (*$self)[i];}
     void __setitem__(int i, double x)  { (*$self)[i]=x;}
  };
};

