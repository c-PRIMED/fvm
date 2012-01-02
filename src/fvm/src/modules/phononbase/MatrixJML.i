%{
#include "MatrixJML.h"
  %}

template<class T>
class MatrixJML
{
 public:

  MatrixJML() {}
  virtual ~MatrixJML() {}
  virtual T& getElement(const int i,const int j);
  virtual void zero()=0;

};

%template(MatrixJMLA) MatrixJML< ATYPE_STR >;
