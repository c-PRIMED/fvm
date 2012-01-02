%{
#include "CometMatrix.h"
#include "MatrixJML.h"
  %}

%include "ArrayBase.i"

template<class T>
class CometMatrix: public MatrixJML<T>
{
 public:
  typedef Array<T> TArray;
  
  CometMatrix(const int order);
  
  T& getElement(const int i,const int j);
  void printElement(const int& i,const int& j);
  void Solve(TArray& bVec);
  
 private:
  int _elements;
  TArray _values;
  int _order;
  
};

%template(CometMatrixA) CometMatrix< ATYPE_STR >;
