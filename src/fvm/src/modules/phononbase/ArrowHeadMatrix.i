%{
#include "ArrowHeadMatrix.h"
#include "MatrixJML.h"
  %}

%include "ArrayBase.i"

template<class T>
class ArrowHeadMatrix: public MatrixJML<T>
{
 public:
  typedef Array<T> TArray;
  
  ArrowHeadMatrix(const int order);
  T& getElement(const int i,const int j);
  void zero();
  void ones();
  void printElement(const int& i,const int& j);
  void Solve(TArray& bVec);
  void DummySolve();
  
 private:
  int _elements;
  TArray _values;
  int _order;
  
};

%template(ArrowHeadMatrixA) ArrowHeadMatrix< ATYPE_STR >;
