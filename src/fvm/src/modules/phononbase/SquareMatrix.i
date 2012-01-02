%{
#include "SquareMatrix.h"
#include "MatrixJML.h"
  %}

%include "ArrayBase.i"

template<class T>
class SquareMatrix: public MatrixJML<T>
{
 public:
  typedef Array<T> TArray;
  typedef shared_ptr<TArray> TArrPtr;
  typedef Array<int> IntArray;
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;

  SquareMatrix(const int N);
  T& getElement(const int i, const int j);
  T& operator()(const int i, const int j);
  void zero();
  void Solve(TArray& bVec);
  void makeCopy(SquareMatrix<T>& o);
  void printMatrix();
  T getTraceAbs();
  void multiply(const TArray& x, TArray& b);
  void testSolve();

%extend{
    void setElement(const int i, const int j, T value)
    {self->getElement(i,j)=value;}

    ArrayBase* getNewArray(const int size)
      {
	Array<T>* newArray=new Array<T>(size);
	newArray->zero();
	return newArray;
      }
  }
   
 private:
  const int _order;
  const int _elements;
  bool _sorted;
  IntArray _pivotRows;
  TArray _maxVals;
  TArray _values;
 
};

%template(SquareMatrixA) SquareMatrix< ATYPE_STR >;
