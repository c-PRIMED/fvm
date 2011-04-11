%{
#include "MatrixOperation.h"
%}


using namespace std;

template <class T, int N>
class SquareMatrix {
 public:
  SquareMatrix();
  ~SquareMatrix();
  SquareMatrix(const T& s);
  T& operator()(int i, int j);
 

  %extend
  {
    T getVal(int i, int j) {return (*$self)(i,j);}

    void setVal(int i, int j, double x) {(*$self)(i,j)=x;}

    void printA() {(*$self).print();}

  };
};


template<class T, int N>
SquareMatrix<T,N> inverseGauss (SquareMatrix<T, N>& a, int size);

template<class T, int N>
SquareMatrix<T,N>  inverse(SquareMatrix<T,N>& b, int size);


%template(SquareMatrixA) SquareMatrix<double, 10>;
%template(inverseGaussA) inverseGauss<double, 10>;
%template(inverseA) inverse<double, 3>;
