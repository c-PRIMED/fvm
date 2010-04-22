%{
#include "MatrixOperation.h"
%}
using namespace std;
template <class T> 
class matrix{
 public:
   matrix();
   ~matrix();	
   double detMatrix (T matrix[10][10], const int size);

%template (matrixA) matrix<ATYPE_STR>;
};