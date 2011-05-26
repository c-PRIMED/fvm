%{
#include "DistFunctFields.h"
  %}

//%include "Vector.i"

%include "std_vector.i"
%include "std_map.i"
%include "std_set.i"


using namespace std;

class DistFunctFields
{
 public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  
  //std::vector<shared_ptr<Field> > DistFunctFieldsPtr;
  std::vector<Field*> dsf;
  
  %extend{
   std::vector<Field*> getDSFVec(){
   
       self->dsf
   }
     
  }
  
};


