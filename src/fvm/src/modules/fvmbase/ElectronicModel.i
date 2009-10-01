%{
#include "ElectronicModel.h"
%}

using namespace std;

%include "FloatVarDict.i"

template <class T>
struct ElectronicBC : public FloatVarDict<T>
{
  string bcType;
}; 

template <class T>
struct ElectronicModelOptions : public FloatVarDict<T>
{
  double electrostaticsTolerance;
  double chargetransportTolerance;
  bool transient;
  bool tunneling;
  bool chargetransport;
  bool electrostatics;
  bool ibm;

  LinearSolver *electrostaticsLinearSolver;
  LinearSolver *chargetransportLinearSolver;
}; 


%template(ElectronicBCA) ElectronicBC<ATYPE_STR>;
%template(ElectronicBCList) std::vector<ElectronicBC<ATYPE_STR>* >;
%template(ElectronicBCsMap) std::map<int,ElectronicBC<ATYPE_STR>* >;

%template(ElectronicModelOptionsA) ElectronicModelOptions<ATYPE_STR>;


%import "Model.i"

%include "boost_shared_ptr.i"
SWIG_SHARED_PTR(ArrayVec3Ptr,Array<Vector<ATYPE_STR,3> >)

%include "ElectronicModel.h"


%template(ElectronicModelA) ElectronicModel<ATYPE_STR>;

