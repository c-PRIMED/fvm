%{
#include "ElectricModel.h"
%}

using namespace std;

%include "FloatVarDict.i"

template <class T>
struct ElectricBC : public FloatVarDict<T>
{
  string bcType;
}; 

template <class T>
struct ElectricVC : public FloatVarDict<T>
{
  string vcType;
}; 

template <class T>
struct ElectricModelOptions : public FloatVarDict<T>
{
  double electrostaticsTolerance;
  double chargetransportTolerance;
  bool transient;
  bool tunneling;
  bool chargetransport;
  bool electrostatics_enable;
  bool ibm;

  LinearSolver *electrostaticsLinearSolver;
  LinearSolver *chargetransportLinearSolver;
}; 


%template(ElectricBCA) ElectricBC<ATYPE_STR>;
%template(ElectricBCList) std::vector<ElectricBC<ATYPE_STR>* >;
%template(ElectricBCsMap) std::map<int,ElectricBC<ATYPE_STR>* >;

%template(ElectricVCA) ElectricVC<ATYPE_STR>;
%template(ElectricVCList) std::vector<ElectricVC<ATYPE_STR>* >;
%template(ElectricVCsMap) std::map<int,ElectricVC<ATYPE_STR>* >;


%template(ElectricModelOptionsA) ElectricModelOptions<ATYPE_STR>;


%import "Model.i"

%include "boost_shared_ptr.i"
SWIG_SHARED_PTR(ArrayVec3Ptr,Array<Vector<ATYPE_STR,3> >)

%include "ElectricModel.h"


%template(ElectricModelA) ElectricModel<ATYPE_STR>;

