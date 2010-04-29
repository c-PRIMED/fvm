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
  bool transient_enable;
  bool tunneling_enable;
  bool chargetransport_enable;
  bool electrostatics_enable;
  bool ibm;
  bool capture_enable;
  bool emission_enable;
  bool drift_enable;
  bool injection_enable;
  bool drift_enable;
  bool diffusion_enable;  

  LinearSolver *electrostaticsLinearSolver;
  LinearSolver *chargetransportLinearSolver;
}; 

template <class T>
struct ElectricModelConstants : public FloatVarDict<T>
{
 
};


%template(ElectricBCA) ElectricBC< ATYPE_STR >;
%template(ElectricBCList) std::vector<ElectricBC< ATYPE_STR >* >;
%template(ElectricBCsMap) std::map<int,ElectricBC< ATYPE_STR >* >;


%template(ElectricVCA) ElectricVC< ATYPE_STR >;
%template(ElectricVCList) std::vector<ElectricVC< ATYPE_STR >* >;
%template(ElectricVCsMap) std::map<int,ElectricVC< ATYPE_STR >* >;


%template(ElectricModelOptionsA) ElectricModelOptions< ATYPE_STR >;
%template(ElectricModelConstantsA) ElectricModelConstants<ATYPE_STR>;

%import "Model.i"

%include "boost_shared_ptr.i"
SWIG_SHARED_PTR(ArrayVec3Ptr,Array<Vector< ATYPE_STR,3> >)

%include "ElectricModel.h"


%template(ElectricModelA) ElectricModel< ATYPE_STR >;

