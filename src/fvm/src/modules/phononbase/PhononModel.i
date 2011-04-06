%{
#include "PhononModel.h"
  %}

%include "Model.i"
%include "FloatVarDict.i"
%include "PhononModel.h"

template <class T>
struct PhononBC : public FloatVarDict<T>
{
  string bcType;
}; 

template<class T>
struct PhononModelOptions : public FloatVarDict<T>
{
  bool printNormalizedResiduals;
  bool transient;
  int timeDiscretizationOrder;
  bool constantcp;
  bool constanttau;
  LinearSolver* PhononLinearSolver;
  double absTolerance;
  double relTolerance;
};

%template(PhononBCA) PhononBC<ATYPE_STR>;
%template(PhononModelOptionsA) PhononModelOptions<ATYPE_STR>;
%template(PhononBCList) std::vector<PhononBC< ATYPE_STR >* >;
%template(PhononBCsMap) std::map<int,PhononBC< ATYPE_STR >* >;
%template(PhononModelA) PhononModel< ATYPE_STR >;
