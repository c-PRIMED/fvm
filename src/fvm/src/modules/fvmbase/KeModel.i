%{
#include "KeModel.h"
%}

using namespace std;

%include "FloatVarDict.i"

template <class T>
struct KeBC : public FloatVarDict<T>
{
  string bcType;
};

template <class T>
struct KeVC : public FloatVarDict<T>
{
  string vcType;
};

template <class T>
struct KeModelOptions : public FloatVarDict<T>
{
  bool transient;
  int timeDiscretizationOrder;
  double relativeTolerance;
  double absoluteTolerance;
  LinearSolver *linearSolver;
  double cmu;
  double c2mu;
  double c1mu;
  double vk;
  double emp;
  double sigmak;
  double sigmae;
  bool useCentralDifference;
  bool printNormalizedResiduals;

};


%template(KeBCA) KeBC< ATYPE_STR >;
%template(KeVCA) KeVC< ATYPE_STR >;
%template(KeBCList) std::vector<KeBC< ATYPE_STR >* >;
%template(KeBCsMap) std::map<int,KeBC< ATYPE_STR >* >;
%template(KeVCList) std::vector<KeVC< ATYPE_STR >* >;
%template(KeVCsMap) std::map<int,KeVC< ATYPE_STR >* >;

%template(KeModelOptionsA) KeModelOptions< ATYPE_STR >;


%import "Model.i"

%include "KeModel.h"


%template(KeModelA) KeModel< ATYPE_STR >;
