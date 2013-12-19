%{
#include "VacancyModel.h"
%}

using namespace std;

%include "FloatVarDict.i"
%include "Model.i"

template <class T>
struct VacancyBC : public FloatVarDict<T>
{
  string bcType;
}; 

template <class T>
struct VacancyVC : public FloatVarDict<T>
{
  string vcType;
}; 

template <class T>
struct VacancyModelOptions : public FloatVarDict<T>
{
  double relativeTolerance;
  double absoluteTolerance;
  LinearSolver *linearSolver;
  bool useCentralDifference;
  bool transient;
}; 


%template(VacancyBCA) VacancyBC< ATYPE_STR >;
%template(VacancyVCA) VacancyVC< ATYPE_STR >;
%template(VacancyBCList) std::vector<VacancyBC< ATYPE_STR >* >;
%template(VacancyBCsMap) std::map<int,VacancyBC< ATYPE_STR >* >;
%template(VacancyVCList) std::vector<VacancyVC< ATYPE_STR >* >;
%template(VacancyVCsMap) std::map<int,VacancyVC< ATYPE_STR >* >;

%template(VacancyModelOptionsA) VacancyModelOptions< ATYPE_STR >;


%import "Model.i"

%include "VacancyModel.h"


%template(VacancyModelA) VacancyModel< ATYPE_STR >;

