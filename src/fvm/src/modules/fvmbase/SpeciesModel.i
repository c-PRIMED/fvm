%{
#include "SpeciesModel.h"
%}

using namespace std;

%include "FloatVarDict.i"

template <class T>
struct SpeciesBC : public FloatVarDict<T>
{
  string bcType;
}; 

template <class T>
struct SpeciesVC : public FloatVarDict<T>
{
  string vcType;
}; 

template <class T>
struct SpeciesModelOptions : public FloatVarDict<T>
{
  double relativeTolerance;
  double absoluteTolerance;
  LinearSolver *linearSolver;
  bool useCentralDifference;
  bool transient;
  bool ButlerVolmer;
  int timeDiscretizationOrder;
}; 


%template(SpeciesBCA) SpeciesBC< ATYPE_STR >;
%template(SpeciesVCA) SpeciesVC< ATYPE_STR >;
%template(SpeciesBCList) std::vector<SpeciesBC< ATYPE_STR >* >;
%template(SpeciesBCsMap) std::map<int,SpeciesBC< ATYPE_STR >* >;
%template(SpeciesVCList) std::vector<SpeciesVC< ATYPE_STR >* >;
%template(SpeciesVCsMap) std::map<int,SpeciesVC< ATYPE_STR >* >;

%template(SpeciesModelOptionsA) SpeciesModelOptions< ATYPE_STR >;


%import "Model.i"

%include "SpeciesModel.h"


%template(SpeciesModelA) SpeciesModel< ATYPE_STR >;

