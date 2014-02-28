%{
#include "FractureModel.h"
%}

using namespace std;

%include "FloatVarDict.i"
%include "Model.i"

template <class T>
struct FractureBC : public FloatVarDict<T>
{
  string bcType;
}; 

template <class T>
struct FractureVC : public FloatVarDict<T>
{
  string vcType;
}; 

template <class T>
struct FractureModelOptions : public FloatVarDict<T>
{
  double relativeTolerance;
  double absoluteTolerance;
  LinearSolver *linearSolver;
  bool transient;
}; 


%template(FractureBCA) FractureBC< ATYPE_STR >;
%template(FractureVCA) FractureVC< ATYPE_STR >;
%template(FractureBCList) std::vector<FractureBC< ATYPE_STR >* >;
%template(FractureBCsMap) std::map<int,FractureBC< ATYPE_STR >* >;
%template(FractureVCList) std::vector<FractureVC< ATYPE_STR >* >;
%template(FractureVCsMap) std::map<int,FractureVC< ATYPE_STR >* >;

%template(FractureModelOptionsA) FractureModelOptions< ATYPE_STR >;


%import "Model.i"

%include "FractureModel.h"


%template(FractureModelA) FractureModel< ATYPE_STR >;

