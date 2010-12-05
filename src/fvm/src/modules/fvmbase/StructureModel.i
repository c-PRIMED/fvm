%{
#include "StructureModel.h"
#include "AMG.h"
%}


using namespace std;
using namespace boost;


%include "FloatVarDict.i"

template <class T>
struct StructureBC : public FloatVarDict<T>
{
  string bcType;
}; 

template <class T>
struct StructureVC : public FloatVarDict<T>
{
}; 

template <class T>
struct StructureModelOptions : public FloatVarDict<T>
{
  bool printNormalizedResiduals;
  bool transient;
  bool incompressible;
  double deformationTolerance;
  int timeDiscretizationOrder;
  LinearSolver *deformationLinearSolver;
  LinearSolver *coupledLinearSolver;
}; 

//%template(Vector3) Vector<ATYPE_STR,3>;

%template(StructureBCA) StructureBC< ATYPE_STR >;
%template(StructureBCList) std::vector<StructureBC< ATYPE_STR >* >;
%template(StructureBCsMap) std::map<int,StructureBC< ATYPE_STR >* >;

%template(StructureVCA) StructureVC< ATYPE_STR >;
%template(StructureVCList) std::vector<StructureVC< ATYPE_STR >* >;
%template(StructureVCsMap) std::map<int,StructureVC< ATYPE_STR >* >;

%template(StructureModelOptionsA) StructureModelOptions< ATYPE_STR >;

%include "StructureModel.h"


%template(StructureModelA) StructureModel< ATYPE_STR >;

