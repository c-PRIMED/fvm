%{
#include "BatteryModel.h"
%}

using namespace std;

%include "FloatVarDict.i"

template <class T>
struct BatterySpeciesBC : public FloatVarDict<T>
{
  string bcType;
}; 

template <class T>
struct BatterySpeciesVC : public FloatVarDict<T>
{
  string vcType;
}; 
template <class T>
struct BatteryPotentialBC : public FloatVarDict<T>
{
  string bcType;
}; 

template <class T>
struct BatteryPotentialVC : public FloatVarDict<T>
{
  string vcType;
}; 

template <class T>
struct BatteryThermalBC : public FloatVarDict<T>
{
  string bcType;
}; 

template <class T>
struct BatteryThermalVC : public FloatVarDict<T>
{
  string vcType;
}; 

template <class T>
struct BatteryModelOptions : public FloatVarDict<T>
{
  double relativeTolerance;
  double absoluteTolerance;
  double relativeSpeciesTolerance;
  double absoluteSpeciesTolerance;
  double relativePotentialTolerance;
  double absolutePotentialTolerance;
  double relativeThermalTolerance;
  double absoluteThermalTolerance;
  double relativePCTolerance;
  double absolutePCTolerance;
  LinearSolver *linearSolver;
  LinearSolver *linearSolverSpecies;
  LinearSolver *linearSolverPotential;
  LinearSolver *linearSolverThermal;
  LinearSolver *linearSolverPC;
  bool useCentralDifference;
  bool transient;
  bool ButlerVolmer;
  bool thermalModelPC;
  int timeDiscretizationOrder;
  int advanceVerbosity;
}; 


%template(BatterySpeciesBCA) BatterySpeciesBC< ATYPE_STR >;
%template(BatterySpeciesVCA) BatterySpeciesVC< ATYPE_STR >;
%template(BatterySpeciesBCList) std::vector<BatterySpeciesBC< ATYPE_STR >* >;
%template(BatterySpeciesBCsMap) std::map<int,BatterySpeciesBC< ATYPE_STR >* >;
%template(BatterySpeciesVCList) std::vector<BatterySpeciesVC< ATYPE_STR >* >;
%template(BatterySpeciesVCsMap) std::map<int,BatterySpeciesVC< ATYPE_STR >* >;

%template(BatteryPotentialBCA) BatteryPotentialBC< ATYPE_STR >;
%template(BatteryPotentialVCA) BatteryPotentialVC< ATYPE_STR >;
%template(BatteryPotentialBCList) std::vector<BatteryPotentialBC< ATYPE_STR >* >;
%template(BatteryPotentialBCsMap) std::map<int,BatteryPotentialBC< ATYPE_STR >* >;
%template(BatteryPotentialVCList) std::vector<BatteryPotentialVC< ATYPE_STR >* >;
%template(BatteryPotentialVCsMap) std::map<int,BatteryPotentialVC< ATYPE_STR >* >;

%template(BatteryThermalBCA) BatteryThermalBC< ATYPE_STR >;
%template(BatteryThermalVCA) BatteryThermalVC< ATYPE_STR >;
%template(BatteryThermalBCList) std::vector<BatteryThermalBC< ATYPE_STR >* >;
%template(BatteryThermalBCsMap) std::map<int,BatteryThermalBC< ATYPE_STR >* >;
%template(BatteryThermalVCList) std::vector<BatteryThermalVC< ATYPE_STR >* >;
%template(BatteryThermalVCsMap) std::map<int,BatteryThermalVC< ATYPE_STR >* >;

%template(BatteryModelOptionsA) BatteryModelOptions< ATYPE_STR >;


%import "Model.i"

%include "BatteryModel.h"


%template(BatteryModelA) BatteryModel< ATYPE_STR >;

