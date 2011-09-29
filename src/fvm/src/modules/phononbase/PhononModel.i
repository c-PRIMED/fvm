%{
#include "PhononModel.h"
  %}

%include "Model.i"
%include "FloatVarDict.i"

template<class T>
class PhononModel : public Model
{

 public:
  
  typedef map<int,PhononBC<T>*> PhononBCMap; 
  typedef Kspace<T> Tkspace;
  PhononModel(const MeshList& meshes,const GeomFields& geomFields,Tkspace& kspace,PhononMacro& macro);
  PhononModelOptions<T>& getOptions();
  PhononBCMap& getBCs();
  void init();
  void callBoundaryConditions();
  void updateTL();
  void COMETupdateTL();
  void updatee0();
  void updateHeatFlux();
  void advance(const int niter);
  void advanceCOMET(const int niter);
  void printTemp();
  T HeatFluxIntegral(const Mesh& mesh, const int faceGroupId);
   
  private:

    const GeomFields& _geomFields;
    Tkspace& _kspace;       //kspace
    PhononMacro& _macro;
    PhononModelOptions<T> _options;
    PhononBCMap _bcMap;
    MFRPtr _initialnorm;
    int _niters;
    
};


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
  LinearSolver* PhononLinearSolver;
  double absTolerance;
  double relTolerance;
  int showResidual;
};

%template(PhononBCA) PhononBC<ATYPE_STR>;
%template(PhononModelOptionsA) PhononModelOptions<ATYPE_STR>;
%template(PhononBCList) std::vector<PhononBC< ATYPE_STR >* >;
%template(PhononBCsMap) std::map<int,PhononBC< ATYPE_STR >* >;
%template(PhononModelA) PhononModel< ATYPE_STR >;
