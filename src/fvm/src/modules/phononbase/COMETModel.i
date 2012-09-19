%{
#include "COMETModel.h"
  %}

%include "Model.i"
%include "FloatVarDict.i"

template<class T>
class COMETModel : public Model
{

 public:
  
  typedef Kspace<T> Tkspace;
  typedef map<int,COMETBC<T>*> COMETBCMap;
  typedef Vector<T,3> VectorT3;
  typedef COMETModel<T> TCOMET;
  typedef shared_ptr<TCOMET> TCOMETPtr;
  typedef Tkspace* TkspPtr;
  typedef vector<TkspPtr> TkspList;
  typedef map<int,int> MeshKspaceMap;

  COMETModel(const MeshList& meshes,const int level,GeomFields& geomFields,
	     TkspList& kspaces,PhononMacro& macro);
  
  void init();
  void updateTL();
  void advance(int);
  T HeatFluxIntegral(const Mesh& mesh, const int faceGroupId);
  T getWallArea(const Mesh& mesh, const int faceGroupId);
  VectorT3 getWallAreaVector(const Mesh& mesh, const int faceGroupId);
  COMETModelOptions<T>& getOptions();
  COMETBCMap& getBCs();
  MeshKspaceMap& getMKMap();
  T getResidual();
  ArrayBase* getValueArray(const Mesh& mesh, const int cell);
  void equilibrate();
  ArrayBase& getLatticeTemp(const Mesh& mesh);
  void setStraightLine(const T T1,const T T2);
  void applyTemperatureBoundaries();
  void initFromOld();
  T calcDomainStats();
  void makeCellColors(const int level);
  void makePlotColors(const int level);

 private:
  
  const int _level; 
  GeomFields& _geomFields;
  TkspList& _kspaces;
  PhononMacro& _macro;
  TCOMET* _finestLevel;
  TCOMET* _coarserLevel;
  TCOMET* _finerLevel;
  PhononBCMap _bcMap;
  COMETModelOptions<T> _options;
  BCfaceList _BFaces;
  BCcellList _BCells;
  T _residual;
  MeshKspaceMap _MeshKspaceMap;
  IClist _IClist;
  MeshICmap _MeshToIC;
  int _rank;
 
};

template <class T>
struct COMETBC : public FloatVarDict<T>
{
  string bcType;
  string InterfaceModel;
}; 

template<class T>
struct COMETModelOptions : public FloatVarDict<T>
{
  bool printNormalizedResiduals;
  bool transient;
  int timeDiscretizationOrder;
  double absTolerance;
  double relTolerance;
  int showResidual;
  int maxLevels;
  string AgglomerationMethod;
  int preSweeps;
  int postSweeps;
  double relFactor;
  bool withNormal;
  double NewtonTol;
  int BinsPerBand;
  string Convection;
  string DomainStats;
};

%template(COMETBCA) COMETBC<ATYPE_STR>;
%template(COMETModelOptionsA) COMETModelOptions<ATYPE_STR>;
%template(COMETBCList) std::vector<COMETBC< ATYPE_STR >* >;
%template(COMETBCsMap) std::map<int,COMETBC< ATYPE_STR >* >;
%template(COMETModelA) COMETModel<ATYPE_STR>;
%template(COMETMKmap) std::map<int,int>;
