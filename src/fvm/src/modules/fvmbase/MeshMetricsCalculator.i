%{
#include "MeshMetricsCalculator.h"
#include "atype.h"
#include "Model.h"
%}

using namespace std;

%include "Model.i"

template<class T>
class MeshMetricsCalculator : public Model
{
public:
  MeshMetricsCalculator(GeomFields& geomFields, const MeshList& meshes, bool transient=False);
  virtual void init();
  void calculateBoundaryNodeNormal();
  void recalculate();
  void updateTime();
  void recalculate_deform();
  void computeIBInterpolationMatrices(const StorageSite& particles, const int option=0);
  void computeIBInterpolationMatricesCells();
  
  void eraseIBInterpolationMatrices(const StorageSite& particles);
  void computeSolidInterpolationMatrices(const StorageSite& particles);
  void computeIBandSolidInterpolationMatrices(const StorageSite& particles);
  void computeGridInterpolationMatrices(const StorageSite& grids, const StorageSite& faces );	
#ifdef USING_ATYPE_TANGENT
  void setTangentCoords(int meshID, int faceZoneID, int dim);
#endif
};


%template(MeshMetricsCalculatorA) MeshMetricsCalculator< ATYPE_STR >;

