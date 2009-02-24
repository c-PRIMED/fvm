%{
#include "MeshMetricsCalculator.h"
#include "atype.h"
%}

using namespace std;
%include "Model.i"
%include "GeomFields.h"
%include "atype.i"


template<class T>
class MeshMetricsCalculator : public Model
{
public:
  MeshMetricsCalculator(GeomFields& geomFields, const MeshList& meshes);
  virtual void init();
  void computeIBInterpolationMatrices(const StorageSite& particles);
#ifdef USING_ATYPE_TANGENT
  void setTangentCoords(int meshID, int faceZoneID, int dim);
#endif
};


%template(MeshMetricsCalculatorA) MeshMetricsCalculator<ATYPE_STR>;

