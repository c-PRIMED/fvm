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
};


%template(MeshMetricsCalculatorA) MeshMetricsCalculator<ATYPE_STR>;

