%{
#include "MeshMetricsCalculator.h"
%}

using namespace std;
%include "Model.i"
template<class T>
class MeshMetricsCalculator : public Model
{
public:
  MeshMetricsCalculator(const MeshList& meshes);
  virtual void init();
};

%template(MeshMetricsCalculator_double) MeshMetricsCalculator<double>;

