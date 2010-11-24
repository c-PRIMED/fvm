%{
#include "MovingMeshModel.h"
%}


using namespace std;
using namespace boost;

%include "FloatVarDict.i"


template <class T>
struct MovingMeshModelOptions : public FloatVarDict<T>
{
  double absTolerance;
  double relativeTolerance;
  int nNodeDisplacementSweeps;
  int timeDiscretizationOrder;
}; 

//%template(Vector3) Vector<ATYPE_STR,3>;

%template(MovingMeshModelOptionsA) MovingMeshModelOptions< ATYPE_STR >;


%include "MovingMeshModel.h"


%template(MovingMeshModelA) MovingMeshModel< ATYPE_STR >;

