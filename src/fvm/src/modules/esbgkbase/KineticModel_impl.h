#include "Mesh.h"

#include "NumType.h"
#include "Array.h"
#include "Field.h"
#include "CRConnectivity.h"
#include "LinearSystem.h"
#include "StorageSite.h"
#include "MultiFieldMatrix.h"

#include "Vector.h"
#include "VectorTranspose.h"

#include "AMG.h"
#include "Linearizer.h"
#include "GradientModel.h"
#include "GenericIBDiscretization.h"
#include "StressTensor.h"

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

template<class T>
class FlowModel<T>::Impl
{
public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef VectorTranspose<T,3> VectorT3T;
  
  typedef Array<VectorT3> VectorT3Array;
  typedef DiagonalTensor<T,3> DiagTensorT3;
  

  
  
  }
private:
  const MeshList _meshes;
  macroFields& _macroFields;
 FlowBCMap _bcMap;
  FlowVCMap _vcMap;
}
