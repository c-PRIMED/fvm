#include <atype.h>


#include "Kspace.h"
template class Kspace<ATYPE>;

#include "DensityOfStates.h"
template class DensityOfStates<ATYPE>;

#include "PhononModel.h"
template class PhononModel<ATYPE>;

#include "PhononBoundary.h"
template class PhononBoundary<ATYPE>;

#include "PhononInterface.h"
template class PhononInterface<ATYPE>;

#include "ArrowHeadMatrix.h"
template class ArrowHeadMatrix<ATYPE>;

#include "SquareMatrix.h"
template class SquareMatrix<ATYPE>;

#include "COMETModel.h"
template class COMETModel<ATYPE>;
