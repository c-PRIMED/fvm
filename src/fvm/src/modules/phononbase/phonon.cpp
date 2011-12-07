#include <atype.h>


#include "Kspace.h"
template class Kspace<ATYPE>;

#include "PhononModel.h"
template class PhononModel<ATYPE>;

#include "PhononBoundary.h"
template class PhononBoundary<ATYPE>;

#include "PhononInterface.h"
template class PhononInterface<ATYPE>;

#include "ArrowHeadMatrix.h"
template class ArrowHeadMatrix<ATYPE>;

#include "COMETModel.h"
template class COMETModel<ATYPE>;
