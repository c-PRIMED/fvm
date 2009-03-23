#include <atype.h>

#include "MeshMetricsCalculator.h"
#include "MeshMetricsCalculator_impl.h"

template class MeshMetricsCalculator<ATYPE>;
template class matrix<ATYPE>;

#include "ThermalModel.h"
#include "ThermalModel_impl.h"
template class ThermalModel<ATYPE>;


#include "FlowModel.h"
#include "FlowModel_impl.h"

template class FlowModel<ATYPE>;

