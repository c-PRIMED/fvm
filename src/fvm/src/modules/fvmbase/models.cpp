#include <atype.h>

#include "MeshMetricsCalculator.h"
#include "MeshMetricsCalculator_impl.h"

template class MeshMetricsCalculator<ATYPE>;
template class matrix<ATYPE>;

#include "ThermalModel.h"
#include "ThermalModel_impl.h"
template class ThermalModel<ATYPE>;

#include "ElectricModel.h"
#include "ElectricModel_impl.h"
template class ElectricModel<ATYPE>;

#include "FlowModel.h"
#include "FlowModel_impl.h"

template class FlowModel<ATYPE>;

#include "IdealGasDensityModel.h"
#include "IdealGasDensityModel_impl.h"

template class IdealGasDensityModel<ATYPE>;

#include "MovingMeshModel.h"

template class MovingMeshModel<ATYPE>;

#include "StructureModel.h"
#include "StructureModel_impl.h"

template class StructureModel<ATYPE>;

