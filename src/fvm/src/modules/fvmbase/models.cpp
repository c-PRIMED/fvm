#include <atype.h>

#include "MeshMetricsCalculator.h"
#include "MeshMetricsCalculator_impl.h"

template class MeshMetricsCalculator<ATYPE>;


#include "ThermalModel.h"
#include "ThermalModel_impl.h"
template class ThermalModel<ATYPE>;

#ifndef USING_ATYPE_PC
#include "ElectricModel.h"
#include "ElectricModel_impl.h"
template class ElectricModel<ATYPE>;

#endif

#ifndef USING_ATYPE_PC
#include "ContactModel.h"
#include "ContactModel_impl.h"
template class ContactModel<ATYPE>;

#endif


#include "FlowModel.h"
#include "FlowModel_impl.h"

template class FlowModel<ATYPE>;

#include "IdealGasDensityModel.h"
#include "IdealGasDensityModel_impl.h"

template class IdealGasDensityModel<ATYPE>;

#include "RosselandModel.h"
#include "RosselandModel_impl.h"

template class RosselandModel<ATYPE>;



//#include "OneDConduction.h"

//template class OneDConduction<ATYPE>;

//#include "ShockTube.h"
//template class ShockTube<ATYPE>;

#ifndef USING_ATYPE_PC
#include "StructureModel.h"
#include "StructureModel_impl.h"

template class StructureModel<ATYPE>;

#include "PlateModel.h"
#include "PlateModel_impl.h"

template class PlateModel<ATYPE>;

#include "MovingMeshModel.h"

template class MovingMeshModel<ATYPE>;


#include "StructureDeformationModel.h"
template class StructureDeformationModel<ATYPE>;

#endif
