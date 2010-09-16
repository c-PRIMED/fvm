#include "StructureFields.h"

StructureFields::StructureFields(const string baseName) :
  deformation(baseName + ".deformation"),
  deformationGradient(baseName + ".deformationGradient"),
  deformationFlux(baseName + ".deformationFlux"),
  eta(baseName + ".eta"),
  eta1(baseName + ".eta1"),
  density(baseName + ".density"),
  deformationN1(baseName + ".deformationN1"),
  deformationN2(baseName + ".deformationN2"),
  tractionX(baseName + ".tractionX"),
  tractionY(baseName + ".tractionY"),
  tractionZ(baseName + ".tractionZ"),
  bodyForce(baseName + ".bodyForce"),
  volume0(baseName + ".volume0")
{}

