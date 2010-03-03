#include "StructureFields.h"

StructureFields::StructureFields(const string baseName) :
  deformation(baseName + ".deformation"),
  deformationGradient(baseName + ".deformationGradient"),
  deformationFlux(baseName + ".deformationFlux"),
  eta(baseName + ".eta"),
  lambda(baseName + ".lambda"),
  density(baseName + ".density"),
  deformationN1(baseName + ".deformationN1"),
  deformationN2(baseName + ".deformationN2")
{}

