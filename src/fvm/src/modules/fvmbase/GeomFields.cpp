#include "GeomFields.h"

GeomFields::GeomFields(const string baseName) :
  coordinate(baseName+"coordinate"),
  coordinateN1(baseName+"coordinateN1"),
  area(baseName+"area"),
  areaN1(baseName+"areaN1"),
  areaMag(baseName+"areaMag"),
  volume(baseName+"volume"),
  volumeN1(baseName+"volumeN1"),
  volumeN2(baseName+"volumeN2"),
  sweptVolDot(baseName+"sweptVolDot"),
  sweptVolDotN1(baseName+"sweptVolDotN1"),
  gridFlux(baseName+"gridFlux"),
  nodeDisplacement(baseName+"nodeDisplacement"),
  boundaryNodeNormal(baseName+"boundaryNodeNormal")
{}

