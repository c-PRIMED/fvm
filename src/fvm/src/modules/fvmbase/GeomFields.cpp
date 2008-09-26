#include "GeomFields.h"

GeomFields::GeomFields(const string baseName) :
  coordinate(baseName+"coordinate"),
  area(baseName+"area"),
  areaMag(baseName+"areaMag"),
  volume(baseName+"volume")
{}

