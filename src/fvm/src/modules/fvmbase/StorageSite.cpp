#include "StorageSite.h"
#include "OneToOneIndexMap.h"


StorageSite::StorageSite(const int selfCount, const int nGhost,
                         const int offset, const StorageSite* parent):
  _count(selfCount+nGhost),
  _selfCount(selfCount),
  _offset(offset),
  _parent(parent),
  _scatterProcID(-1),
  _gatherProcID(-1),
  _tag(-1)
{
  logCtorVerbose("of size %d", _count);
}

StorageSite::~StorageSite()
{
  logDtorVerbose("of size %d" ,_count);
}


void
StorageSite::clearGatherScatterMaps()
{
  _gatherMap.clear();
  _scatterMap.clear();
}
