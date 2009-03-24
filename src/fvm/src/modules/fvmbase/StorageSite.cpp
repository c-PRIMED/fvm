#include "StorageSite.h"
#include "OneToOneIndexMap.h"


StorageSite::StorageSite(const int selfCount, const int nGhost,
                         const int offset, const StorageSite* parent):
  _count(selfCount+nGhost),
  _selfCount(selfCount),
  _offset(offset),
  _parent(parent),
  _mappers()
{
  logCtorVerbose("of size %d", _count);
}

StorageSite::~StorageSite()
{
  logDtorVerbose("of size %d" ,_count);
}


void
StorageSite::scatterGatherMaps( )
{

       MappersMap::const_iterator it_mapper;
       //loop over interfaces
        for ( it_mapper = _mappers.begin(); it_mapper != _mappers.end(); it_mapper++ ){
            const StorageSite *site = it_mapper->first;

            int size = site->getCount();
            shared_ptr< Array<int> > gather_values ( new Array<int> (size) );
            shared_ptr< Array<int> > scatter_values( new Array<int> (size) );

            const OneToOneIndexMap&  indexMap = *(it_mapper->second);
            for ( int i = 0; i < size; i++){
                 (*gather_values)[i]  = indexMap.getFromIndices()[i];
                 (*scatter_values)[i] = indexMap.getToIndices()[i];
            }
            _scatterMap.insert( pair< const StorageSite*, shared_ptr< Array<int> > > ( site, scatter_values ) );
            _gatherMap.insert ( pair< const StorageSite*, shared_ptr< Array<int> > > ( site, gather_values  ) );
         }


}