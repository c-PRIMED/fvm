#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

#include "Field.h"
#include "ArrayBase.h"
#include "StorageSite.h"
#include "Array.h"
#include "OneToOneIndexMap.h"
#include <iostream>


#include <vector>
#include "Vector.h"

using namespace std;

Field::Field(const string& name):
  IContainer(),
  _name(name),
  _arrays()
{
  logCtor();
}  

Field::~Field()
{
  logDtor();
}


bool
Field::hasArray(const StorageSite& s) const
{
  return ((_arrays.find(&s) != _arrays.end()) ||
          (s.getParent() && (_arrays.find(s.getParent()) != _arrays.end())));
}

const ArrayBase&
Field::operator[](const StorageSite& s) const
{
  ArrayMap::const_iterator pos = _arrays.find(&s);
  if (pos != _arrays.end())
    return *pos->second;
  return const_cast<Field&>(*this)._create(s);
}

ArrayBase&
Field::operator[](const StorageSite& s)
{
  ArrayMap::iterator pos = _arrays.find(&s);
  if (pos != _arrays.end())
    return *pos->second;
  
  return _create(s);
}

shared_ptr<ArrayBase>
Field::getArrayPtr(const StorageSite& s)
{
  if (_arrays.find(&s) == _arrays.end())
    _create(s);
  
  return _arrays[&s];
}

void
Field::addArray(const StorageSite& s, shared_ptr<ArrayBase> a)
{
  removeArray(s);
  _arrays[&s]=a;
}

void
Field::removeArray(const StorageSite& s)
{
  if (_childSitesMap.find(&s) != _childSitesMap.end())
  {
      vector<const StorageSite*>& childSites = *_childSitesMap[&s];
      foreach (const StorageSite* site,childSites)
      {
          removeArray(*site);
      }
      delete &childSites;
      _childSitesMap.erase(&s);
      
  }
  _arrays.erase(&s);
}

void
Field::removeArrays(const StorageSiteList& sites)
{
  foreach(const StorageSite* s, sites)
  {
      removeArray(*s);
  }
}

void
Field::clear()
{
  _arrays.clear();
  _childSitesMap.clear();
  _ghostArrays.clear();
}

ArrayBase& 
Field::_create(const StorageSite& s)
{
  const StorageSite* parentSite = s.getParent();
  if (parentSite)
  {
      ArrayBase& parentArray = operator[](*parentSite);
      shared_ptr<ArrayBase> anew = parentArray.createOffsetArray(s.getOffset(),
                                                                 s.getCount());
      _arrays[&s]=anew; 

      if (_childSitesMap.find(parentSite) == _childSitesMap.end())
        _childSitesMap[parentSite] = new vector<const StorageSite*>();

      vector<const StorageSite*>& childSites = *_childSitesMap[parentSite];
      childSites.push_back(&s);
      return *anew;
  }
  else
  {
      ostringstream e;
      e << "Field::operator[] No array found and no creator defined for "
        << &s << endl;
      throw CException(e.str());
  }
}

shared_ptr<IContainer>
Field::newClone() const
{
  shared_ptr<Field> c(new Field(_name));
  
  foreach(const ArrayMap::value_type& pos, _arrays)
  {
      c->addArray(*(pos.first),
                  dynamic_pointer_cast<ArrayBase>(pos.second->newClone()));
  }
  return c;
}


shared_ptr<IContainer>
Field::newCopy() const
{
  shared_ptr<Field> c(new Field(_name));
  foreach(const ArrayMap::value_type& pos, _arrays)
  {
      c->addArray(*(pos.first),
                  dynamic_pointer_cast<ArrayBase>(pos.second->newCopy()));
  }
  return c;
}


void
Field::zero()
{
  foreach(ArrayMap::value_type& pos, _arrays)
  {
      pos.second->zero();
  }
}

Field&
Field::operator=(const Field& oField)
{
  foreach(ArrayMap::value_type& pos, _arrays)
  {
      *pos.second = oField[*pos.first];
  }
  return *this;
}

void
Field::copyFrom(const IContainer& oc)
{
  const Field& other = dynamic_cast<const Field&>(oc);
  operator=(other);
}


void
Field::createSyncGatherArrays(const StorageSite& site)
{
  ArrayBase& thisArray = operator[](site);

  const StorageSite::GatherMap& gatherMap = site.getGatherMap();

  foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap)
  {
      const StorageSite& oSite = *mpos.first;
      EntryIndex e(&oSite, &site);
      const Array<int>& toIndices = *(mpos.second);
      if (_ghostArrays.find(e) == _ghostArrays.end())
      {
          _ghostArrays[e] = thisArray.newSizedClone(toIndices.getLength());
      }	
  }
}

void
Field::syncScatter(const StorageSite& site)
{
  const ArrayBase& thisArray = operator[](site);
      	
  const StorageSite::ScatterMap& scatterMap = site.getScatterMap();
  
  foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap)
  {
      const StorageSite& oSite = *mpos.first;
      EntryIndex e(&site, &oSite);	
      const Array<int>& fromIndices = *(mpos.second);
      if (_ghostArrays.find(e) == _ghostArrays.end()){
        _ghostArrays[e] = thisArray.newSizedClone( fromIndices.getLength() );
      }

      ArrayBase& ghostArray = *_ghostArrays[e];	
      thisArray.scatter(ghostArray,fromIndices);

  }
}

void
Field::syncGather(const StorageSite& site)
{
  ArrayBase& thisArray = operator[](site);
      
  const StorageSite::GatherMap& gatherMap = site.getGatherMap();
  
  foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap)
  {
      const StorageSite& oSite = *mpos.first;
      const Array<int>& toIndices = *(mpos.second);
      EntryIndex e(&oSite, &site);

      if (_ghostArrays.find(e) == _ghostArrays.end())
      {
         ostringstream e;
         e << "Field::syncScatter: ghost array not found for"
           << &oSite << endl;
         throw CException(e.str());
      }

      const ArrayBase& ghostArray = *_ghostArrays[e];
      thisArray.gather(ghostArray,toIndices);
  }
}

void
Field::syncLocal()
{  
   // scatter first (prepare ship packages)
   foreach(ArrayMap::value_type& pos, _arrays)
      syncScatter(*pos.first);

   foreach(ArrayMap::value_type& pos, _arrays)
      createSyncGatherArrays(*pos.first);

#ifdef FVM_PARALLEL
   //SENDING
   MPI::Request   request_send[ get_request_size() ];
   MPI::Request   request_recv[ get_request_size() ];
   int indxSend = 0;
   int indxRecv = 0;
   foreach(ArrayMap::value_type& pos, _arrays){
      const StorageSite& site = *pos.first;
      const StorageSite::ScatterMap& scatterMap = site.getScatterMap();
      foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
          const StorageSite&  oSite = *mpos.first;
          //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
             EntryIndex e(&site,&oSite);
             ArrayBase& sendArray = *_ghostArrays[e];
             int to_where  = oSite.getGatherProcID();
             if ( to_where != -1 ){
                int mpi_tag = oSite.getTag();
                request_send[indxSend++] =  
                     MPI::COMM_WORLD.Isend( sendArray.getData(), sendArray.getDataSize(), MPI::BYTE, to_where, mpi_tag );
             }
      }
      //RECIEVING
      //getting values from other meshes to fill g
      const StorageSite::GatherMap& gatherMap = site.getGatherMap();
      foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap){
         const StorageSite&  oSite = *mpos.first;
         //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
         EntryIndex e(&oSite,&site);
         ArrayBase& recvArray = *_ghostArrays[e];
         int from_where       = oSite.getGatherProcID();
         if ( from_where != -1 ){
             int mpi_tag = oSite.getTag();
             request_recv[indxRecv++] =  
                    MPI::COMM_WORLD.Irecv( recvArray.getData(), recvArray.getDataSize(), MPI::BYTE, from_where, mpi_tag );
         }
       }
   }

   int count  = get_request_size();
   MPI::Request::Waitall( count, request_recv );
   MPI::Request::Waitall( count, request_send );
#endif

  // gather 
  foreach(ArrayMap::value_type& pos, _arrays)
    syncGather(*pos.first);

}

int
Field::get_request_size()
{
   int indx =  0;
   foreach(ArrayMap::value_type& pos, _arrays){
      const StorageSite& site = *pos.first;
      const StorageSite::ScatterMap& scatterMap = site.getScatterMap();
      foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
          const StorageSite&  oSite = *mpos.first;
          //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
          if ( oSite.getGatherProcID() != -1 )
             indx++;
      }
   }
   return indx;

}
