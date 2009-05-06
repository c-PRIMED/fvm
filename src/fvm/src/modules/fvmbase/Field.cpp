#include "Field.h"
#include "ArrayBase.h"
#include "StorageSite.h"
#include "Array.h"
#include "OneToOneIndexMap.h"
#include <iostream>

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif


#include <vector>
#include "Vector.h"

using namespace std;

Field::Field(const string& name):
  IContainer(),
  _name(name),
  _arrays(),
  MPI_FIELD_TAG(2009),
  _syncGatherArrays(false)
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
  _syncGatherArrays = true;
  ArrayBase& thisArray = operator[](site);

  const StorageSite::GatherMap& gatherMap = site.getGatherMap();

  foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap)
  {
      const StorageSite& oSite = *mpos.first;

      if (_ghostGatherArrays.find(&oSite) == _ghostGatherArrays.end())
      {
          _ghostGatherArrays[&oSite] = thisArray.newSizedClone(oSite.getCount());
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

      const Array<int>& fromIndices = *(mpos.second);
      if (_ghostScatterArrays.find(&oSite) == _ghostScatterArrays.end()){
        _ghostScatterArrays[&oSite] = thisArray.newSizedClone(oSite.getCount());
      }

      ArrayBase& ghostScatterArray = *_ghostScatterArrays[&oSite];
      thisArray.scatter(ghostScatterArray,fromIndices);


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

      if (_ghostGatherArrays.find(&oSite) == _ghostGatherArrays.end())
      {
          _ghostGatherArrays[&oSite] = thisArray.newSizedClone(oSite.getCount());
//          ostringstream e;
//          e << "Field::syncScatter: ghost array not found for"
//            << &oSite << endl;
//          throw CException(e.str());
      }
      
      const ArrayBase& ghostGatherArray = *_ghostGatherArrays[&oSite];
      thisArray.gather(ghostGatherArray,toIndices);
  }
}

void
Field::syncLocal()
{
   // scatter first (prepare ship packages)
   foreach(ArrayMap::value_type& pos, _arrays)
      syncScatter(*pos.first);

  //gather arrays  are allocated (once)
   if ( !_syncGatherArrays )
      foreach(ArrayMap::value_type& pos, _arrays)
         createSyncGatherArrays(*pos.first);

#ifdef FVM_PARALLEL
   //SENDING
   MPI::Request   request_send[ _ghostScatterArrays.size() ];
   int indx = 0;
   foreach( const ArrayMap::value_type& mpos, _ghostScatterArrays){
       const StorageSite&  site = *mpos.first;
       ArrayBase& sendArray = *mpos.second;
       int to_where  = site.getGatherProcID();
       request_send[indx++] =  
             MPI::COMM_WORLD.Isend( sendArray.getData(), sendArray.getDataSize(), MPI::BYTE, to_where, MPI_FIELD_TAG);

   }

   //RECIEVING
   MPI::Request   request_recv[ _ghostGatherArrays.size() ];
   //getting values from other meshes to fill g
   indx = 0;
   foreach( const ArrayMap::value_type& mpos, _ghostGatherArrays){
       const StorageSite&  site = *mpos.first;
       ArrayBase& recvArray = *mpos.second;
       int from_where  = site.getGatherProcID();
       request_recv[indx++] = 
             MPI::COMM_WORLD.Irecv( recvArray.getData(), recvArray.getDataSize(), MPI::BYTE, from_where, MPI_FIELD_TAG );

   }

   int count_recv = _ghostGatherArrays.size();
   MPI::Request::Waitall( count_recv, request_recv);

#endif


  // gather 
  foreach(ArrayMap::value_type& pos, _arrays)
    syncGather(*pos.first);
}