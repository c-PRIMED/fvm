// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

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
  _ghostArraysLevel1.clear();
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
Field::createSyncGatherArraysLevel1(const StorageSite& site)
{
  ArrayBase& thisArray = operator[](site);

  const StorageSite::GatherMap& gatherMap = site.getGatherMapLevel1();

  foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap)
  {
      const StorageSite& oSite = *mpos.first;
      EntryIndex e(&oSite, &site);
      const Array<int>& toIndices = *(mpos.second);
      if (_ghostArraysLevel1.find(e) == _ghostArraysLevel1.end())
      {
          _ghostArraysLevel1[e] = thisArray.newSizedClone(toIndices.getLength());
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
Field::syncScatterLevel1(const StorageSite& site)
{
  const ArrayBase& thisArray = operator[](site);
      	
  const StorageSite::ScatterMap& scatterMap = site.getScatterMapLevel1();
  
  foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap)
  {
      const StorageSite& oSite = *mpos.first;
      EntryIndex e(&site, &oSite);	
      const Array<int>& fromIndices = *(mpos.second);
      if (_ghostArraysLevel1.find(e) == _ghostArraysLevel1.end()){
        _ghostArraysLevel1[e] = thisArray.newSizedClone( fromIndices.getLength() );
      }

      ArrayBase& ghostArray = *_ghostArraysLevel1[e];	
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
Field::syncGatherLevel1(const StorageSite& site)
{
  ArrayBase& thisArray = operator[](site);
      
  const StorageSite::GatherMap& gatherMap = site.getGatherMapLevel1();
  
  foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap)
  {
      const StorageSite& oSite = *mpos.first;
      const Array<int>& toIndices = *(mpos.second);
      EntryIndex e(&oSite, &site);

      if (_ghostArraysLevel1.find(e) == _ghostArraysLevel1.end())
      {
         ostringstream e;
         e << "Field::syncScatter: ghost array not found for"
           << &oSite << endl;
         throw CException(e.str());
      }

      const ArrayBase& ghostArray = *_ghostArraysLevel1[e];
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


   //sycnLocal1
   syncLocalLevel1();

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



void
Field::syncLocalLevel1()
{  

   // scatter first (prepare ship packages)
   foreach(ArrayMap::value_type& pos, _arrays){
      if ( pos.second->getLength() == pos.first->getCountLevel1() )
         syncScatterLevel1(*pos.first);
  }

   foreach(ArrayMap::value_type& pos, _arrays){
      if ( pos.second->getLength() == pos.first->getCountLevel1() )
         createSyncGatherArraysLevel1(*pos.first);
   }

#ifdef FVM_PARALLEL
   //SENDING
   int countScatter   = get_request_size_scatter_level1();
   int countGather    = get_request_size_gather_level1();
   MPI::Request   request_send[ countScatter ];
   MPI::Request   request_recv[ countGather  ];
   int indxSend = 0;
   int indxRecv = 0;
   foreach(ArrayMap::value_type& pos, _arrays){
      const StorageSite& site = *pos.first;
      if ( pos.second->getLength() == pos.first->getCountLevel1() ){
         const StorageSite::ScatterMap& scatterMap = site.getScatterMapLevel1();
         foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
             const StorageSite&  oSite = *mpos.first;
             //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
             EntryIndex e(&site,&oSite);
             ArrayBase& sendArray = *_ghostArraysLevel1[e];
             int to_where  = oSite.getGatherProcID();
             if ( to_where != -1 ){
                int mpi_tag = oSite.getTag();
                request_send[indxSend++] =  
                     MPI::COMM_WORLD.Isend( sendArray.getData(), sendArray.getDataSize(), MPI::BYTE, to_where, mpi_tag );
             }
         }
         //RECIEVING
         //getting values from other meshes to fill g
         const StorageSite::GatherMap& gatherMap = site.getGatherMapLevel1();
         foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap){
            const StorageSite&  oSite = *mpos.first;
            //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
            EntryIndex e(&oSite,&site);
            ArrayBase& recvArray = *_ghostArraysLevel1[e];
            int from_where       = oSite.getGatherProcID();
            if ( from_where != -1 ){
               int mpi_tag = oSite.getTag();
               request_recv[indxRecv++] =  
                   MPI::COMM_WORLD.Irecv( recvArray.getData(), recvArray.getDataSize(), MPI::BYTE, from_where, mpi_tag );
            }
         }

      }
   }

   MPI::Request::Waitall( countGather, request_recv );
   MPI::Request::Waitall( countScatter, request_send );
#endif

   // gather 
   foreach(ArrayMap::value_type& pos, _arrays){
      if ( pos.second->getLength() == pos.first->getCountLevel1() )
        syncGatherLevel1(*pos.first);
   }
}

int
Field::get_request_size_scatter_level1()
{
   int indx =  0;
   foreach(ArrayMap::value_type& pos, _arrays){
      const StorageSite& site = *pos.first;
      const StorageSite::ScatterMap& scatterMap = site.getScatterMapLevel1();
      if ( pos.second->getLength() == pos.first->getCountLevel1() ){
         foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
            const StorageSite&  oSite = *mpos.first;
            //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
            if ( oSite.getGatherProcID() != -1 )
               indx++;
         }
      }
    }
   return indx;

}

int
Field::get_request_size_gather_level1()
{
   int indx =  0;
   foreach(ArrayMap::value_type& pos, _arrays){
      const StorageSite& site = *pos.first;
      const StorageSite::GatherMap& gatherMap = site.getGatherMapLevel1();
      if ( pos.second->getLength() == pos.first->getCountLevel1() ){
         foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap){
            const StorageSite&  oSite = *mpos.first;
            //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
            if ( oSite.getGatherProcID() != -1 )
               indx++;
         }
      }
   }
   return indx;

}


//only for one field need to be communicated
void
Field::createSyncGatherArraysVectorFields(const StorageSite& site, Field& field, const size_t numDir)
{
  ArrayBase& thisArray = field[site];
  GhostArrayMap& ghostArrays = field.getGhostArrayMap();

  const StorageSite::GatherMap& gatherMap = site.getGatherMap();

  foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap)
  {
      const StorageSite& oSite = *mpos.first;
      EntryIndex e(&oSite, &site);
      const Array<int>& toIndices = *(mpos.second);
      if (ghostArrays.find(e) == ghostArrays.end())
      {
          ghostArrays[e] = thisArray.newSizedClone(toIndices.getLength() * numDir);
      }	
  }
}

void
Field::createSyncGatherArraysVectorFieldsLevel1(const StorageSite& site, Field& field, const size_t numDir)
{
  ArrayBase& thisArray = field[site];
  GhostArrayMap& ghostArrays = field.getGhostArrayMapLevel1();

  const StorageSite::GatherMap& gatherMap = site.getGatherMapLevel1();

  foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap)
  {
      const StorageSite& oSite = *mpos.first;
      EntryIndex e(&oSite, &site);
      const Array<int>& toIndices = *(mpos.second);
      if (ghostArrays.find(e) == ghostArrays.end())
      {
          ghostArrays[e] = thisArray.newSizedClone(toIndices.getLength() * numDir);
      }	
  }
}


void
Field::syncScatterVectorFields(const StorageSite& site, std::vector<Field*> & dsf)
{
   
  const size_t numDir = dsf.size();
  Field& field0 = *dsf[0];
  ArrayBase& thisArray0 = field0[site];
  GhostArrayMap& ghostArrays = field0.getGhostArrayMap();
      
  const StorageSite::ScatterMap& scatterMap = site.getScatterMap();
  
  foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap)
  {
      const StorageSite& oSite = *mpos.first;
      EntryIndex e(&site, &oSite);	
      const Array<int>& fromIndices = *(mpos.second);
      if (ghostArrays.find(e) == ghostArrays.end()){
        ghostArrays[e] = thisArray0.newSizedClone( fromIndices.getLength() * numDir );
      }

      const int chunkSize = fromIndices.getLength();	 
      ArrayBase& ghostArray = *ghostArrays[e];	
      for (size_t dir=0; dir < numDir; dir++){
	Field& field = *dsf[dir];
	const ArrayBase& thisArray = field[site];
	const int offset = chunkSize * dir;
        thisArray.scatter(ghostArray,fromIndices, offset);
      } 
  }
}

void
Field::syncScatterVectorFieldsLevel1(const StorageSite& site, std::vector<Field*> & dsf)
{
   
  const size_t numDir = dsf.size();
  Field& field0 = *dsf[0];
  ArrayBase& thisArray0 = field0[site];
  GhostArrayMap& ghostArrays = field0.getGhostArrayMapLevel1();
      
  const StorageSite::ScatterMap& scatterMap = site.getScatterMapLevel1();
  
  foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap)
  {
      const StorageSite& oSite = *mpos.first;
      EntryIndex e(&site, &oSite);	
      const Array<int>& fromIndices = *(mpos.second);
      if (ghostArrays.find(e) == ghostArrays.end()){
        ghostArrays[e] = thisArray0.newSizedClone( fromIndices.getLength() * numDir );
      }

      const int chunkSize = fromIndices.getLength();	 
      ArrayBase& ghostArray = *ghostArrays[e];	
      for (size_t dir=0; dir < numDir; dir++){
	Field& field = *dsf[dir];
	const ArrayBase& thisArray = field[site];
	const int offset = chunkSize * dir;
        thisArray.scatter(ghostArray,fromIndices, offset);
      } 
  }
}



void
Field::syncGatherVectorFields(const StorageSite& site,  std::vector<Field*>& dsf)
{
  
  const size_t  numDir = dsf.size();
  Field& field0 = *dsf[0];
  GhostArrayMap& ghostArrays = field0.getGhostArrayMap();
     
  const StorageSite::GatherMap& gatherMap = site.getGatherMap();
  
  foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap)
  {
      const StorageSite& oSite = *mpos.first;
      const Array<int>& toIndices = *(mpos.second);
      EntryIndex e(&oSite, &site);

      if (ghostArrays.find(e) == ghostArrays.end())
      {
         ostringstream e;
         e << "Field::syncScatter: ghost array not found for"
           << &oSite << endl;
         throw CException(e.str());
      }

      const ArrayBase& ghostArray = *ghostArrays[e];
      const int chunkSize = toIndices.getLength();
      for (size_t dir = 0; dir < numDir; dir++){
	 Field& field = *dsf[dir];
	 ArrayBase& thisArray = field[site];
	 const size_t offset = dir * chunkSize;
	 thisArray.gather(ghostArray,toIndices, offset);
     }
  }
}

void
Field::syncGatherVectorFieldsLevel1(const StorageSite& site,  std::vector<Field*>& dsf)
{
  
  const size_t  numDir = dsf.size();
  Field& field0 = *dsf[0];
  GhostArrayMap& ghostArrays = field0.getGhostArrayMapLevel1();
     
  const StorageSite::GatherMap& gatherMap = site.getGatherMapLevel1();
  
  foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap)
  {
      const StorageSite& oSite = *mpos.first;
      const Array<int>& toIndices = *(mpos.second);
      EntryIndex e(&oSite, &site);

      if (ghostArrays.find(e) == ghostArrays.end())
      {
         ostringstream e;
         e << "Field::syncScatter: ghost array not found for"
           << &oSite << endl;
         throw CException(e.str());
      }

      const ArrayBase& ghostArray = *ghostArrays[e];
      const int chunkSize = toIndices.getLength();
      for (size_t dir = 0; dir < numDir; dir++){
	 Field& field = *dsf[dir];
	 ArrayBase& thisArray = field[site];
	 const size_t offset = dir * chunkSize;
	 thisArray.gather(ghostArray,toIndices, offset);
     }
  }
}
void 
Field::syncLocalVectorFields(std::vector<Field*> & dsf)
{
       // scatter first (prepare ship packages)
   Field& field0 = *dsf[0];
   ArrayMap& arrays = field0.getArrayMap();
   foreach(ArrayMap::value_type& pos, arrays)
      Field::syncScatterVectorFields(*pos.first, dsf);

   const size_t numDir = dsf.size();
   foreach(ArrayMap::value_type& pos, arrays)
      Field::createSyncGatherArraysVectorFields(*pos.first, field0, numDir );
   
   GhostArrayMap& ghostArrays = field0.getGhostArrayMap();
   
#ifdef FVM_PARALLEL
   //SENDING
   MPI::Request   request_send[ Field::get_request_size(field0) ];
   MPI::Request   request_recv[ Field::get_request_size(field0) ];
   int indxSend = 0;
   int indxRecv = 0;
   foreach(ArrayMap::value_type& pos, arrays){
      const StorageSite& site = *pos.first;
      const StorageSite::ScatterMap& scatterMap = site.getScatterMap();
      foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
          const StorageSite&  oSite = *mpos.first;
          //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
             EntryIndex e(&site,&oSite);
             ArrayBase& sendArray = *ghostArrays[e];
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
         ArrayBase& recvArray = *ghostArrays[e];
         int from_where       = oSite.getGatherProcID();
         if ( from_where != -1 ){
             int mpi_tag = oSite.getTag();
             request_recv[indxRecv++] =  
                    MPI::COMM_WORLD.Irecv( recvArray.getData(), recvArray.getDataSize(), MPI::BYTE, from_where, mpi_tag );
         }
      }
   }

   int count  = Field::get_request_size(field0);
   MPI::Request::Waitall( count, request_recv );
   MPI::Request::Waitall( count, request_send );
#endif

  // gather 
  foreach(ArrayMap::value_type& pos, arrays)
    Field::syncGatherVectorFields(*pos.first, dsf);


   //sycnLocal1
   
  Field::syncLocalVectorFieldsLevel1(dsf);
} 
 
void 
Field::syncLocalVectorFieldsLevel1(std::vector<Field*> & dsf)
{
   // scatter first (prepare ship packages)
   Field& field0 = *dsf[0];
   ArrayMap& arrays = field0.getArrayMap();
   foreach(ArrayMap::value_type& pos, arrays){ 
      if ( pos.second->getLength() == pos.first->getCountLevel1() )
          Field::syncScatterVectorFieldsLevel1(*pos.first, dsf);
   }
   
   const size_t numDir = dsf.size();
   foreach(ArrayMap::value_type& pos, arrays){
      if ( pos.second->getLength() == pos.first->getCountLevel1() )   
         Field::createSyncGatherArraysVectorFieldsLevel1(*pos.first, field0, numDir );
   }
     
   GhostArrayMap& ghostArrays = field0.getGhostArrayMapLevel1();
   
#ifdef FVM_PARALLEL
   //SENDING
   MPI::Request   request_send[ Field::get_request_size_level1(field0) ];
   MPI::Request   request_recv[ Field::get_request_size_level1(field0) ];
   int indxSend = 0;
   int indxRecv = 0;
   foreach(ArrayMap::value_type& pos, arrays){
      const StorageSite& site = *pos.first;
      const StorageSite::ScatterMap& scatterMap = site.getScatterMapLevel1();
      foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
          const StorageSite&  oSite = *mpos.first;
          //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
             EntryIndex e(&site,&oSite);
             ArrayBase& sendArray = *ghostArrays[e];
             int to_where  = oSite.getGatherProcID();
             if ( to_where != -1 ){
                int mpi_tag = oSite.getTag();
                request_send[indxSend++] =  
                     MPI::COMM_WORLD.Isend( sendArray.getData(), sendArray.getDataSize(), MPI::BYTE, to_where, mpi_tag );
             }
      }
      //RECIEVING
      //getting values from other meshes to fill g
      const StorageSite::GatherMap& gatherMap = site.getGatherMapLevel1();
      foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap){
         const StorageSite&  oSite = *mpos.first;
         //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
         EntryIndex e(&oSite,&site);
         ArrayBase& recvArray = *ghostArrays[e];
         int from_where       = oSite.getGatherProcID();
         if ( from_where != -1 ){
             int mpi_tag = oSite.getTag();
             request_recv[indxRecv++] =  
                    MPI::COMM_WORLD.Irecv( recvArray.getData(), recvArray.getDataSize(), MPI::BYTE, from_where, mpi_tag );
         }
      }
   }

   int count  = Field::get_request_size_level1(field0);
   MPI::Request::Waitall( count, request_recv );
   MPI::Request::Waitall( count, request_send );
#endif

  // gather 
  foreach(ArrayMap::value_type& pos, arrays){
     if ( pos.second->getLength() == pos.first->getCountLevel1() )
         Field::syncGatherVectorFieldsLevel1(*pos.first, dsf);
  }	 


  
}  


int
Field::get_request_size(Field& field)
{
   ArrayMap& arrays = field.getArrayMap();
   int indx =  0;
   foreach(ArrayMap::value_type& pos, arrays){
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


int
Field::get_request_size_level1(Field& field)
{
   ArrayMap& arrays = field.getArrayMap();
   int indx =  0;
   foreach(ArrayMap::value_type& pos, arrays){
      const StorageSite& site = *pos.first;
      const StorageSite::ScatterMap& scatterMap = site.getScatterMapLevel1();
      foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
          const StorageSite&  oSite = *mpos.first;
          //checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
          if ( oSite.getGatherProcID() != -1 )
             indx++;
      }
   }
   return indx;

}

