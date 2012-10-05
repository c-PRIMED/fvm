// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

#include "MultiField.h"
#include "OneToOneIndexMap.h"
#include "MultiFieldReduction.h"



#include <iostream>

MultiField::MultiField():
  IContainer(),
  _length(0),
  _arrays(),
  _arrayIndices(),
  _arrayMap()
{
  logCtor();
}

MultiField::~MultiField()
{
  logDtor();
}

bool
MultiField::hasArray(const ArrayIndex& i) const
{
  return _arrayMap.find(i) != _arrayMap.end();
}

const ArrayBase&
MultiField::operator[](const ArrayIndex& i) const
{
  ArrayMap::const_iterator pos = _arrayMap.find(i);
  if (pos != _arrayMap.end())
  {
      return *_arrays[pos->second];
  }

  ostringstream e;
  e << "MultiField::operator[] No array found";
   throw CException(e.str());
}

ArrayBase&
MultiField::operator[](const ArrayIndex& i)
{
  ArrayMap::const_iterator pos = _arrayMap.find(i);
  if (pos != _arrayMap.end())
  {
      return *_arrays[pos->second];
  }

  ostringstream e;
  e << "MultiField::operator[] No array found";
  throw CException(e.str());
}

shared_ptr<ArrayBase>
MultiField::getArrayPtr(const ArrayIndex& i)
{
  ArrayMap::const_iterator pos = _arrayMap.find(i);
  if (pos != _arrayMap.end())
  {
      return _arrays[pos->second];
  }

  ostringstream e;
  e << "MultiField::operator[] No array found";
  throw CException(e.str());
}

shared_ptr<IContainer>
MultiField::newClone() const
{
  shared_ptr<MultiField> c(new MultiField());
  for(int i=0; i<_length; i++)
  {
      const ArrayBase& myArray = *_arrays[i];
      shared_ptr<IContainer> aCloneI(myArray.newClone());
      shared_ptr<ArrayBase> aClone(dynamic_pointer_cast<ArrayBase>(aCloneI));
      if (aClone)
        c->addArray(_arrayIndices[i],aClone);
      else
      {
          throw CException("array clone failed");
      }
  }
  return c;
}


shared_ptr<IContainer>
MultiField::newCopy() const
{
  shared_ptr<MultiField> c(new MultiField());
  for(int i=0; i<_length; i++)
  {
      const ArrayBase& myArray = *_arrays[i];
      shared_ptr<ArrayBase> aCopy(dynamic_pointer_cast<ArrayBase>(myArray.newCopy()));
      c->addArray(_arrayIndices[i],aCopy);
  }
  return c;
}


void
MultiField::zero()
{
  for(int i=0; i<_length; i++)
    _arrays[i]->zero();

  foreach(GhostArrayMap::value_type& pos, _ghostArrays)
  {
      pos.second->zero();
  }
}

void
MultiField::copyFrom(const IContainer& oc)
{
  const MultiField& ofield = dynamic_cast<const MultiField&>(oc);
  foreach(ArrayMap::value_type& pos, _arrayMap)
  {
      ArrayBase& myArray = operator[](pos.second);
      const ArrayBase& oArray = ofield[pos.first];
      myArray.copyFrom(oArray);
  }
}

MultiField&
MultiField::operator+=(const MultiField& ofield)
{
  for(int i=0; i<_length; i++)
  {
      ArrayBase& myArray = *_arrays[i];
      const ArrayBase& oArray = ofield[_arrayIndices[i]];
      myArray += oArray;
  }
  return *this;
}

MultiField&
MultiField::operator-=(const MultiField& ofield)
{
  for(int i=0; i<_length; i++)
  {
      ArrayBase& myArray = *_arrays[i];
      const ArrayBase& oArray = ofield[_arrayIndices[i]];
      myArray -= oArray;
  }
  return *this;
}

MultiField&
MultiField::saxpy(const MultiFieldReduction& alphaMF, const MultiField& xMF)
{
  for(int i=0; i<_length; i++)
  {
      ArrayBase& myArray = *_arrays[i];
      const ArrayBase& alpha = alphaMF[*_arrayIndices[i].first];
      const ArrayBase& x = xMF[_arrayIndices[i]];
      
      myArray.saxpy(alpha,x);
  }
  return *this;
}

MultiField&
MultiField::msaxpy(const MultiFieldReduction& alphaMF, const MultiField& xMF)
{
  for(int i=0; i<_length; i++)
  {
      ArrayBase& myArray = *_arrays[i];
      const ArrayBase& alpha = alphaMF[*_arrayIndices[i].first];
      const ArrayBase& x = xMF[_arrayIndices[i]];
      
      myArray.msaxpy(alpha,x);
  }
  return *this;
}

 
MultiField&
MultiField::operator/=(const MultiFieldReduction& alphaMF)
{
  for(int i=0; i<_length; i++)
  {
      ArrayBase& myArray = *_arrays[i];
      const ArrayBase& alpha = alphaMF[*_arrayIndices[i].first];
      myArray /= alpha;
  }
  return *this;
}

MultiField&
MultiField::operator*=(const MultiFieldReduction& alphaMF)
{
  for(int i=0; i<_length; i++)
  {
      ArrayBase& myArray = *_arrays[i];
      const ArrayBase& alpha = alphaMF[*_arrayIndices[i].first];
      myArray *= alpha;
  }
  return *this;
}

shared_ptr<MultiFieldReduction>
MultiField::getOneNorm() const
{
  shared_ptr<MultiField> c(new MultiField());
  for(int i=0; i<_length; i++)
  {
      const  StorageSite& thisSite = *(_arrayIndices[i].second);
      const ArrayBase& myArray = *_arrays[i];
      c->addArray(_arrayIndices[i],
                  myArray.getOneNorm(thisSite.getSelfCount()));
  }
  shared_ptr<MultiFieldReduction> r(c->reduceSum());
  return r;
}

shared_ptr<MultiFieldReduction>
MultiField::dotWith(const MultiField& ofield) const
{
  shared_ptr<MultiField> dotpField(new MultiField());

  for(int i=0; i<_length; i++)
  {
      const  StorageSite& thisSite = *(_arrayIndices[i].second);
      const ArrayBase& myArray = *_arrays[i];
      const ArrayBase& otherArray = ofield[_arrayIndices[i]];
      dotpField->addArray(_arrayIndices[i],
                          myArray.dotWith(otherArray, thisSite.getSelfCount()));
  }

  shared_ptr<MultiFieldReduction> r(dotpField->reduceSum());
  return r;
}

shared_ptr<MultiFieldReduction>
MultiField::reduceSum() const
{
  shared_ptr<MultiFieldReduction> sum(new MultiFieldReduction());
  
  for(int i=0; i<_length; i++)
  {
      const ArrayBase& myArray = *_arrays[i];
      const Field& fIndex = *_arrayIndices[i].first;
      if (sum->hasArray(fIndex))
        (*sum)[fIndex] += myArray;
      else
        sum->addArray(fIndex,dynamic_pointer_cast<ArrayBase>(myArray.newCopy()));
  }

#ifdef FVM_PARALLEL 
  sum->sync(); //global reduction operation for residual check;
#endif
  return sum;
}

void
MultiField::addArray(const MultiField::ArrayIndex& i, shared_ptr<ArrayBase> a)
{
  _arrayMap[i] = _length;
  _arrays.push_back(a);
  _arrayIndices.push_back(i);
  ++_length;
}


void
MultiField::removeArray(const ArrayIndex& aIndex)
{
  int loc = _arrayMap[aIndex];

  _arrays.erase(_arrays.begin()+loc);
  _arrayIndices.erase(_arrayIndices.begin()+loc);
  _arrayMap.erase(aIndex);

  for(ArrayMap::iterator pos=_arrayMap.begin();
      pos!=_arrayMap.end();
      ++pos)
  {
      int currentValue = pos->second;
      if (currentValue == loc)
        cerr << "Error in removing array" << endl;
      else if (currentValue > loc)
        pos->second = currentValue-1;
  }
  _length = _arrays.size();
}



shared_ptr<MultiField>
MultiField::extract(const ArrayIndexList& indices)
{
  shared_ptr<MultiField> r(new MultiField);
  foreach(ArrayIndex i,indices)
  {
      r->addArray(i,getArrayPtr(i));
      removeArray(i);
  }
  return r;
}

void
MultiField::merge(const MultiField& other)
{
  foreach(const ArrayMap::value_type& pos, other._arrayMap)
  {
      addArray(pos.first, other._arrays[pos.second]);
  }
}

void
MultiField::syncScatter(const ArrayIndex& i)
{

  const  ArrayBase& thisArray  = *_arrays[_arrayMap[i]];
  const  StorageSite& thisSite = *i.second;

  const StorageSite::ScatterMap& scatterMap = thisSite.getScatterMap();

  foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap)
  {
      const StorageSite& oSite = *mpos.first;

      ArrayIndex oIndex(i.first,&oSite);

      // skip if the other site is not in this Multifield
      //if (!hasArray(oIndex))
      //  continue;
      
      // arrays are stored with (Sender,receiver) as the key
      EntryIndex eIndex(i,oIndex);
      const Array<int>& fromIndices = *(mpos.second);
      if (_ghostArrays.find(eIndex) == _ghostArrays.end()){
        _ghostArrays[eIndex] = thisArray.newSizedClone(fromIndices.getLength());
      }

      ArrayBase& ghostArray = *_ghostArrays[eIndex];
      thisArray.scatter(ghostArray,fromIndices);
  }

}


void
MultiField::syncScatterLevel1(const ArrayIndex& i)
{

  const  ArrayBase& thisArray  = *_arrays[_arrayMap[i]];
  const  StorageSite& thisSite = *i.second;
  if ( thisArray.getLength() == thisSite.getCountLevel1() ){
     const StorageSite::ScatterMap& scatterMap = thisSite.getScatterMapLevel1();

     foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap)
     {
        const StorageSite& oSite = *mpos.first;

        ArrayIndex oIndex(i.first,&oSite);
        // arrays are stored with (Sender,receiver) as the key
        EntryIndex eIndex(i,oIndex);
        const Array<int>& fromIndices = *(mpos.second);
        if (_ghostArraysLevel1.find(eIndex) == _ghostArraysLevel1.end()){
          _ghostArraysLevel1[eIndex] = thisArray.newSizedClone(fromIndices.getLength());
        }

        ArrayBase& ghostArray = *_ghostArraysLevel1[eIndex];
        thisArray.scatter(ghostArray,fromIndices);
     }
  }

}

void
MultiField::createSyncGatherArrays(const ArrayIndex& i)
{

  ArrayBase& thisArray = *_arrays[_arrayMap[i]];
  const StorageSite& thisSite = *i.second;

  const StorageSite::GatherMap& gatherMap = thisSite.getGatherMap();

  foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap)
  {
      const StorageSite& oSite = *mpos.first;

      ArrayIndex oIndex(i.first,&oSite);
      EntryIndex eIndex(oIndex,i);

     const Array<int>& toIndices = *(mpos.second);
      if (_ghostArrays.find(eIndex) == _ghostArrays.end()){
         _ghostArrays[eIndex]  = thisArray.newSizedClone(toIndices.getLength());
      }
  }

}

void
MultiField::createSyncGatherArraysLevel1(const ArrayIndex& i)
{

  ArrayBase& thisArray = *_arrays[_arrayMap[i]];
  const StorageSite& thisSite = *i.second;
  
  if ( thisArray.getLength() == thisSite.getCountLevel1() ){
     const StorageSite::GatherMap& gatherMap = thisSite.getGatherMapLevel1();

     foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap)
     {
        const StorageSite& oSite = *mpos.first;

        ArrayIndex oIndex(i.first,&oSite);
        EntryIndex eIndex(oIndex,i);

        const Array<int>& toIndices = *(mpos.second);
        if (_ghostArraysLevel1.find(eIndex) == _ghostArraysLevel1.end()){
          _ghostArraysLevel1[eIndex]  = thisArray.newSizedClone(toIndices.getLength());
        }
     }
  }

}

void
MultiField::syncGather(const ArrayIndex& i)
{

  ArrayBase& thisArray = *_arrays[_arrayMap[i]];
  const StorageSite& thisSite = *i.second;

  const StorageSite::GatherMap& gatherMap = thisSite.getGatherMap();
 
  foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap)
  {
      const StorageSite& oSite = *mpos.first;
      ArrayIndex oIndex(i.first,&oSite);

      const Array<int>& toIndices = *(mpos.second);
      EntryIndex eIndex(oIndex,i);

      if (_ghostArrays.find(eIndex) != _ghostArrays.end()){	
          const ArrayBase& ghostArray = *_ghostArrays[eIndex];
          thisArray.gather(ghostArray,toIndices);
      }
  }

}

void
MultiField::syncGatherLevel1(const ArrayIndex& i)
{

  ArrayBase& thisArray = *_arrays[_arrayMap[i]];
  const StorageSite& thisSite = *i.second;

  if ( thisArray.getLength() == thisSite.getCountLevel1() ){
     const StorageSite::GatherMap& gatherMap = thisSite.getGatherMapLevel1();
 
     foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap)
     {
        const StorageSite& oSite = *mpos.first;
        ArrayIndex oIndex(i.first,&oSite);

        const Array<int>& toIndices = *(mpos.second);
        EntryIndex eIndex(oIndex,i);

        if (_ghostArraysLevel1.find(eIndex) != _ghostArraysLevel1.end()){	
           const ArrayBase& ghostArray = *_ghostArraysLevel1[eIndex];
           thisArray.gather(ghostArray,toIndices);
        }  
     }
  }

}


void
MultiField::sync()
{
  foreach(ArrayIndex i, _arrayIndices)
    syncScatter(i);

  foreach(ArrayIndex i, _arrayIndices)
        createSyncGatherArrays(i);

#ifdef FVM_PARALLEL
  //this communication is based on assumption that MPI communication will be only done inside the same field
  //SENDING
   MPI::Request   request_send[ get_request_size() ];
   MPI::Request   request_recv[ get_request_size() ];
   int sendIndx = 0;
   int recvIndx = 0;
   foreach ( ArrayIndex i, _arrayIndices ){
       const  StorageSite& thisSite = *i.second;
       const StorageSite::ScatterMap& scatterMap = thisSite.getScatterMap();
       foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
           const StorageSite& oSite = *mpos.first;
           ArrayIndex oIndex(i.first,&oSite);
          // arrays are stored with (Sender,receiver) as the key
             EntryIndex eIndex(i,oIndex);
             ArrayBase& sendArray = *_ghostArrays[eIndex];
             int to_where  = oSite.getGatherProcID();
             if ( to_where != -1 ){
                int mpi_tag = oSite.getTag();
                request_send[sendIndx++] =  
                       MPI::COMM_WORLD.Isend( sendArray.getData(), sendArray.getDataSize(), MPI::BYTE, to_where, mpi_tag );
             }
       }

  
     //RECIEVING
      const StorageSite::GatherMap& gatherMap = thisSite.getGatherMap();
      foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap){
           const StorageSite& oSite = *mpos.first;
           ArrayIndex oIndex(i.first,&oSite);
             EntryIndex eIndex(oIndex,i);
             ArrayBase& recvArray = *_ghostArrays[eIndex];
             int from_where  = oSite.getGatherProcID();
             if ( from_where != -1 ){
                 int mpi_tag = oSite.getTag();
                 request_recv[recvIndx++] =  
                       MPI::COMM_WORLD.Irecv( recvArray.getData(), recvArray.getDataSize(), MPI::BYTE, from_where, mpi_tag );
             }
      }

   }

   int count   = get_request_size();
   MPI::Request::Waitall( count, request_recv );
   MPI::Request::Waitall( count, request_send );

#endif

  foreach(ArrayIndex i, _arrayIndices)
    syncGather(i);

  //syncLevel1
  syncLevel1();

}


int
MultiField::get_request_size()
{
   int indx = 0;
   foreach ( ArrayIndex i, _arrayIndices ){
       const  StorageSite& thisSite = *i.second;
       const StorageSite::ScatterMap& scatterMap = thisSite.getScatterMap();
       foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
           const StorageSite& oSite = *mpos.first;
           ArrayIndex oIndex(i.first,&oSite);
          // arrays are stored with (Sender,receiver) as the key
           if ( oSite.getGatherProcID() != -1 )
              indx++;
       }
    }
   return indx;

}


void
MultiField::syncLevel1()
{
  foreach(ArrayIndex i, _arrayIndices)
    syncScatterLevel1(i);

  foreach(ArrayIndex i, _arrayIndices)
        createSyncGatherArraysLevel1(i);

#ifdef FVM_PARALLEL
  //this communication is based on assumption that MPI communication will be only done inside the same field
  //SENDING
   int countScatter   = get_request_size_scatter_level1();
   int countGather    = get_request_size_gather_level1();
   MPI::Request   request_send[ countScatter ];
   MPI::Request   request_recv[ countGather  ];
   int sendIndx = 0;
   int recvIndx = 0;
   foreach ( ArrayIndex i, _arrayIndices ){
       const  StorageSite& thisSite = *i.second;
       const StorageSite::ScatterMap& scatterMap = thisSite.getScatterMapLevel1();
       const ArrayBase& thisArray = *_arrays[_arrayMap[i]];
       if ( thisArray.getLength() == thisSite.getCountLevel1() ){
          foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
             const StorageSite& oSite = *mpos.first;
             ArrayIndex oIndex(i.first,&oSite);
             // arrays are stored with (Sender,receiver) as the key
             EntryIndex eIndex(i,oIndex);
             ArrayBase& sendArray = *_ghostArraysLevel1[eIndex];
             int to_where  = oSite.getGatherProcID();
             if ( to_where != -1 ){
                int mpi_tag = oSite.getTag();
                request_send[sendIndx++] =  
                       MPI::COMM_WORLD.Isend( sendArray.getData(), sendArray.getDataSize(), MPI::BYTE, to_where, mpi_tag );
             }
          }

         //RECIEVING
         const StorageSite::GatherMap& gatherMap = thisSite.getGatherMapLevel1();
         foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap){
            const StorageSite& oSite = *mpos.first;
            ArrayIndex oIndex(i.first,&oSite);
            EntryIndex eIndex(oIndex,i);
            ArrayBase& recvArray = *_ghostArraysLevel1[eIndex];
            int from_where  = oSite.getGatherProcID();
            if ( from_where != -1 ){
               int mpi_tag = oSite.getTag();
               request_recv[recvIndx++] =  
                     MPI::COMM_WORLD.Irecv( recvArray.getData(), recvArray.getDataSize(), MPI::BYTE, from_where, mpi_tag );
           }
         }
 
     }
   }

   MPI::Request::Waitall( countGather, request_recv );
   MPI::Request::Waitall( countScatter , request_send );

#endif

  foreach(ArrayIndex i, _arrayIndices)
    syncGatherLevel1(i);

}


int
MultiField::get_request_size_scatter_level1()
{
   int indx = 0;
   foreach ( ArrayIndex i, _arrayIndices ){
       const  StorageSite& thisSite = *i.second;
       const StorageSite::ScatterMap& scatterMap = thisSite.getScatterMapLevel1();
       const ArrayBase& thisArray = *_arrays[_arrayMap[i]];
       if ( thisArray.getLength() == thisSite.getCountLevel1() ){
          foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
             const StorageSite& oSite = *mpos.first;
             ArrayIndex oIndex(i.first,&oSite);
             // arrays are stored with (Sender,receiver) as the key
             if ( oSite.getGatherProcID() != -1 )
               indx++;
          }
       }
    }
   return indx;

}

int
MultiField::get_request_size_gather_level1()
{
   int indx = 0;
   foreach ( ArrayIndex i, _arrayIndices ){
       const  StorageSite& thisSite = *i.second;
       const StorageSite::GatherMap& gatherMap = thisSite.getGatherMapLevel1();
       const ArrayBase& thisArray = *_arrays[_arrayMap[i]];
       if ( thisArray.getLength() == thisSite.getCountLevel1() ){
          foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap){
             const StorageSite& oSite = *mpos.first;
             ArrayIndex oIndex(i.first,&oSite);
             // arrays are stored with (Sender,receiver) as the key
             if ( oSite.getGatherProcID() != -1 )
                indx++;
          }
       }
    }
   return indx;

}
