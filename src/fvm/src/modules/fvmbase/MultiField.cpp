#include "MultiField.h"
#include "OneToOneIndexMap.h"
#include "MultiFieldReduction.h"

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif


#include <iostream>

MultiField::MultiField():
  IContainer(),
  _length(0),
  _arrays(),
  _arrayIndices(),
  _arrayMap(),
  MPI_MULTIFIELD_TAG(3009),
  _syncGatherArrays( false )
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
MultiField::dotWith(const MultiField& ofield)
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

      const Array<int>& fromIndices = *(mpos.second);
      if (_ghostScatterArrays.find(oIndex) == _ghostScatterArrays.end()){
        _ghostScatterArrays[oIndex] = thisArray.newSizedClone(oSite.getCount());
      }

      ArrayBase& ghostScatterArray = *_ghostScatterArrays[oIndex];
/*      cout << " proc_id = " << MPI::COMM_WORLD.Get_rank() << " ghostScatterArray.getDataSize() = " <<
            ghostScatterArray.getDataSize() << " fromIndices.getDataSize() =  " << fromIndices.getDataSize() << 
            " thisArray.getDatasize() = " << thisArray.getDataSize() << endl;*/
      thisArray.scatter(ghostScatterArray,fromIndices);
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

      if (_ghostGatherArrays.find(oIndex) == _ghostGatherArrays.end())
         _ghostGatherArrays[oIndex]  = thisArray.newSizedClone(oSite.getCount());
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
      if (_ghostGatherArrays.find(oIndex) == _ghostGatherArrays.end())
         _ghostGatherArrays[oIndex]  = thisArray.newSizedClone(oSite.getCount());
       ArrayBase& ghostGatherArray = *_ghostGatherArrays[oIndex];
       thisArray.gather(ghostGatherArray,toIndices);
  }


}

void
MultiField::sync()
{

  foreach(ArrayIndex i, _arrayIndices)
    syncScatter(i);

  if ( !_syncGatherArrays )
  {
     foreach(ArrayIndex i, _arrayIndices)
        createSyncGatherArrays(i);
  }

#ifdef FVM_PARALLEL
  //SENDING
   MPI::Request   request_send[ _ghostScatterArrays.size() ];
   int indx = 0;
   foreach( const GhostArrayMap::value_type& mpos, _ghostScatterArrays ){
       const ArrayIndex& arrayIndex = mpos.first;
       const StorageSite&  site =  *arrayIndex.second;
       ArrayBase& sendArray = *mpos.second;
       int to_where  = site.getGatherProcID();
       request_send[indx++] =  
             MPI::COMM_WORLD.Isend( sendArray.getData(), sendArray.getDataSize(), MPI::BYTE, to_where, MPI_MULTIFIELD_TAG );
   }


   //RECIEVING
   MPI::Request   request_recv[ _ghostGatherArrays.size() ];
   //getting values from other meshes to fill g
   indx = 0;
   foreach( const GhostArrayMap::value_type& mpos, _ghostGatherArrays){
       const ArrayIndex& arrayIndex = mpos.first;
       const StorageSite&  site = *arrayIndex.second;
       ArrayBase& recvArray = *mpos.second;
       int from_where  = site.getGatherProcID();
       request_recv[indx++] = 
             MPI::COMM_WORLD.Irecv( recvArray.getData(), recvArray.getDataSize(), MPI::BYTE, from_where, MPI_MULTIFIELD_TAG );
   }
   int count_recv = _ghostGatherArrays.size();
   MPI::Request::Waitall( count_recv, request_recv );

#endif

  foreach(ArrayIndex i, _arrayIndices)
    syncGather(i);

}
