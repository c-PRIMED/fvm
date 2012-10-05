// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _MULTIFIELD_H_
#define _MULTIFIELD_H_

#include "misc.h"
#include "IContainer.h"

#include "ArrayBase.h"
#include "StorageSite.h"

class MultiFieldReduction;
class Field;

class MultiField : public IContainer
{
public:

  typedef pair<const Field*,  const StorageSite*> ArrayIndex;
  typedef map<ArrayIndex,int> ArrayMap;
  typedef vector<shared_ptr<ArrayBase> > ArrayList;
  typedef vector<ArrayIndex> ArrayIndexList;
  typedef pair<ArrayIndex,ArrayIndex>  EntryIndex;
  typedef map <EntryIndex,shared_ptr<ArrayBase> > GhostArrayMap;

  MultiField();
  
  virtual ~MultiField();
  DEFINE_TYPENAME("MultiField");

  const ArrayBase& operator[](const ArrayIndex&) const;
  ArrayBase& operator[](const ArrayIndex&);

  shared_ptr<ArrayBase> getArrayPtr(const ArrayIndex&);

  const ArrayBase& operator[](const int i) const {return *_arrays[i];}
  ArrayBase& operator[](const int i) {return *_arrays[i];}

  bool hasArray(const ArrayIndex&) const;
  
  virtual void zero();
  virtual void copyFrom(const IContainer& oc);

  virtual MultiField& operator+=(const MultiField& o);
  virtual MultiField& operator-=(const MultiField& o);
  virtual MultiField& operator/=(const MultiFieldReduction& alpha);
  virtual MultiField& operator*=(const MultiFieldReduction& alpha);

  virtual shared_ptr<IContainer> newCopy() const;
  virtual shared_ptr<IContainer> newClone() const;

  int getLength() const {return _length;}
  const ArrayIndex getArrayIndex(const int i) const {return _arrayIndices[i];}

  void addArray(const ArrayIndex& aIndex, shared_ptr<ArrayBase> a);
  void removeArray(const ArrayIndex& aIndex);

  MultiField& saxpy(const MultiFieldReduction& alphaMF, const MultiField& xMF);
  MultiField& msaxpy(const MultiFieldReduction& alphaMF, const MultiField& xMF);

  shared_ptr<MultiFieldReduction> reduceSum() const;
  shared_ptr<MultiFieldReduction> getOneNorm() const;

  shared_ptr<MultiFieldReduction> dotWith(const MultiField& ofield) const;

  const ArrayIndexList& getArrayIndices() const {return _arrayIndices;}

  shared_ptr<MultiField> extract(const ArrayIndexList& indices);
  void merge(const MultiField& other);

  void  syncScatter(const ArrayIndex& i);
  void syncGather(const ArrayIndex& i);
  void sync();
  
private:

  void  createSyncGatherArrays(const ArrayIndex& i);
  void  syncScatterLevel1(const ArrayIndex& i);
  void  createSyncGatherArraysLevel1(const ArrayIndex& i);
  void syncGatherLevel1(const ArrayIndex& i);

  int  get_request_size();
  int  get_request_size_scatter_level1();
  int  get_request_size_gather_level1();

  void syncLevel1();


  int _length;
  ArrayList _arrays;
  ArrayIndexList _arrayIndices;
  ArrayMap _arrayMap;
  GhostArrayMap _ghostArrays;
  GhostArrayMap _ghostArraysLevel1;
};

#endif
