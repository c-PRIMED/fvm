#ifndef _FIELD_H_
#define _FIELD_H_

#include "IContainer.h"
#include "StorageSite.h"
#include "ArrayBase.h"



class Field : public IContainer
{
public:
  typedef map<const StorageSite*, shared_ptr<ArrayBase> > ArrayMap;
  typedef map<const StorageSite*, vector<const StorageSite*>* > ChildSitesMap;
  typedef pair<const StorageSite*, const StorageSite*> EntryIndex;
  typedef map<EntryIndex, shared_ptr<ArrayBase> > GhostArrayMap;

  Field(const string& name);
  
  virtual ~Field();

  DEFINE_TYPENAME("Field");

  Field& operator=(const Field& oField);
  
  const ArrayBase& operator[](const StorageSite&) const;
  ArrayBase& operator[](const StorageSite&);

  shared_ptr<ArrayBase> getArrayPtr(const StorageSite&);
  void addArray(const StorageSite&, shared_ptr<ArrayBase> a);
  void removeArray(const StorageSite&);
  void removeArrays(const StorageSiteList& sites);
  
  virtual void copyFrom(const IContainer& a);
  virtual void zero();
  
  virtual shared_ptr<IContainer> newCopy() const;
  virtual shared_ptr<IContainer> newClone() const;

  bool hasArray(const StorageSite& s) const;
  
  void syncLocal();
  void syncLocalLevel1();

  const string getName() const {return _name;}

  void clear();
  
private:

  Field(const Field&);
  void  createSyncGatherArrays(const StorageSite& site);
  void  syncScatter(const StorageSite& site);
  void  syncGather(const StorageSite& site);
  void  createSyncGatherArraysLevel1(const StorageSite& site);
  void  syncScatterLevel1(const StorageSite& site);
  void  syncGatherLevel1(const StorageSite& site);
  int  get_request_size();
  int  get_request_size_level1();


  //ArrayBase& getGhostArray(const StorageSite&);
  
  const string _name;
  ArrayMap _arrays;
  GhostArrayMap _ghostArrays;
  GhostArrayMap _ghostArraysLevel1;
  
  ChildSitesMap _childSitesMap;

  ArrayBase& _create(const StorageSite& site);
  
};

#endif
