#ifndef _FIELD_H_
#define _FIELD_H_

#include "IContainer.h"
#include "StorageSite.h"
#include "ArrayBase.h"



class Field : public IContainer
{
public:
  typedef map<const StorageSite*, ArrayBase*> ArrayMap;
  typedef map<const StorageSite*, vector<const StorageSite*>* > ChildSitesMap;
  
  Field();
  
  virtual ~Field();

  const ArrayBase& operator[](const StorageSite&) const;
  ArrayBase& operator[](const StorageSite&);

  void addArray(const StorageSite&, ArrayBase* a);
  void removeArray(const StorageSite&);
  void removeArrays(const StorageSiteList& sites);
  
  virtual IContainer& operator=(const IContainer& a);
  virtual void zero();
  virtual void copyFrom(const IContainer& oc);
  
  virtual shared_ptr<IContainer> newCopy() const;
  virtual shared_ptr<IContainer> newClone() const;

  bool hasArray(const StorageSite& s) const;
  
  void syncGather(const StorageSite& s);
  //void syncScatter(const StorageSite& s);
  
private:
  ArrayMap& _arrays;
  ChildSitesMap _childSitesMap;
};

#endif
