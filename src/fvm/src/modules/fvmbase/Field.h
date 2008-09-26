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
  
  //  void syncGather(const StorageSite& s);
  //void syncScatter(const StorageSite& s);

  const string& getName() const {return _name;}
  
private:
  Field(const Field&);
  const string _name;
  ArrayMap _arrays;
  ChildSitesMap _childSitesMap;

  ArrayBase& _create(const StorageSite& site);
};

#endif
