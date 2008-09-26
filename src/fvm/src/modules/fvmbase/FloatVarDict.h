#ifndef _FLOATVARDICT_H_
#define _FLOATVARDICT_H_

#include <map>
#include "misc.h"

template<class T>
class FloatVarDict : public map<string, T>
{
public:

  typedef map<string,T> T_Parent;
  
  T operator[](const string varName) const
  {
    typename T_Parent::const_iterator pos = this->find(varName);
    if (pos != this->end())
      return pos->second;
    throw CException("uknown var" + varName);
  }

protected:
  void defineVar(const string varName, const T defaultValue)
  {
    this->insert(make_pair(varName,defaultValue));
  }
};

#endif
