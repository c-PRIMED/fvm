%{
#include "CException.h"
#include <rlog/StdioNode.h>
#include <rlog/RLogChannel.h>
#include "ArrayBase.h"
#include "Field.h"
#include "Mesh.h"
  %}

class Field
{
public:
  Field(const string& name);


  void syncLocal();

  const string getName() const;

  void clear();
  
  void removeArray(const StorageSite&);
  
  %extend
  {
    boost::shared_ptr<ArrayBase> __getitem__(const StorageSite* s)
    {
        return self->getArrayPtr(*s);
    }
    void __setitem__(const StorageSite* s, boost::shared_ptr<ArrayBase> a)
    {
      self->addArray(*s,a);
    }
    
    bool __contains__(const StorageSite* s)
    {
        return self->hasArray(*s);
    }
  }
private:
};
