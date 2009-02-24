%{
#include "CException.h"
#include <rlog/StdioNode.h>
#include <rlog/RLogChannel.h>
#include "ArrayBase.h"
#include "Field.h"
#include "Mesh.h"
  %}

%include "ArrayBase.i"

class Field
{
public:
  Field(const string& name);

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
  }
private:
};
