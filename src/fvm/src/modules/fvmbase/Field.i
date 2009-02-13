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
    ArrayBase* __getitem__(const StorageSite* s)
    {
        return &(self->operator[](*s));
    }
  }
private:
};
