%module fvmbaseExt
%{
#include "CException.h"
#include <rlog/StdioNode.h>
#include <rlog/RLogChannel.h>
#include "AMG.h"
#include "BCGStab.h"
rlog::StdioNode *stdLog;
%}

%include "std_string.i"
%include "std_vector.i"
%include "std_except.i"
%include "std_map.i"
%include "shared_ptr.i"

namespace std{
%template(vectorStr) vector<string>;
 }


%exception {
    try {
        $action
    }
    catch (CException e) {
        SWIG_exception(SWIG_RuntimeError,e.what());
    }
    catch(...) {
        SWIG_exception(SWIG_RuntimeError,"Unknown exception");
    }
}

using namespace std;

%include "Mesh.i"
%include "AMG.i"
%include "BCGStab.i"

%inline %{
  

  void enableDebug(const string channel)
  {
          stdLog->subscribeTo(rlog::GetGlobalChannel(channel.c_str()));
  }


  %}

%init %{
stdLog = new rlog::StdioNode();

%}
