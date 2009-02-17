%module fvmbaseExt
%{
#define SWIG_FILE_WITH_INIT
  %}
%{
#include "CException.h"
#include <rlog/StdioNode.h>
#include <rlog/RLogChannel.h>
#include "AMG.h"
#include "BCGStab.h"
#include "ArrayBase.h"
  
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
%include "ArrayBase.i"
%include "Mesh.i"
%include "AMG.i"
%include "BCGStab.i"

%include "CellMark.i"
%include "MPM_Particles.i"
%include "Octree.i"

%include "Field.i"
%include "CRConnectivity.i"


%inline %{
  void enableDebug(const string channel)
  {
          stdLog->subscribeTo(rlog::GetGlobalChannel(channel.c_str()));
  }


  %}

%init %{

  import_array();
stdLog = new rlog::StdioNode();

%}
