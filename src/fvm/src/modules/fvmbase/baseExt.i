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
#include "Array.h"
  
rlog::StdioNode *stdLog;
%}


%include "std_string.i"
%include "std_vector.i"
%include "std_except.i"
%include "std_map.i"


namespace std{
%template(vectorStr) vector<string>;
%template(vectorInt) vector<int>;
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
%include "Vector.i"

typedef Vector<double,3> VecD3;
%template(VecD3) Vector<double,3>;

%include "Field.i"
%include "CRConnectivity.i"
%include "Mesh.i"
%include "LinearSolver.i"
%include "AMG.i"
%include "BCGStab.i"

%include "CellMark.i"
%include "MPM_Particles.i"
%include "Octree.i"
%include "Grid.i"
%include "FVMParticles.i"
%include "MPMCoupling.i"


%inline %{
  void enableDebug(const string channel)
  {
          stdLog->subscribeTo(rlog::GetGlobalChannel(channel.c_str()));
  }

  boost::shared_ptr<ArrayBase> newIntArray(const int size)
  {
    return shared_ptr<ArrayBase>(new Array<int>(size));
  }
  
  %}

%init %{

  import_array();
stdLog = new rlog::StdioNode();

%}
