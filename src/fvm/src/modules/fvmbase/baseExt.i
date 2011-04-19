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
#include "CG.h"
#include "ILU0Solver.h"
#include "JacobiSolver.h"
#include "DirectSolver.h"
#include "SpikeStorage.h"
#include "SpikeSolver.h"
#include "ArrayBase.h"
#include "Array.h"
#include "MatrixOperation.h"
#include <signal.h>
#include <fenv.h>
  
rlog::StdioNode *stdLog;

 void sigFPEHandler(int) throw (CException)
{
  // convert to our exception
  throw CException("floating point exception");
}

%}


%include "std_string.i"
%include "std_vector.i"
%include "std_except.i"
%include "std_map.i"
%include "std_set.i"


namespace std{
%template(vectorStr) vector<string>;
%template(vectorInt) vector<int>;
%template(vectorDouble) vector<double>;
%template(IntSet)    set<int>;
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
%include "CG.i"
%include "JacobiSolver.i"
%include "DirectSolver.i"
%include "SpikeSolver.i"
%include "SpikeStorage.i"




%include "MeshAssembler.i"
%include "MeshDismantler.i"
%include "ILU0Solver.i"
%include "AABB.i"
%include "KSearchTree.i"
%include "IBManager.i"
%include "MatrixOperation.i"

#ifdef FVM_PARALLEL

%include "StorageSiteMerger.i"
%include "LinearSystemMerger.i"
#endif


%inline %{
  void enableDebug(const string channel)
  {
          stdLog->subscribeTo(rlog::GetGlobalChannel(channel.c_str()));
  }
  void enableFPETrapping()
  {
    feenableexcept
      (
       FE_DIVBYZERO
       | FE_INVALID
       | FE_OVERFLOW
       );
    
    signal(SIGFPE, sigFPEHandler);
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
