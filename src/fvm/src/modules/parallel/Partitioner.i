%module fvmparallel
%{
#define SWIG_FILE_WITH_INIT
%}
%{
  #include "PartMesh.h"
  #include "MeshPartitioner.h"
%}

%include "PartMesh.i"
%include "MeshPartitioner.i"
