%{
#include "MeshAssembler.h"
%}

class MeshAssembler
{
public:
    
   MeshAssembler( const MeshList& meshList );
   ~MeshAssembler();
 
  const MeshList&  meshList();
   void  debug_print();

};
