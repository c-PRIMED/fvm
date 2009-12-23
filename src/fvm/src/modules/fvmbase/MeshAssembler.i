%{
#include "MeshAssembler.h"
%}

class MeshAssembler
{
public:
    
   MeshAssembler( const MeshList& meshList );
   ~MeshAssembler();

   void  debug_print();

};
