%{
#include "MeshDismantler.h"
%}

class MeshDismantler
{
public:
    
   MeshDismantler( const MeshList& meshList );
   ~MeshDismantler();
 
  const MeshList&  meshList();
   void  debug_print();

};
