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
   void  debug_sites();
   void  debug_localToGlobal_mappers();
   void  debug_globalCellToMeshID_mappers();
   void  debug_sync_localToGlobal_mappers();
   void  debug_faceCells();
   void  debug_localNodeToGlobal();

};
