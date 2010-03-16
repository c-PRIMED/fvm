%{
#include "MeshDismantler.h"
%}

class MeshDismantler
{
public:
    
   MeshDismantler( const MeshList& meshList );
   ~MeshDismantler();
 
   const MeshList&  meshList();

   //debug files
   void debug_print();
   void debug_cell_site();
   void debug_face_site();
   void debug_node_site();
   void debug_cells_mapper();
   void debug_face_cells();
   void debug_nodes_mapper();
   void debug_face_nodes();
   void debug_scatter_mappers();
   void debug_gather_mappers();

};
