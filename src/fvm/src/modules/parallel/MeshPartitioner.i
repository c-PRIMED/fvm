%{
#include "MeshPartitioner.h"
#include "mpi.h"
%}

%include "std_vector.i"
%include "std_string.i"
%import "Mesh.i"
using namespace std;


class MeshPartitioner{

public:

    enum ETYPE{ TRI = 1, TETRA = 2, HEXA = 3, QUAD = 4 };
    enum WTYPE{ NOWEIGHTS = 0, WEIGHTS_ONLY_EDGES = 1, WEIGTHS_ONLY_VERTICES  = 2,
                WEIGHTS_BOTH_VERTICES_EDGES = 3};
    enum NUMFLAG{ C_STYLE = 0, FORTRAN_STYLE = 1 };
    enum CELLTYPE{ INTERIOR = 1, GHOST_BOUNDARY_CELL = 2, GHOST_INTERFACE_CELL};

    MeshPartitioner(const  MeshList& mesh_list, vector<int> npart,
             vector<int> eType);
    void partition();
    void fiedler_order( const string& fname );
    void mesh();
    const MeshList&  meshList();
    void dumpTecplot();
    void mesh_xdmfplot();
    void isCleanup(bool clean_up);
    void isDebug(bool debug);
    void DEBUG_compute_elem_dist();
    void DEBUG_parmetis_mesh();
    void DEBUG_fiedler_partition();
    void DEBUG_elem_connectivity();
    void DEBUG_map_part_elms();
    void DEBUG_count_elems_part();
    void DEBUG_exchange_part_elems();
    void DEBUG_mapBounIDAndCell();
    void DEBUG_resize_elem();
    void DEBUG_CRConnectivity_cellParts();
    void DEBUG_CRConnectivity_faceParts();
    void DEBUG_interfaces();
    void DEBUG_faceCells_faceNodes();
    void DEBUG_non_interior_cells();
    void DEBUG_preserve_cell_order();
    void DEBUG_order_faceCells_faceNodes();
    void DEBUG_coordinates();
    void DEBUG_exchange_interface_meshes();
    void DEBUG_mesh();
    void DEBUG_local_global();
    void DEBUG_cellcells_global();
    void DEBUG_globalCellID_procID_map();
    void DEBUG_gatherCellsLevel1_partID_map();
    void DEBUG_level1_scatter_gather_cells();
    void DEBUG_CRConnectivity_cellCells2();
    // set property methods
    void setWeightType(MeshPartitioner::WTYPE weight_type);
    void setNumFlag(MeshPartitioner::NUMFLAG num_flag);
 
};


