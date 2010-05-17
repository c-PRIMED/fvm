%{
#include "MeshPartitioner.h"
#include "mpi.h"
%}

%include "std_vector.i"
%include "std_string.i"
%import "Mesh.i"
using namespace std;
%template(IntVector) vector<int>;


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
    // set property methods
    void setWeightType(MeshPartitioner::WTYPE weight_type);
    void setNumFlag(MeshPartitioner::NUMFLAG num_flag);
 
};


