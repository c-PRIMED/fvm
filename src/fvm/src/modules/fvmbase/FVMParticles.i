%{
#include "FVMParticles.h"
%}

%include "std_vector.i"
%import "Mesh.i"
%import "ArrayBase.i"
%template(IntVector) vector<int>;

using namespace std;


class FVMParticles{

public:

     FVMParticles( const MeshList& meshList );
    ~FVMParticles();

     void   setParticles(int nsweep);
     void   setSweepIter( int sweep);

      const ArrayBase&  getCellIDs( int mesh_id ) const;
      int  getNumOfFluidParticles ( int mesh_id ) const;
};
