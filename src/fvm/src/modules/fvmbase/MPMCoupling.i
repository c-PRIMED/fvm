%{
#include "MPMCoupling.h"
%}

%include "std_vector.i"
%import "Mesh.i"
%import "ArrayBase.i"
%template(IntVector) vector<int>;

using namespace std;


class MPMCoupling{

public:

    MPMCoupling( MPM& mpm, Field& coordinate, Field& velocity, StorageSite& cell_site );
   ~MPMCoupling();
       
    void updateMPM();
    void acceptMPM();

};
