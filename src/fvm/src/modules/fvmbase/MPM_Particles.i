%{

#include "Array.h"
#include "StorageSite.h"
#include "Vector.h"
#include "Field.h"
#include "CRConnectivity.h"
  %}

using namespace std;

%include "ArrayBase.i"
%include "Vector.i" 
class MPM
{
 public: 
  MPM();

  ~MPM();

 const StorageSite& getParticles() const {return _particles;}
  
 const shared_ptr<Array<VecD3> >& getCoordinates() {return  _coordinates;}
   
 const shared_ptr<Array<VecD3> >& getVelocities() {return  _velocities;}


};



