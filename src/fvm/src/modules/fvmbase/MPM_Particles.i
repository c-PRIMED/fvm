%{

#include "Array.h"
#include "StorageSite.h"
#include "Vector.h"
#include "Field.h"
#include "CRConnectivity.h"
  %}

using namespace std;



class MPM
{
 public: 
  MPM();

  ~MPM();

  const StorageSite& getParticles() const {return _particles;}
  
  const boost::shared_ptr<ArrayBase> getCoordinates() {return  _coordinates;}
   
  const boost::shared_ptr<ArrayBase> getVelocities() {return  _velocities;}

  void Impl(string fileName);

  void setandwriteParticles(const char *file);	
};



