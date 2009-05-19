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
  MPM(string fileName);
  MPM();

  ~MPM();
  
  const StorageSite& getParticles(int num_particles ) {return _particles;}
  
  void setCoordinates(const boost::shared_ptr<ArrayBase> x);
 
  void setVelocities (const boost::shared_ptr<ArrayBase> v);

  void setTypes      (const boost::shared_ptr<ArrayBase> type);

  const StorageSite& getParticles() const {return _particles;}
  
  
  const boost::shared_ptr<ArrayBase> getCoordinates() {return  _coordinates;}
   
  const boost::shared_ptr<ArrayBase> getVelocities() {return  _velocities;}

 
  void setandwriteParticles(const char *file);	
};



