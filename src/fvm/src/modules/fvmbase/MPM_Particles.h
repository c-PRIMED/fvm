#ifndef _MPM_H_
#define _MPM_H_

#include "Array.h"
#include "StorageSite.h"
#include "Vector.h"
#include "Field.h"
#include "CRConnectivity.h"


class MPM
{
 public: 
  MPM();

  ~MPM();

  typedef Vector<double,3> VecD3;

  const StorageSite& getParticles() const {return _particles;}
  
  shared_ptr<Array<VecD3> > getCoordinates() {return  _coordinates;}
   
  shared_ptr<Array<VecD3> > getVelocities() {return  _velocities;}
  
  
  void setCoordinates(shared_ptr<Array<VecD3> > x) {_coordinates = x;}
 
  void setVelocities(shared_ptr<Array<VecD3> > x) {_velocities = x;}
 
  void setandwriteParticles(char *file);

  const shared_ptr<Array<VecD3> > readVelocities(char *file);

  const shared_ptr<Array<VecD3> > readCoordinates(char *file);
   
  void Init (shared_ptr<Array<VecD3> > coordinates,
	     shared_ptr<Array<VecD3> > velocities );

 
 protected:
  StorageSite _particles;
  shared_ptr<Array<VecD3> > _coordinates;
  shared_ptr<Array<VecD3> > _velocities;


};
#endif
