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
  typedef Array<VecD3> VecD3Array;

  const StorageSite& getParticles() const {return _particles;}
  
  const shared_ptr<Array<VecD3> >& getCoordinates() {return  _coordinates;}
   
  const shared_ptr<Array<VecD3> >& getVelocities() {return  _velocities;}

  const shared_ptr<Array<int> >& getTypes() {return  _types;}
  
  
  void setCoordinates(const shared_ptr<Array<VecD3> > x) {_coordinates = x;}
 
  void setVelocities(const shared_ptr<Array<VecD3> > x) {_velocities = x;}

  void setTypes(const shared_ptr<Array<int> > x) {_types = x;}
 
  void setandwriteParticles(const char *file);

  const shared_ptr<Array<VecD3> > readVelocities(const char *file);

  const shared_ptr<Array<VecD3> > readCoordinates(const char *file);
  
  const shared_ptr<Array<int> > readTypes(const char *file);
  
  void Init (const shared_ptr<Array<VecD3> > coordinates,
	     const shared_ptr<Array<VecD3> > velocities,
	     const shared_ptr<Array<int> > types);
  void Impl(string fileName);
 
 protected:
  StorageSite _particles;
  shared_ptr<Array<VecD3> > _coordinates;
  shared_ptr<Array<VecD3> > _velocities;
  shared_ptr<Array<int> > _types;      //1 for surface particles or 0 for internal particles


};
#endif
