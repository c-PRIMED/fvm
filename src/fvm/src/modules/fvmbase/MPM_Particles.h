// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

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
  MPM(string fileName);

  MPM();

  ~MPM();

  typedef Vector<double,3> VecD3;
  typedef Array<VecD3> VecD3Array;

  const StorageSite& getParticles() const {return _particles;}
  
  const StorageSite& getParticles(int num_particles)  {
       _particles.setCount( num_particles);
     return _particles;}
    
  const shared_ptr<Array<VecD3> >& getCoordinates() {return  _coordinates;}
   
  const shared_ptr<Array<VecD3> >& getVelocities() {return  _velocities;}

  const shared_ptr<Array<double> >& getTemperatures() {return  _temperatures;}
  
  const shared_ptr<Array<int> >& getTypes() {return  _types;}
  
  
  void setCoordinates(const shared_ptr<ArrayBase> x) 
  {_coordinates = dynamic_pointer_cast< Array<VecD3> > ( x );}
 
  void setVelocities(const shared_ptr<ArrayBase > v ) 
  {_velocities = dynamic_pointer_cast< Array<VecD3> > (v);}

  void setTypes(const shared_ptr<ArrayBase > type) 
  {_types =  dynamic_pointer_cast< Array<int> > (type);}
 
  void setTemperatures(const shared_ptr<ArrayBase > t) 
  {_temperatures =  dynamic_pointer_cast< Array<double> > (t);}

  void setandwriteParticles(const char *file);

  const shared_ptr<Array<VecD3> > readVelocities(const char *file);

  const shared_ptr<Array<VecD3> > readCoordinates(const char *file);
  
  const shared_ptr<Array<int> > readTypes(const char *file);

  const shared_ptr<Array<double> > readTemperatures(const char *file);

  void Init (const shared_ptr<Array<VecD3> > coordinates,
	     const shared_ptr<Array<VecD3> > velocities,
	     const shared_ptr<Array<int> > types,
	     const shared_ptr<Array<double> > temperatures);
  void Impl(string fileName);
 
 protected:
  StorageSite _particles;
  shared_ptr<Array<VecD3> > _coordinates;
  shared_ptr<Array<VecD3> > _velocities;
  shared_ptr<Array<int> > _types;      //1 for surface particles or 0 for internal particles
  shared_ptr<Array<double> > _temperatures;

};
#endif
