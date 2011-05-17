#ifndef _CONTACTMODEL_H_
#define _CONTACTMODEL_H_

#include "Model.h"
#include "GeomFields.h"
#include "ContactFields.h"
#include "Mesh.h"
#include "Array.h"
#include "Vector.h"


template<class T>
class ContactModel : public Model
{
 public: 

  class Impl;
  
  ContactModel(const GeomFields& geomFields,
	       ContactFields& contactFields,
	       const MeshList& meshes);

  virtual ~ContactModel();

  virtual void init();

  void computeSolidSurfaceForce(const StorageSite& particles);

  void computeSolidSurfaceForcePerUnitArea(const StorageSite& particles); 
  ContactModelConstants<T>& getConstants();

  class NearestCell;
	
 private:
  shared_ptr<Impl> _impl;

};


#endif
