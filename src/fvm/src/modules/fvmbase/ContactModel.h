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


  struct NearestCell
  {
     NearestCell():
        mesh(0),
        cell(-1),
	distanceSquared(0)
     {}

     const Mesh* mesh;
     int cell;
     double distanceSquared;
     set<int> neighbors;
  };

 private:
  shared_ptr<Impl> _impl;

};


#endif
