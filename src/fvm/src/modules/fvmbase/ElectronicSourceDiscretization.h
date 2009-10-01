#ifndef _ELECTRONICSOURCEDISCRETIZATION_H_
#define _ELECTRONICSOURCEDISCRETIZATION_H_

#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "GeomFields.h"
#include "CRMatrixRect.h"

template<class X>
class ElectronicSourceDiscretization : public Discretization
{
public:

  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  typedef Array<X> XArray;

  typedef Array<T_Scalar> TArray;
 
  ElectronicSourceDiscretization(const MeshList& meshes,
                                         const GeomFields& geomFields,
                                         ElectronicFields& electronicFields)  :
    Discretization(meshes),
    _geomFields(geomFields),
    _electronicFields(electronicFields)
   {}

                          
  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField&, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();

    const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);

    const TArray& source = dynamic_cast<const TArray&>(_electronicFields.source[cells]);

    const MultiField::ArrayIndex sIndex(&_electronicFields.potential,&cells);

    TArray& rCell = dynamic_cast<TArray&>(rField[sIndex]);   
    
    const int nCells = cells.getSelfCount();
    
    for(int c=0; c<nCells; c++)
    {
        rCell[c ] += cellVolume[c]*source[c];
    }
  }
    

private:
  const GeomFields& _geomFields;
  ElectronicFields& _electronicFields;
};

#endif
