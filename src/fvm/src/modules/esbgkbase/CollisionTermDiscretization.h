#ifndef _COLLISIONTERMDISCRETIZATION_H_
#define _COLLISIONTERMDISCRETIZATION_H_

#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "GeomFields.h"
#include "CRMatrixRect.h"
#include "CRMatrix.h"

template<class X,class Diag, class OffDiag>
class CollisionTermDiscretization : public Discretization
{
public:

  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;

  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;

  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;
 
  CollisionTermDiscretization(const MeshList& meshes,
		       const GeomFields& geomFields,
		       const Field& varField,
		       const Field& sourceField,
		       const Field& collisionFrequency)  :
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
      _sourceField(sourceField),
      _collisionFrequency(collisionFrequency)
   {}

                          
  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField&, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();

    const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);

    const XArray& source = dynamic_cast<const XArray&>(_sourceField[cells]);

    const MultiField::ArrayIndex cVarIndex(&_varField, &cells);
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));
    DiagArray& diag = matrix.getDiag();

    const XArray& x = dynamic_cast<const XArray&>(_varField[cells]);
    TArray& rCell = dynamic_cast<TArray&>(rField[cVarIndex]);   
    
    const XArray& nue = dynamic_cast<const XArray&>(_collisionFrequency[cells]);
    const int nCells = cells.getSelfCount();
    
    for(int c=0; c<nCells; c++)
    {
      rCell[c ] -= cellVolume[c]*nue[c]*(x[c]-source[c]);
      diag[c]-=cellVolume[c]*nue[c];
    }
    cout << "diag[0] = " << diag[0]    
  }
    

private:
  const GeomFields& _geomFields;
  const Field& _varField;
  const Field& _sourceField;
  const Field& _collisionFrequency;
};

#endif

  
    
  

   
        
  
	
