#ifndef _CAPTUREDISCRETIZATION_H_
#define _CAPTUREDISCRETIZATION_H_

#include "PhysicsConstant.h"
#include "ElectricBC.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "DiagonalMatrix.h"
#include "CRMatrix.h"

//Diag type: 2x2Tensor
/***  | d0,  d1 |
      | d2,  d3 |    ***/
//OffDiag type: 2x2Tensor
/***  | o0,  o1 |
      | o2,  o3 |    ***/
//X type: VectorT2
/***  | x0 |
      | x1 |         ***/


// in capture model, d1, d3 and x0, x1 are modified
// the rest remains unchanged

template <class X, class Diag, class OffDiag>
class CaptureDiscretization : public Discretization
{

 public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  typedef Array<T_Scalar> TArray;
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::OffDiagArray OffDiagArray;
  typedef Array<X> XArray;

  CaptureDiscretization(const MeshList& meshes,
			const GeomFields& geomFields,
			const Field& varField,
			const Field& varN1Field,
			const Field& totaltrapsField,
			const Field& capturecrossField,
			const ElectricModelConstants<T_Scalar>& constants):
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _varN1Field(varN1Field),
    _totaltrapsField(totaltrapsField),
    _capturecrossField(capturecrossField),
    _constants(constants)
      {}

  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();
    
    const int nCells = cells.getSelfCount();

    const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    const TArray& electron_totaltraps = dynamic_cast<const TArray&> (_totaltrapsField[cells]);

    const TArray& free_electron_capture_cross = dynamic_cast<const TArray&> (_capturecrossField[cells]);

    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));
    
    const XArray& xCell = dynamic_cast<const XArray&>(xField[cVarIndex]);

    const XArray& xN1Cell = dynamic_cast<const XArray&>(_varN1Field[cells]);
    
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);
       
    DiagArray& diag = matrix.getDiag();

    OffDiagArray& offdiag = matrix.getOffDiag();

    const T_Scalar& electron_effmass = _constants["electron_effmass"];
    const T_Scalar& temperature = _constants["OP_temperature"];
    
    const T_Scalar velocity = sqrt(8.0 * K_SI * temperature / (PI * ME * electron_effmass));
    
    for(int c=0; c<nCells; c++){
      
      T_Scalar fluxCoeff = cellVolume[c] * velocity * free_electron_capture_cross[c];
     
      rCell[c][0] += fluxCoeff * (xN1Cell[c][1] * electron_totaltraps[c] - xN1Cell[c][1] * xN1Cell[c][0]); 
      //diag[c][0] -= fluxCoeff * xCell[c][1];  
      
      rCell[c][1] -= fluxCoeff * (xN1Cell[c][1] * electron_totaltraps[c] - xN1Cell[c][1] * xN1Cell[c][0]); 
      //diag[c][1] -= fluxCoeff * (electron_totaltraps[c]-xCell[c][0]); 

    }    
  }

 private:
  const GeomFields& _geomFields;
  const Field& _varField;
  const Field& _varN1Field;
  const Field& _totaltrapsField;
  const Field& _capturecrossField;
  const ElectricModelConstants<T_Scalar>& _constants;
};


#endif
