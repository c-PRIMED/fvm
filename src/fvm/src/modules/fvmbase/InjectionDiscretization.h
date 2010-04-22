#ifndef _INJECTIONDISCRETIZATION_H_
#define _INJECTIONDISCRETIZATION_H_

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

// in injection model, only x1 is modified
// the rest remain unchanged

template <class X, class Diag, class OffDiag>
class InjectionDiscretization : public Discretization
{

 public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  typedef Array<T_Scalar> TArray;
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::OffDiagArray OffDiagArray;
  typedef Array<X> XArray;

  InjectionDiscretization(const MeshList& meshes,
			  const GeomFields& geomFields,
			  const Field& varField,
			  const Field& electricField,
			  const Field& conductionbandField,
			  const ElectricModelConstants<T_Scalar>& constants,
			  Field& transmissionField) :
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _electricField(electricField),
    _conductionbandField(conductionbandField),
    _constants(constants),
    _transmissionField(transmissionField)  
      {}

  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();
    
    const int nCells = cells.getCount();
    
    const TArray& electric_field = dynamic_cast<const TArray&> (_electricField[cells]);

    const TArray& conduction_band = dynamic_cast<const TArray&> (_conductionbandField[cells]);

    TArray& transmission = dynamic_cast<TArray&> (_transmissionField[cells]);

    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));
    
    const XArray& xCell = dynamic_cast<const XArray&>(xField[cVarIndex]);
    
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);
       
    DiagArray& diag = matrix.getDiag();

    OffDiagArray& offdiag = matrix.getOffDiag();

    T_Scalar currentTime = 1.0;
    const T_Scalar volt = getMembraneVoltage(currentTime);
    const T_Scalar& fermilevelsubstrate = _constants["substrate_function"];
    const T_Scalar& fermilevelmembrane = _constants["substrate_function"] - volt;
    const T_Scalar& electron_effmass = _constants["electron_effmass"];
    const T_Scalar& temperature = _constants["OP_temperature"];


    const T_Scalar& alpha = 4.0 * PI * electron_effmass * ME * QE / pow(H_SI, 3.0);

    T_Scalar fermilevel, supplyfunction, flux;

    
    for(int c=0; c<nCells; c++){
      
      const T_Scalar en = conduction_band[c];
      
      const T_Scalar energystep = electric_field[c]; //need to change for 3D
      
      ElectronTransmissionCoefficient(en, transmission);

      //injection substrate
      
      fermilevel = fermilevelsubstrate;
      
      supplyfunction = ElectronSupplyFunction(en, fermilevel, temperature);

      flux = alpha * transmission[c] * energystep * supplyfunction;
      
      rCell[c][1] += flux;   //sign?

      //injection membrane

      fermilevel = fermilevelmembrane;
      
      supplyfunction = ElectronSupplyFunction(en, fermilevel, temperature);
       
      flux = alpha * transmission[c] * energystep * supplyfunction;
      
      rCell[c][1] += flux;   //sign?
    }
  }


 private:
  const GeomFields& _geomFields;
  const Field& _varField;
  const Field& _electricField;
  const Field& _conductionbandField;
  const ElectricModelConstants<T_Scalar>& _constants;
  Field& _transmissionField;

};





#endif
