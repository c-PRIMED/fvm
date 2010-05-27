
#ifndef _TIMEDERIVATIVEDISCRETIZATION_KMODEL_H_
#define _TIMEDERIVATIVEDISCRETIZATION_KMODEL_H_

#include <math.h>
#include "CRMatrix.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "GeomFields.h"

template<class X, class Diag, class OffDiag>
class TimeDerivativeDiscretization_Kmodel : public Discretization
{
public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;

  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;


 TimeDerivativeDiscretization_Kmodel(const MeshList& meshes,
				     const GeomFields& geomFields,
				     Field& varField,Field& varN1Field,Field& varN2Field,
				     const T_Scalar dT,
				     const T_Scalar nonDimLength):
				   
  Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _varN1Field(varN1Field),
    _varN2Field(varN2Field),
    _dT(dT),
    _nonDimLength(nonDimLength)
  {}
  
  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {    
    const StorageSite& cells = mesh.getCells();

    //const TArray& dsff = 
    //  dynamic_cast<const TArray&>(_dsff[cells]);

    const TArray& cellVolume = 
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);

    // Field& _varField =  *_dsff.dsf[_direction];  
    //Field& _varN1Field =  *_dsffN1.dsf[_direction]; 
    //Field& _varN2Field =  *_dsffN2.dsf[_direction]; 
  
    
    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    CCMatrix& matrix = 
      dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));

    DiagArray& diag = matrix.getDiag();

    const XArray& x = dynamic_cast<const XArray&>(_varField[cells]);
    const XArray& xN1 = dynamic_cast<const XArray&>(_varN1Field[cells]);
   
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);

    
    const int nCells = cells.getSelfCount();


    if (_varN2Field.hasArray(cells))
      {
        // second order
        const XArray& xN2 = dynamic_cast<const XArray&>(_varN2Field[cells]);
	
        T_Scalar onePointFive(1.5);
        T_Scalar two(2.0);
        T_Scalar pointFive(0.5);
	
	
	for(int c=0; c<nCells; c++)
	  {
	    const T_Scalar fbydT = cellVolume[c]/_dT/pow(_nonDimLength,3);
	    rCell[c] -= fbydT*(onePointFive*x[c]- two*xN1[c]
			       + pointFive*xN2[c]);
	    diag[c] -= fbydT*onePointFive;
	  }
      }
    
    else
      {
	for(int c=0; c<nCells; c++)
	  {
	    const T_Scalar fbydT = cellVolume[c]/_dT/pow(_nonDimLength,3);
	    rCell[c] -= fbydT*(x[c]- xN1[c]);
	    diag[c] -= fbydT;
	  }
      }
    
  }
private:
  const GeomFields& _geomFields;
  const Field& _varField;
  const Field& _varN1Field;
  const Field& _varN2Field;
  //const DistFunctFields<X>& _dsff;  
  //const DistFunctFields<X>& _dsffN1;  
  //const DistFunctFields<X>& _dsffN2; 
  //const int _direction;
  const T_Scalar _dT;
  const T_Scalar _nonDimLength;
};

#endif
