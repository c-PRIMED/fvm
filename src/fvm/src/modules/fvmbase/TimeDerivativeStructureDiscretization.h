
#ifndef _TIMEDERIVATIVESTRUCTUREDISCRETIZATION_H_
#define _TIMEDERIVATIVESTRUCTUREDISCRETIZATION_H_

#include "CRMatrix.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"

template<class X, class Diag, class OffDiag>
class TimeDerivativeStructureDiscretization : public Discretization
{
public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;

  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;


  TimeDerivativeStructureDiscretization(const MeshList& meshes,
					const GeomFields& geomFields,
					Field& varField,
					Field& varN1Field,
					Field& varN2Field,
					const Field& densityField,
					const Field& volume0Field,
					const T_Scalar dT) :
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _varN1Field(varN1Field),
    _varN2Field(varN2Field),
    _densityField(densityField),
    _volume0Field(volume0Field),
    _dT(dT)
  {}
  
  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {

    const StorageSite& cells = mesh.getCells();

    const TArray& density =
      dynamic_cast<const TArray&>(_densityField[cells]);

    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    
    const TArray& cellVolume0 =
      dynamic_cast<const TArray&>(_volume0Field[cells]);

    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    CCMatrix& matrix =
      dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));

    DiagArray& diag = matrix.getDiag();

    const XArray& x = dynamic_cast<const XArray&>(_varField[cells]);
    const XArray& xN1 = dynamic_cast<const XArray&>(_varN1Field[cells]);
    
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);

    
    const int nCells = cells.getSelfCount();


    // second order
    const XArray& xN2 = dynamic_cast<const XArray&>(_varN2Field[cells]);

    T_Scalar two(2.0);
    const T_Scalar _dT2 = _dT*_dT;
    if (_geomFields.volumeN1.hasArray(cells))
    {
        const TArray& cellVolumeN1 = 
	  dynamic_cast<const TArray&>(_geomFields.volumeN1[cells]);
	const TArray& cellVolumeN2 = 
	  dynamic_cast<const TArray&>(_geomFields.volumeN2[cells]);
	for(int c=0; c<nCells; c++)
	{
	    const T_Scalar rhoVbydT2 = density[c]*cellVolume[c]/_dT2;
	    const T_Scalar rhobydT2 = density[c]/_dT2;
	    const T_Scalar term1 = cellVolume[c];
	    const T_Scalar term2 = two*cellVolumeN1[c];
	    const T_Scalar term3 = cellVolumeN2[c];
	    rCell[c] -= rhobydT2*(term1*x[c]- term2*xN1[c]
				 + term3*xN2[c]);
	    diag[c] -= rhoVbydT2;
	}
    }
    else
    {
        for(int c=0; c<nCells; c++)
	{
	    const T_Scalar rhoVbydT2 = density[c]*cellVolume0[c]/_dT2;
	    rCell[c] -= rhoVbydT2*(x[c]- two*xN1[c]
				  + xN2[c]);
	    diag[c] -= rhoVbydT2;
	}
    }
  }
private:
  const GeomFields& _geomFields;
  const Field& _varField;
  const Field& _varN1Field;
  const Field& _varN2Field;
  const Field& _densityField;
  const Field& _volume0Field;
  const T_Scalar _dT;
};

#endif
