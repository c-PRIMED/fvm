// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.


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
					Field& varN3Field,
					const Field& densityField,
					const Field& volume0Field,
					const bool& variableTimeStep,
					const T_Scalar dT,
					const T_Scalar dTN1,
					const T_Scalar dTN2) :
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _varN1Field(varN1Field),
    _varN2Field(varN2Field),
    _varN3Field(varN3Field),
    _densityField(densityField),
    _volume0Field(volume0Field),
    _variableTimeStep(variableTimeStep),
    _dT(dT),
    _dTN1(dTN1),
    _dTN2(dTN2)
  {}
  
  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {

    const StorageSite& cells = mesh.getCells();

    const TArray& density =
      dynamic_cast<const TArray&>(_densityField[cells]);

    //const TArray& cellVolume =
    //  dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    
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

    if(!_variableTimeStep)
    {
	if (_varN3Field.hasArray(cells))
	{
	    const XArray& xN3 = dynamic_cast<const XArray&>(_varN3Field[cells]);
	    T_Scalar five(5.0);
	    T_Scalar four(4.0);
	    for(int c=0; c<nCells; c++)
	    {
		const T_Scalar rhoVbydT2 = density[c]*cellVolume0[c]/_dT2;
		rCell[c] -= rhoVbydT2*(two*x[c] - five*xN1[c] + four*xN2[c]
				       - xN3[c]);
		diag[c] -= two*rhoVbydT2;
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
    else
    {
        T_Scalar a = (_dT + _dTN1)/_dT;
        T_Scalar b = (_dT + _dTN1 + _dTN2)/_dT;
        T_Scalar one(1.0);
        if (_varN3Field.hasArray(cells))
        {
            T_Scalar c1 = (two*a*b*(pow(a,two)-pow(b,two))+two*b*(pow(b,two)-one)-two*a*(pow(a,two)-one))/
              (a*b*(a-one)*(b-one)*(a-b));
            T_Scalar c2 = -two*(a+b)/((a-1)*(b-1));
            T_Scalar c3 = -two*(b+one)/(a*(a-b)*(a-one));
            T_Scalar c4 = two*(a+one)/(b*(a-b)*(b-one));
            const XArray& xN3 = dynamic_cast<const XArray&>(_varN3Field[cells]);
            for(int c=0; c<nCells; c++)
            {
                const T_Scalar rhoVbydT2 = density[c]*cellVolume0[c]/_dT2;
                rCell[c] -= rhoVbydT2*(c1*x[c] + c2*xN1[c] + c3*xN2[c]
                                           + c4*xN3[c]);
                diag[c] -= c1*rhoVbydT2;
            }
        }
        else
        {
            T_Scalar c1 = two/a;
            T_Scalar c2 = -two/(a-one);
            T_Scalar c3 = two/(a*(a-one));
            for(int c=0; c<nCells; c++)
            {
                const T_Scalar rhoVbydT2 = density[c]*cellVolume0[c]/_dT2;
                rCell[c] -= rhoVbydT2*(c1*x[c] + c2*xN1[c]
                                           + c3*xN2[c]);
                diag[c] -= c1*rhoVbydT2;
            }
        }
    }

    
  }

  void explicitAdvance(const Mesh& mesh,
                       MultiField& xField,
                       MultiField& rField,
                       MultiField& deltaField)
  {
    
    const StorageSite& cells = mesh.getCells();

    const TArray& density =
      dynamic_cast<const TArray&>(_densityField[cells]);

    //const TArray& cellVolume =
    //  dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    
    const TArray& cellVolume0 =
      dynamic_cast<const TArray&>(_volume0Field[cells]);

    const MultiField::ArrayIndex  cVarIndex(&_varField,&cells);


    const XArray& x = dynamic_cast<const XArray&>(_varField[cells]);
    const XArray& xN1 = dynamic_cast<const XArray&>(_varN1Field[cells]);
    
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);
    XArray& deltaCell = dynamic_cast<XArray&>(deltaField[cVarIndex]);

    
    const int nCells = cells.getSelfCount();


    // second order
    const XArray& xN2 = dynamic_cast<const XArray&>(_varN2Field[cells]);

    T_Scalar two(2.0);
    const T_Scalar _dT2 = _dT*_dT;

    if (_varN3Field.hasArray(cells))
    {
        const XArray& xN3 = dynamic_cast<const XArray&>(_varN3Field[cells]);
	T_Scalar five(5.0);
	T_Scalar four(4.0);
	for(int c=0; c<nCells; c++)
	{
	    const T_Scalar rhoVbydT2 = density[c]*cellVolume0[c]/_dT2;
	    rCell[c] -= rhoVbydT2*(two*x[c] - five*xN1[c] + four*xN2[c]
				   - xN3[c]);
            
	    T_Scalar diag = -two*rhoVbydT2;
            deltaCell[c] = -rCell[c]/diag;
	}
    }
    else
    {
        for(int c=0; c<nCells; c++)
	{
	    const T_Scalar rhoVbydT2 = density[c]*cellVolume0[c]/_dT2;
	    rCell[c] -= rhoVbydT2*(x[c]- two*xN1[c]
				   + xN2[c]);
	    T_Scalar diag = -rhoVbydT2;
            deltaCell[c] = -rCell[c]/diag;
	}
    }
    
  }

private:
  const GeomFields& _geomFields;
  const Field& _varField;
  const Field& _varN1Field;
  const Field& _varN2Field;
  const Field& _varN3Field;
  const Field& _densityField;
  const Field& _volume0Field;
  const bool _variableTimeStep;
  const T_Scalar _dT;
  const T_Scalar _dTN1;
  const T_Scalar _dTN2;
};

#endif
