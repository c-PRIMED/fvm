// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _SOURCEDISCRETIZATIONFORFRACTURE_H_
#define _SOURCEDISCRETIZATIONFORFRACTURE_H_

#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "GeomFields.h"
#include "CRMatrixRect.h"
#include "Vector.h"

template<class T, class Diag, class OffDiag>
class SourceDiscretizationforFracture : public Discretization
{
public:

  //typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  //typedef Array<X> XArray;

  typedef Array<T> TArray;
  
  typedef Vector<T,3> VectorT3;
  
  typedef CRMatrix<Diag,OffDiag,T> CCMatrix;
  
  typedef typename CCMatrix::DiagArray DiagArray;
 
  SourceDiscretizationforFracture(const MeshList& meshes,
		       const GeomFields& geomFields,
		       const Field& varField,
		       const Field& sourceField,
		       const Field& sourcecoefField)  :
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _sourceField(sourceField),
    _sourcecoefField(sourcecoefField)
   {}

                          
  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField&, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();

    const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);

    const TArray& source = dynamic_cast<const TArray&>(_sourceField[cells]);
    
    const TArray& sourcecoef = dynamic_cast<const TArray&>(_sourcecoefField[cells]);

    const MultiField::ArrayIndex cVarIndex(&_varField, &cells);

    TArray& rCell = dynamic_cast<TArray&>(rField[cVarIndex]);   
    
    const int nCells = cells.getSelfCount();
    
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));
    
    DiagArray& diag = matrix.getDiag();
    
    for(int c=0; c<nCells; c++)
    {
        rCell[c] += cellVolume[c]*(1.0+source[c]);
		diag[c] += cellVolume[c]*sourcecoef[c];

	  //if (c >= 400){
	  //cout << c << ": " << cellVolume[c] << " " << diag[c] << "\n" << endl;}
    }
  }
    

private:
  const GeomFields& _geomFields;
  const Field& _varField;
  const Field& _sourceField;
  const Field& _sourcecoefField;
};

#endif
