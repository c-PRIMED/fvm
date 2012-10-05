// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _ELECDIFFUSIONDISCRETIZATION_H_
#define _ELECDIFFUSIONDISCRETIZATION_H_

#include "CRMatrix.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "DiagonalMatrix.h"
#include "Gradient.h"
#include "DiagonalTensor.h"
/*
template<class T>
inline T harmonicAverage(const T& x0, const T& x1)
{
  const T sum = x0+x1;
  if (x0+x1 != NumTypeTraits<T>::getZero())
    return 2.0*x0*x1/sum;
  else
    return sum;
}
*/
  
template<class X, class Diag, class OffDiag>
class ElecDiffusionDiscretization : public Discretization
{
public:

  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;

  typedef Gradient<X> XGrad;
  
  typedef Array<int> IntArray;
  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  typedef Array<XGrad> GradArray;
  
  ElecDiffusionDiscretization(const MeshList& meshes,
                          const GeomFields& geomFields,
                          Field& varField,
                          const Field& diffusivityField,
                          const Field& varGradientField) :
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _diffusivityField(diffusivityField),
    _varGradientField(varGradientField)
  {}

  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();
    const StorageSite& faces = mesh.getFaces();
    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,
                                                             cVarIndex));

    const VectorT3Array& faceArea =
      dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
    
    const TArray& faceAreaMag =
      dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);

    const VectorT3Array& cellCentroid =
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);

    const VectorT3Array& faceCentroid =
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[faces]);

    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    
    const CRConnectivity& faceCells = mesh.getAllFaceCells();

    CCAssembler& assembler = matrix.getPairWiseAssembler(faceCells);
    DiagArray& diag = matrix.getDiag();

    const XArray& xCell = dynamic_cast<const XArray&>(xField[cVarIndex]);
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);

    const GradArray& xGradCell =
      dynamic_cast<const GradArray&>(_varGradientField[cells]);

    const TArray& diffCell =
      dynamic_cast<const TArray&>(_diffusivityField[cells]);

    const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);

    const int nFaces = faces.getCount();
    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);

        T_Scalar vol0 = cellVolume[c0];
        T_Scalar vol1 = cellVolume[c1];
        
        VectorT3 ds=cellCentroid[c1]-cellCentroid[c0];

        // for ib faces ignore the solid cell and use the face centroid for diff metric
        if (((ibType[c0] == Mesh::IBTYPE_FLUID)
             && (ibType[c1] == Mesh::IBTYPE_BOUNDARY)) ||
            ((ibType[c1] == Mesh::IBTYPE_FLUID)
             && (ibType[c0] == Mesh::IBTYPE_BOUNDARY)))
        {
            if (ibType[c0] == Mesh::IBTYPE_FLUID)
            {
                vol1 = 0.;
                ds = faceCentroid[f]-cellCentroid[c0];
            }
            else
            {
                vol0 = 0.;
                ds = cellCentroid[c1]-faceCentroid[f];
            }
        }
        
        T_Scalar faceDiffusivity(1.0);
        if (vol0 == 0.)
          faceDiffusivity = diffCell[c1];
        else if (vol1 == 0.)
          faceDiffusivity = diffCell[c0];
        else
          faceDiffusivity = harmonicAverage(diffCell[c0],diffCell[c1]);
	
	const T_Scalar diffMetric = faceAreaMag[f]*faceAreaMag[f]/dot(faceArea[f],ds);
        const T_Scalar diffCoeff = faceDiffusivity*diffMetric;
        const VectorT3 secondaryCoeff = faceDiffusivity*(faceArea[f]-ds*diffMetric);
        
        const XGrad gradF = (xGradCell[c0]*vol0+xGradCell[c1]*vol1)/(vol0+vol1);

        const X dFluxSecondary = gradF*secondaryCoeff;
        
        const X dFlux = diffCoeff*(xCell[c1]-xCell[c0]) + dFluxSecondary;

        rCell[c0][1] += dFlux[1];
        rCell[c1][1] -= dFlux[1];
        
        assembler.getCoeff01(f)[3] +=diffCoeff;
        assembler.getCoeff10(f)[3] +=diffCoeff;

        diag[c0][3] -= diffCoeff;
        diag[c1][3] -= diffCoeff;
    }
  }
private:
  const GeomFields& _geomFields;
  Field& _varField;
  const Field& _diffusivityField; 
  const Field& _varGradientField;
};

#endif
