// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _MOMENTUMPRESSUREGRADIENTDISCRETIZATION_H_
#define _MOMENTUMPRESSUREGRADIENTDISCRETIZATION_H_

#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "GeomFields.h"
#include "FlowFields.h"
#include "DiagonalMatrix.h"
#include "Gradient.h"
#include "GradientMatrix.h"
#ifdef PV_COUPLED
#include "CRMatrixRect.h"
#endif

#include "GradientModel.h"

template<class X>
class MomentumPressureGradientDiscretization : public Discretization
{
public:

  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  typedef typename MatrixTraitHelper<X>::T_Matrix CCMatrix;

  typedef Gradient<X> XGrad;
  typedef GradientMatrix<X> GradMatrix;
  typedef typename GradMatrix::Coord VPCoeff;

  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  typedef Array<XGrad> GradArray;

#ifdef PV_COUPLED
  typedef CRMatrixRect<VPCoeff,X,VectorT3> VPMatrix;
#endif

  MomentumPressureGradientDiscretization(const MeshList& meshes,
                                         const GeomFields& geomFields,
                                         FlowFields& flowFields,
                                         GradientModel<X>& pressureGradientModel)  :
    Discretization(meshes),
    _geomFields(geomFields),
    _flowFields(flowFields),
    _pressureGradientModel(pressureGradientModel)
  {}

                          
  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField&, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();

    const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);

    const MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
    VectorT3Array& rCell = dynamic_cast<VectorT3Array&>(rField[vIndex]);
    GradArray& pGradCell = dynamic_cast<GradArray&>(_flowFields.pressureGradient[cells]);
    
    const int nCells = cells.getSelfCount();

    
    const StorageSite& faces = mesh.getFaces();
    const CRConnectivity& faceCells = mesh.getAllFaceCells();
    const VectorT3Array& faceArea =
      dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
    const TArray& facePressure = dynamic_cast<const TArray&>(_flowFields.pressure[faces]);

    pGradCell.zero();

    const int nFaces = faces.getCount();
    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);
        pGradCell[c0].accumulate(faceArea[f],facePressure[f]);
        pGradCell[c1].accumulate(faceArea[f],-facePressure[f]);
    }

    for(int c=0; c<nCells; c++)
      pGradCell[c] /= cellVolume[c];

    // copy boundary values from adjacent cells
    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
        const FaceGroup& fg = *fgPtr;
        const StorageSite& faces = fg.site;
        const CRConnectivity& faceCells = mesh.getFaceCells(faces);
        const int faceCount = faces.getCount();

        if (fg.groupType == "symmetry")
        {
            const VectorT3Array& faceArea =
              dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
            const TArray& faceAreaMag =
              dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);
            for(int f=0; f<faceCount; f++)
            {
                const int c0 = faceCells(f,0);
                const int c1 = faceCells(f,1);
                const VectorT3 en = faceArea[f]/faceAreaMag[f];
                reflectGradient(pGradCell[c1], pGradCell[c0], en);
            }
        }
        else
        {
            for(int f=0; f<faceCount; f++)
            {
                const int c0 = faceCells(f,0);
                const int c1 = faceCells(f,1);
                
                pGradCell[c1] = pGradCell[c0];
            }
        }
    }
  

    for(int c=0; c<nCells; c++)
    {
        rCell[c][0] -= cellVolume[c]*pGradCell[c][0];
        rCell[c][1] -= cellVolume[c]*pGradCell[c][1];
        rCell[c][2] -= cellVolume[c]*pGradCell[c][2];
        
    }
    

#ifdef PV_COUPLED
    const MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
    if (mfmatrix.hasMatrix(vIndex,pIndex))
    {
        VPMatrix& vpMatrix =
          dynamic_cast<VPMatrix&>(mfmatrix.getMatrix(vIndex,pIndex));

        const GradientMatrix<X>& gradMatrix =
          _pressureGradientModel.getGradientMatrix(mesh, _geomFields);

        const Array<VPCoeff>& gradMatrixCoeffs = gradMatrix.getCoeffs();
        Array<VPCoeff>& vpDiag = vpMatrix.getDiag();
        Array<VPCoeff>& vpOffDiag = vpMatrix.getOffDiag();
 
        const CRConnectivity& conn = vpMatrix.getConnectivity();
        const Array<int>& row = conn.getRow();
        
        // all of this assumes that the connectivity of the gradient
        // matrix and vp matrix is the same
        for(int c=0; c<nCells; c++)
        {
            for(int nnb=row[c]; nnb<row[c+1]; nnb++)
            {
                const VPCoeff& coeff = gradMatrixCoeffs[nnb];
                vpDiag[c] += coeff*cellVolume[c];
                vpOffDiag[nnb] -= coeff*cellVolume[c];
            }
        }                
    }
#endif
  }

private:
  const GeomFields& _geomFields;
  FlowFields& _flowFields;
  GradientModel<X>& _pressureGradientModel;
};

#endif
