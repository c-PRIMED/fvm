// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _GENERICIBDISCRETIZATION_H_
#define _GENERICIBDISCRETIZATION_H_

#include "CRMatrix.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "DiagonalMatrix.h"
#include "DiagonalTensor.h"
#include "FlowFields.h"
#include "GeomFields.h"
#include "Vector.h"
#include "NumType.h"

  
template <class X, class Diag, class OffDiag>
class GenericIBDiscretization: public Discretization
{
public:

  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::OffDiagArray OffDiagArray;

  typedef Array<int> IntArray;
  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  GenericIBDiscretization (const MeshList& meshes,
		     const GeomFields& geomFields,
		     Field& phiField) :
    Discretization(meshes),
    _geomFields(geomFields),
    _phiField(phiField)
{}

  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {
    const StorageSite& ibFaces = mesh.getIBFaces();

    if (ibFaces.getCount() == 0)
      return;
    
    const StorageSite& cells = mesh.getCells();
    const StorageSite& faces = mesh.getFaces();
    
    const CRConnectivity& faceCells = mesh.getAllFaceCells();

    const MultiField::ArrayIndex vIndex(&_phiField,&cells);
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(vIndex,vIndex));
    CCAssembler& assembler = matrix.getPairWiseAssembler(faceCells);

    XArray& cellPhi =
      dynamic_cast<XArray&> (_phiField[cells]);

    XArray& rCell = dynamic_cast<XArray&>(rField[vIndex]); 
    
    const XArray& ibPhi =
      dynamic_cast<const XArray&>(_phiField[mesh.getIBFaces()]);

    const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);

   const IntArray& ibFaceIndex = dynamic_cast<const IntArray&>(_geomFields.ibFaceIndex[faces]);


    // used for computing the average value in IBTYPE_BOUNDARY cells.
    // we can't use the cellPhi storage during the loop below since a
    // cell may have more than one ib face and the current value of
    // phi in the cell is required to correctly impose the dirichlet
    // condition. so we accumulate the boundary cell value and counts
    // in the following arrays and after all the faces have been
    // visited, overwrite any boundary cell values with the average of
    // all the ib faces
    const int nCells = cells.getCount();
    XArray xB(nCells);
    Array<int> wB(nCells);

    xB.zero();
    wB.zero();
      
    
    const int nFaces = faces.getCount();
    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);

        if (((ibType[c0] == Mesh::IBTYPE_FLUID) && (ibType[c1] == Mesh::IBTYPE_BOUNDARY)) ||
            ((ibType[c1] == Mesh::IBTYPE_FLUID) && (ibType[c0] == Mesh::IBTYPE_BOUNDARY)))
        {
            // this is an iBFace, determine which cell is interior and which boundary

            const int ibFace = ibFaceIndex[f];
            if (ibFace < 0)
              throw CException("invalid ib face index");
            
            const X& facePhi = ibPhi[ibFace];
            if (ibType[c0] == Mesh::IBTYPE_FLUID)
            {
                rCell[c0] += assembler.getCoeff01(f)*(facePhi-cellPhi[c1]);
                rCell[c1] = NumTypeTraits<X>::getZero();
                assembler.getCoeff01(f) = OffDiag(0);
                matrix.setDirichlet(c1);
                xB[c1] += facePhi;
                wB[c1]++;
            }
            else
            {
                rCell[c1] += assembler.getCoeff10(f)*(facePhi-cellPhi[c0]);
		rCell[c0] = NumTypeTraits<X>::getZero();
                assembler.getCoeff10(f) = OffDiag(0);
                matrix.setDirichlet(c0);
                xB[c0] += facePhi;
                wB[c0]++;
            }
        }
        else if ((ibType[c0] == Mesh::IBTYPE_FLUID) &&
            (ibType[c1] == Mesh::IBTYPE_FLUID))
        {
            // leave as  is
        }
        else
        {
            // setup to get zero corrections
            rCell[c0] = NumTypeTraits<X>::getZero();
            rCell[c1] = NumTypeTraits<X>::getZero();
            matrix.setDirichlet(c0);
            matrix.setDirichlet(c1);
        }
    }

    // set the phi for boundary cells as average of the ib face values
    for(int c=0; c<nCells; c++)
    {
        if (wB[c] > 0)
          cellPhi[c] = xB[c] / T_Scalar(wB[c]);
    }
  }

private:
  const GeomFields& _geomFields;
  Field& _phiField;
};

#endif
