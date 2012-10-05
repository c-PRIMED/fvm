// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _CONVECTIONDISCRETIZATION_H_
#define _CONVECTIONDISCRETIZATION_H_

#include "CRMatrix.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "Gradient.h"

template<class X, class Diag, class OffDiag>
class ConvectionDiscretization : public Discretization
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

  ConvectionDiscretization(const MeshList& meshes,
                           const GeomFields& geomFields,
                           Field& varField,
                           const Field& convectingFluxField,
                           const Field& continuityResidualField,
                           const Field& varGradientField,
                           const bool useCentralDifference=false) :
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _convectingFluxField(convectingFluxField),
    _continuityResidualField(continuityResidualField),
    _varGradientField(varGradientField),
    _useCentralDifference(useCentralDifference)
  {}
  
  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {

    const StorageSite& cells = mesh.getCells();
    const StorageSite& faces = mesh.getFaces();


    // should there be some other checks ?
    if (!_convectingFluxField.hasArray(faces))
      return;

    const TArray& convectingFlux =
      dynamic_cast<const TArray&>(_convectingFluxField[faces]);
    const TArray& continuityResidual =
      dynamic_cast<const TArray&>(_continuityResidualField[cells]);

    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    CCMatrix& matrix =
      dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));

    const CRConnectivity& faceCells = mesh.getAllFaceCells();

    CCAssembler& assembler = matrix.getPairWiseAssembler(faceCells);
    DiagArray& diag = matrix.getDiag();

    const XArray& xCell = dynamic_cast<const XArray&>(xField[cVarIndex]);
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);

    //const GradArray& xGradCell = dynamic_cast<GradArray>(_varGradientField[cells]);

    const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
    
    const int nFaces = faces.getCount();
    if (_geomFields.gridFlux.hasArray(faces))
    {
        shared_ptr<TArray> gridFluxPtr(new TArray(nFaces));
	TArray& gridFlux = *gridFluxPtr;
        gridFlux = dynamic_cast<const TArray&>(_geomFields.gridFlux[faces]);

	for(int f=0; f<nFaces; f++)
	{
            const int c0 = faceCells(f,0);
	    const int c1 = faceCells(f,1);
	    const T_Scalar faceCFlux = convectingFlux[f] - gridFlux[f];

	    X varFlux;
	    if (faceCFlux > T_Scalar(0))
	    {
	        varFlux = faceCFlux*xCell[c0];
                diag[c0] -= faceCFlux;
                assembler.getCoeff10(f) += faceCFlux;
	    }
	    else
	    {
                varFlux = faceCFlux*xCell[c1];
                diag[c1] += faceCFlux;
                assembler.getCoeff01(f)-= faceCFlux;
	    }
        
	    rCell[c0] -= varFlux;
	    rCell[c1] += varFlux;
	}
    }
    else
    {
      if (_useCentralDifference){
	for(int f=0; f<nFaces; f++)
          {
              const int c0 = faceCells(f,0);
              const int c1 = faceCells(f,1);
              const T_Scalar faceCFlux = convectingFlux[f];
              bool isIBFace = (((ibType[c0] == Mesh::IBTYPE_FLUID)
                                && (ibType[c1] == Mesh::IBTYPE_BOUNDARY)) ||
                               ((ibType[c1] == Mesh::IBTYPE_FLUID)
                                && (ibType[c0] == Mesh::IBTYPE_BOUNDARY)));
            

              X varFlux =0.5*faceCFlux*(xCell[c0] + xCell[c0]);

              rCell[c0] -= varFlux;
              rCell[c1] += varFlux;

              if (isIBFace)
              {
                  // linearize the actual flux as calculated
                  // above. this will ensure that the Ib
                  // discretization will be able to fix the value
                  // correctly using the ib face value
                  
                  diag[c0] -= 0.5*faceCFlux;
                  assembler.getCoeff10(f) -= 0.5*faceCFlux;
                  diag[c1] += 0.5*faceCFlux;
                  assembler.getCoeff01(f) += 0.5*faceCFlux;
              }
              else
              {
                  // linearize as upwind flux so that linear system
                  // remains diagonally dominant
                  if (faceCFlux > T_Scalar(0))
                  {
                      diag[c0] -= faceCFlux;
                      assembler.getCoeff10(f) += faceCFlux;
                  }
                  else
                  {
                      diag[c1] += faceCFlux;
                      assembler.getCoeff01(f)-= faceCFlux;
                  }
              }
          }
       }
       else
          for(int f=0; f<nFaces; f++)
          {
              const int c0 = faceCells(f,0);
              const int c1 = faceCells(f,1);
              const T_Scalar faceCFlux = convectingFlux[f];

              X varFlux;
            
              if (faceCFlux > T_Scalar(0))
              {
                  varFlux = faceCFlux*xCell[c0];
                  diag[c0] -= faceCFlux;
                  assembler.getCoeff10(f) += faceCFlux;
              }
              else
              {
                  varFlux = faceCFlux*xCell[c1];
                  diag[c1] += faceCFlux;
                  assembler.getCoeff01(f)-= faceCFlux;
              }

              rCell[c0] -= varFlux;
              rCell[c1] += varFlux;
              //cout << "convflux" << varFlux << endl;

          }
    }

    const int nCells = cells.getSelfCount();
    for(int c=0;c<nCells;c++)
    {
        const T_Scalar cImb = continuityResidual[c];
        diag[c] += cImb;
    }

  }
private:
  const GeomFields& _geomFields;
  const Field& _varField;
  const Field& _convectingFluxField; 
  const Field& _continuityResidualField;
  const Field& _varGradientField;
  const bool _useCentralDifference;
};

#endif
