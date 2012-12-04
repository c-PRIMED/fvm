// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _CONVECTIONDISCRETIZATION_KMODEL_H_
#define _CONVECTIONDISCRETIZATION_KMODEL_H_

#include "CRMatrix.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "Gradient.h"

template<class X, class Diag, class OffDiag>
class ConvectionDiscretization_Kmodel : public Discretization
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

  ConvectionDiscretization_Kmodel(const MeshList& meshes,
				  const GeomFields& geomFields,
				  Field& varField,
				  const double cx,
				  const double cy,
				  const double cz,
				  //const double nondim_length,
				  //const double Lx,
				  //const double Ly,
				  //const double Lz,
				  bool useCentralDifference) :
    //const bool useCentralDifference=false) :

  Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField), 
    _cx(cx),
    _cy(cy),
    _cz(cz),
    _useCentralDifference(useCentralDifference)
  {}
  
  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {

    const StorageSite& cells = mesh.getCells();
    const StorageSite& faces = mesh.getFaces();


    // should there be some other checks ?
    //if (!_convectingFluxField.hasArray(faces))
    //  return;

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

    const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);    
    //const X nondim_length=_options["nonDimLt"];
    //const X Lx=_options["nonDimLx"];
    //const X Ly=_options["nonDimLy"];
    //const X Lz=_options["nonDimLz"];
    
    if (_geomFields.gridFlux.hasArray(faces))
    {
        shared_ptr<TArray> gridFluxPtr(new TArray(nFaces));
	TArray& gridFlux = *gridFluxPtr;
        gridFlux = dynamic_cast<const TArray&>(_geomFields.gridFlux[faces]);

	for(int f=0; f<nFaces; f++)
	{
            const int c0 = faceCells(f,0);
	    const int c1 = faceCells(f,1);
	    //const T_Scalar faceCFlux = convectingFlux[f] - gridFlux[f];
	    const T_Scalar faceCFlux = faceArea[f][0]*_cx+faceArea[f][1]*_cy+faceArea[f][2]*_cz - gridFlux[f]; 
	    //const T_Scalar faceCFlux = faceArea[f][0]*_cx*nondim_length/Lx+faceArea[f][1]*_cy*nondim_length/Ly+faceArea[f][2]*_cz*nondim_length/Lz - gridFlux[f];

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
	      const T_Scalar faceCFlux = faceArea[f][0]*_cx+faceArea[f][1]*_cy+faceArea[f][2]*_cz;
              //const T_Scalar faceCFlux = faceArea[f][0]*_cx*nondim_length/Lx+faceArea[f][1]*_cy*nondim_length/Ly+faceArea[f][2]*_cz*nondim_length/Lz;
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
              const T_Scalar faceCFlux = faceArea[f][0]*_cx+faceArea[f][1]*_cy+faceArea[f][2]*_cz;
	      //const T_Scalar faceCFlux = faceArea[f][0]*_cx*nondim_length/Lx+faceArea[f][1]*_cy*nondim_length/Ly+faceArea[f][2]*_cz*nondim_length/Lz;


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

    //const int nCells = cells.getSelfCount();
    //for(int c=0;c<nCells;c++)
    //{
    //    const T_Scalar cImb = continuityResidual[c];
    //    diag[c] += cImb;
    //}

  }
private:
  const GeomFields& _geomFields;
  const Field& _varField;
  const double _cx;
  const double _cy;
  const double _cz;
  const bool _useCentralDifference;  
  KineticModelOptions<X> _options;
};

#endif
