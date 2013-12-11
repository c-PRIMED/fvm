// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _LINEARIZEINTERFACEJUMP_H_
#define _LINEARIZEINTERFACEJUMP_H_ 

#include "Mesh.h"
#include "NumType.h"
#include "Array.h"
#include "Vector.h"
#include "Field.h"
#include "CRConnectivity.h"
#include "StorageSite.h"
#include "MultiFieldMatrix.h"
#include "CRMatrix.h"
#include "FluxJacobianMatrix.h"
#include "DiagonalMatrix.h"




template<class X, class Diag, class OffDiag>
class LinearizeInterfaceJump
{
 public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef Array<T_Scalar> TArray;
  typedef Array<X> XArray;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef SquareTensor<T_Scalar,2> SquareTensorT2;
 

 LinearizeInterfaceJump(const T_Scalar A_coeff, 
			   const X B_coeff,
			   Field& varField):
    _varField(varField), 
    _A_coeff(A_coeff),
    _B_coeff(B_coeff)
    {}


    void discretize(const Mesh& mesh, const Mesh& parentMesh, 
		    const Mesh& otherMesh, MultiFieldMatrix& mfmatrix, 
		    MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();
    const StorageSite& faces = mesh.getParentFaceGroupSite();
    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);

    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex)); 

    const XArray& xCell = dynamic_cast<const XArray&>(xField[cVarIndex]);

    const CRConnectivity& cellCells = mesh.getCellCells();
    
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);

    DiagArray& diag = matrix.getDiag();

    // In the following, parent is assumed to be on the left, and 
    // the other mesh is assumed to be on the right when implimenting
    // the (Phi_R = A*Phi_L + B) equation

    // parent mesh info
    const CRConnectivity& parentFaceCells = parentMesh.getFaceCells(faces);
    const StorageSite& parentCells = parentMesh.getCells();
    const MultiField::ArrayIndex cVarIndexParent(&_varField,&parentCells);
    XArray& rParentCell = dynamic_cast<XArray&>(rField[cVarIndexParent]);
    CCMatrix& parentmatrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndexParent,cVarIndexParent)); 
    DiagArray& parentdiag = parentmatrix.getDiag();

    // other mesh info
    const StorageSite& otherFaces = mesh.getOtherFaceGroupSite();
    const CRConnectivity& otherFaceCells = otherMesh.getFaceCells(otherFaces);
    const StorageSite& otherCells = otherMesh.getCells();
    const MultiField::ArrayIndex cVarIndexOther(&_varField,&otherCells);
    XArray& rOtherCell = dynamic_cast<XArray&>(rField[cVarIndexOther]);
    CCMatrix& othermatrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndexOther,cVarIndexOther)); 
    DiagArray& otherdiag = othermatrix.getDiag();

    for (int f=0; f<faces.getCount(); f++)
    {
       //get parent mesh fluxes and coeffs
       int c0p = parentFaceCells(f,0);
       int c1p = parentFaceCells(f,1);
       if (c1p < parentCells.getSelfCount())
       { 
	 // c0 is ghost cell and c1 is boundry cell, so swap cell numbers
	 // so that c1p refers to the ghost cell in the following
	 int temp = c0p;
	 c0p = c1p;
	 c1p = temp;
	 }

       const X leftFlux = rParentCell[c1p]; // inward shell flux on the left
       const OffDiag dRC0dXC3 = parentmatrix.getCoeff(c1p,  c0p);
       const Diag dRC0dXC0 = parentdiag[c1p];

       //get other mesh fluxes and coeffs
       int c0o = otherFaceCells(f,0);
       int c1o = otherFaceCells(f,1);
       if (c1o < otherCells.getSelfCount())
       { 
	 // c0 is ghost cell and c1 is boundry cell, so swap cell numbers
	 // so that c1o refers to the ghost cell in the following
	 int temp = c0o;
	 c0o = c1o;
	 c1o = temp;
	 }

       const X rightFlux = rOtherCell[c1o]; // inward shell flux on the right
       const OffDiag dRC0dXC2 = othermatrix.getCoeff(c1o,  c0o);
       const OffDiag dRC0dXC1 = otherdiag[c1o];


       //now put flux information from meshes into shell cells
       const int c0 = f;
       const int c1 = cellCells(f,0);
       const int c2 = cellCells(f,1);
       const int c3 = cellCells(f,2);

       // left shell cell - 3 neighbors
       OffDiag& offdiagC0_C1 = matrix.getCoeff(c0,  c1);
       OffDiag& offdiagC0_C2 = matrix.getCoeff(c0,  c2);
       OffDiag& offdiagC0_C3 = matrix.getCoeff(c0,  c3);

       rCell[c0] = rightFlux + leftFlux;
       offdiagC0_C1 = dRC0dXC1;
       offdiagC0_C3 = dRC0dXC3;
       offdiagC0_C2 = dRC0dXC2;
       diag[c0] = dRC0dXC0;

       // right shell cell - 1 neighbor
       OffDiag& offdiagC1_C0 = matrix.getCoeff(c1,  c0);

       rCell[c1] = _A_coeff*xCell[c0] + _B_coeff - xCell[c1];
       offdiagC1_C0 = _A_coeff*NumTypeTraits<OffDiag>::getUnity();
       diag[c1] = NumTypeTraits<Diag>::getNegativeUnity();

       // set other coeffs to zero for right shell cell
       OffDiag& offdiagC1_C2 = matrix.getCoeff(c1,  c2);
       OffDiag& offdiagC1_C3 = matrix.getCoeff(c1,  c3);
       offdiagC1_C2 = NumTypeTraits<OffDiag>::getZero();
       offdiagC1_C3 = NumTypeTraits<OffDiag>::getZero();

       if (c0==0)
	 {
	   cout << xCell[c0] << " " << xCell[c1] << " " << xCell[c2] << " " << xCell[c3] << endl;
	 }

   }

  }
private:
  
  Field& _varField;
  const T_Scalar _A_coeff;
  const X _B_coeff;
  
};

#endif
