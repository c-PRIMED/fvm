// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _LINEARIZEINTERFACEJUMPUNCONNECTED_H_
#define _LINEARIZEINTERFACEJUMPUNCONNECTED_H_ 

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
class LinearizeInterfaceJumpUnconnected
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
 

 LinearizeInterfaceJumpUnconnected(const T_Scalar A_coeff, 
			   const X B_coeff,
			   Field& varField):
    _varField(varField), 
    _A_coeff(A_coeff),
    _B_coeff(B_coeff)
    {}


    void discretize(const Mesh& mesh, const Mesh& parentMesh, 
		    const Mesh& otherMesh, MultiFieldMatrix& mfmatrix, 
		    MultiField& xField, MultiField& rField, LinearSystem& lsShell)
  {
    cout << "LINERIZEINTERFACEJUMPUNCONNECTED" << endl;

    const StorageSite& cells = mesh.getCells();
    const StorageSite& faces = mesh.getParentFaceGroupSite();
    const CRConnectivity& cellCells = mesh.getCellCells();
    XArray& varCell = dynamic_cast<XArray&>(_varField[cells]);

    //lsShell info
    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    CCMatrix& matrix = dynamic_cast<CCMatrix&>((lsShell.getMatrix()).getMatrix(cVarIndex,cVarIndex)); 
    const XArray& xCell = dynamic_cast<const XArray&>((lsShell.getX())[cVarIndex]);
    XArray& rCell = dynamic_cast<XArray&>((lsShell.getB())[cVarIndex]);
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
    const XArray& varCellParent = dynamic_cast<const XArray&>(_varField[parentCells]);


    // other mesh info
    const StorageSite& otherFaces = mesh.getOtherFaceGroupSite();
    const CRConnectivity& otherFaceCells = otherMesh.getFaceCells(otherFaces);
    const StorageSite& otherCells = otherMesh.getCells();
    const MultiField::ArrayIndex cVarIndexOther(&_varField,&otherCells);
    XArray& rOtherCell = dynamic_cast<XArray&>(rField[cVarIndexOther]);
    CCMatrix& othermatrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndexOther,cVarIndexOther)); 
    DiagArray& otherdiag = othermatrix.getDiag();
    const XArray& varCellOther = dynamic_cast<const XArray&>(_varField[otherCells]);


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

       // copy sln varialbe values from interior cells of meshes to ghost cells of shell mesh
       varCell[c3] = varCellParent[c0p];
       varCell[c2] = varCellOther[c0o];
       
       // calculate residuals and jacobian values for R_1 equation
       const X r0 = -1*(rightFlux + leftFlux);
       const X r1 = -1*(_A_coeff*varCell[c0] + _B_coeff - varCell[c1]);
       const OffDiag dRC1dXC2 = NumTypeTraits<OffDiag>::getZero();
       const Diag dRC1dXC1 = NumTypeTraits<Diag>::getNegativeUnity();
       const OffDiag dRC1dXC3 = NumTypeTraits<OffDiag>::getZero();
       const OffDiag dRC1dXC0 = _A_coeff*NumTypeTraits<OffDiag>::getUnity();

       // set variables to match hand written derivation for ease of coding
       const OffDiag A = dRC1dXC0; // A_coeff // 1
       const Diag B = dRC1dXC1; // -1
       const OffDiag C = dRC1dXC2;  // 0
       const OffDiag D = dRC1dXC3; // 0
       const X E = _B_coeff; // 0
       const Diag F = dRC0dXC0;
       const OffDiag G = dRC0dXC1;
       const OffDiag H = dRC0dXC2;
       const OffDiag J = dRC0dXC3;
       const X K = NumTypeTraits<X>::getZero(); // 0
       
       // manipulate parent matrix R and coeffs to incorporate above
       OffDiag& offDiagParentC0_C1 = parentmatrix.getCoeff(c0p,  c1p);
       Diag& diagParentC0 = parentdiag[c0p];
       X& parentR0 = rParentCell[c0p];
       
       if (f==0)
	 {
	   cout << "orignialParentOffDiag " << offDiagParentC0_C1 << endl;
	   cout <<  "originalDiagparent" << diagParentC0 << endl;
	 }
       const OffDiag originalParentOffDiag = offDiagParentC0_C1;
       offDiagParentC0_C1 = originalParentOffDiag*(B*H-G*C)/(A*G-B*F);
       diagParentC0 += originalParentOffDiag*(B*J-D*G)/(A*G-B*F);
       //parentR0 += originalParentOffDiag*(G*r1 - B*r0 + B*K - E*G)/(A*G-B*F);
       parentR0 += originalParentOffDiag*(G*r1 - B*r0)/(A*G-B*F);
       if (f==0)
	 {
	   cout << A << " " << B << " " << C << " " << D << " " << E << " " << F << " " << G << " " << H << " " << J << " " << K << " " <<  originalParentOffDiag << " " << r0 << " " << r1 << endl;
	   cout << (B*H-G*C) << " " << (B*J-D*G) << " " << (A*G-B*F) << endl;
	   cout << offDiagParentC0_C1 << " " << diagParentC0 << endl;
	 }

       // manipulate other matrix R and coeffs to incorporate above
       OffDiag& offDiagOtherC0_C1 = othermatrix.getCoeff(c0o,  c1o);
       Diag& diagOtherC0 = otherdiag[c0o];
       X& otherR0 = rOtherCell[c0o];

       const OffDiag originalOtherOffDiag = offDiagOtherC0_C1;
       offDiagOtherC0_C1 = originalOtherOffDiag*(D*F-A*J)/(A*G-B*F);
       diagOtherC0 += originalOtherOffDiag*(C*F-A*H)/(A*G-B*F);
       //otherR0 += originalOtherOffDiag*(A*r0 - F*r1 + E*F - A*K)/(A*G-B*F);
       otherR0 += originalOtherOffDiag*(A*r0 - F*r1)/(A*G-B*F);
       
       //editing ghost cells (not sure if I need to, seems not to matter as expected)
       // but makes looking for differences easier, so keep for now.
       // should also affect interface flux calculation
       //parent
       
       OffDiag& offDiagParentC1_C0 = parentmatrix.getCoeff(c1p,  c0p);
       Diag& diagParentC1 = parentdiag[c1p];
       X& parentR1 = rParentCell[c1p];

       const Diag originalParentDiagC1 = diagParentC1;//offDiagParentC0_C1;
       offDiagParentC1_C0 += diagParentC1*(B*J-D*G)/(A*G-B*F);
       diagParentC1 *= (B*H-G*C)/(A*G-B*F);
       parentR1 += originalParentDiagC1*(G*r1 - B*r0)/(A*G-B*F);

       //other
       OffDiag& offDiagOtherC1_C0 = othermatrix.getCoeff(c1o,  c0o);
       Diag& diagOtherC1 = otherdiag[c1o];
       X& otherR1 = rOtherCell[c1o];

       const Diag originalOtherDiagC1 = diagOtherC1; //offDiagOtherC0_C1;
       offDiagOtherC1_C0 += diagOtherC1*(C*F-A*H)/(A*G-B*F);
       diagOtherC1 *= (D*F-A*J)/(A*G-B*F);
       otherR1 += originalOtherDiagC1*(A*r0 - F*r1)/(A*G-B*F);




       //now put flux information from meshes into shell cells
       // left shell cell
       OffDiag& offdiagC0_C1 = matrix.getCoeff(c0,  c1);
       OffDiag& offdiagC0_C2 = matrix.getCoeff(c0,  c2);
       OffDiag& offdiagC0_C3 = matrix.getCoeff(c0,  c3);

       rCell[c0] = -r0;
       offdiagC0_C1 = dRC0dXC1;
       offdiagC0_C3 = dRC0dXC3;
       offdiagC0_C2 = dRC0dXC2;
       diag[c0] = dRC0dXC0;

       // right shell cell
       OffDiag& offdiagC1_C0 = matrix.getCoeff(c1,  c0);
       OffDiag& offdiagC1_C2 = matrix.getCoeff(c1,  c2);
       OffDiag& offdiagC1_C3 = matrix.getCoeff(c1,  c3);

       rCell[c1] = -r1;
       offdiagC1_C0 = dRC1dXC0;
       offdiagC1_C2 = dRC1dXC2;
       offdiagC1_C3 = dRC1dXC3;
       diag[c1] = dRC1dXC1;
       
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
