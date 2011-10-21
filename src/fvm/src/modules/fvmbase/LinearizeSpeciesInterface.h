#ifndef _LINEARIZESPECIESINTERFACE_H_
#define _LINEARIZESPECIESINTERFACE_H_ 

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
class LinearizeSpeciesInterface
{
 public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef Array<T_Scalar> TArray;
  typedef Array<X> XArray;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
 

 LinearizeSpeciesInterface(const T_Scalar A_coeff, 
			   const T_Scalar B_coeff,
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
    const CRConnectivity& otherFaceCells = otherMesh.getFaceCells(faces);
    const StorageSite& otherCells = otherMesh.getCells();
    const MultiField::ArrayIndex cVarIndexOther(&_varField,&otherCells);
    XArray& rOtherCell = dynamic_cast<XArray&>(rField[cVarIndexOther]);
    CCMatrix& othermatrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndexOther,cVarIndexOther)); 
    DiagArray& otherdiag = othermatrix.getDiag();

    for (int f=0; f<faces.getCount(); f++)
    {
       //get parent mesh fluxes and coeffs
       const int c0p = parentFaceCells(f,0);
       const int c1p = parentFaceCells(f,1);
       const X leftFlux = rParentCell[c0p]; // inward shell flux on the left
       OffDiag dRC0dXC1 = parentmatrix.getCoeff(c1p,  c0p);
       Diag dRC0dXC0 = parentdiag[c1p];
       if (c1p < parentCells.getSelfCount())
       {
	 // c0 is ghost cell and c1 is boundry cell, so fix coeffs
	 dRC0dXC1 = parentmatrix.getCoeff(c0p,  c1p);
	 dRC0dXC0 = parentdiag[c0p];
       }

       //get other mesh fluxes and coeffs
       const int c0o = otherFaceCells(f,0);
       const int c1o = otherFaceCells(f,1);
       const X rightFlux = rOtherCell[c1o]; // inward shell flux on the right
       OffDiag dRC0dXC3 = othermatrix.getCoeff(c1o,  c0o);
       OffDiag dRC0dXC2 = otherdiag[c1o];
       if (c1o < otherCells.getSelfCount())
       {
	 // c0 is ghost cell and c1 is boundry cell, so fix coeffs
	 dRC0dXC3 = othermatrix.getCoeff(c0o,  c1o);
	 dRC0dXC2 = otherdiag[c0o];
       }

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
       OffDiag& offdiagC2_C0 = matrix.getCoeff(c2,  c0);

       rCell[c2] = _A_coeff*xCell[c0] + _B_coeff - xCell[c2];
       offdiagC2_C0 = _A_coeff;
       diag[c2] = -1;
   }

  }
private:
  
  Field& _varField;
  const T_Scalar _A_coeff;
  const T_Scalar _B_coeff;
  
};

#endif
