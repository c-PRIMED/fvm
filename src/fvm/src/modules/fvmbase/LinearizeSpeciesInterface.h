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
 

 LinearizeSpeciesInterface( const GeomFields& geomFields,
		      const Field& mass_diffusivity,
		      const T_Scalar A_coeff, 
                      const T_Scalar B_coeff,
		      Field& varField):
    _geomFields(geomFields),
    _varField(varField), 
    _mass_diffusivity(mass_diffusivity),
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
    //const TArray& faceAreaMag =
    //  dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);
    
    const XArray& xCell = dynamic_cast<const XArray&>(xField[cVarIndex]);

    //const VectorT3Array& cellCentroid =
    // dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);

    const CRConnectivity& cellCells = mesh.getCellCells();

    //const TArray& mdCell =
    //  dynamic_cast<const TArray&>(_mass_diffusivity[cells]);
    
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);

    DiagArray& diag = matrix.getDiag();
    
    // In the following, parent is assumed to be on the left, and 
    // the other mesh is assumed to be on the right

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
       OffDiag& dRC0dXC1 = parentmatrix.getCoeff(c0p,  c1p);

       //get other mesh fluxes and coeffs
       const int c0o = otherFaceCells(f,0);
       const int c1o = otherFaceCells(f,1);
       const X rightFlux = rOtherCell[c1o]; // inward shell flux on the right
       OffDiag& dRC0dXC3 = othermatrix.getCoeff(c1o,  c0o);
       OffDiag& dRC0dXC2 = otherdiag[c1o];

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
       diag[c0] = parentdiag[c0p];

       // right shell cell - 1 neighbor
       OffDiag& offdiagC2_C0 = matrix.getCoeff(c2,  c0);

       rCell[c2] = _A_coeff*xCell[c0] + _B_coeff - xCell[c2];
       offdiagC2_C0 = _A_coeff;
       diag[c2] = -1;


       cout << rCell[c0] << "   " << rCell[c2] << endl;
   }

    /*
    const int nCells = cells.getSelfCount();

    for (int c=0; c<(nCells/2); c++)
      {
	const int c0 = c;
	const int nb = cellCells.getCount(c);
	if (nb!=3) 
	  throw CException("invalid connectivity in DoubleShellMesh");
	
	
	  const int c1 = cellCells(c,0);
	  const int c2 = cellCells(c,1);
	  const int c3 = cellCells(c,2);

	  OffDiag& offdiagC1 = matrix.getCoeff(c0,  c1);
	  OffDiag& offdiagC2 = matrix.getCoeff(c0,  c2);
	  OffDiag& offdiagC3 = matrix.getCoeff(c0,  c3);

	  VectorT3 ds1=cellCentroid[c0]-cellCentroid[c1];
	  VectorT3 ds2=cellCentroid[c3]-cellCentroid[c2];

	  T_Scalar ds1Mag = mag(ds1);
	  T_Scalar ds2Mag = mag(ds2);
	  
	  // apply flux balance to left shell cell

	  //left flux (between cells c1 and c0)
	
	  const T_Scalar diffCoeffLeft = mdCell[c1]*faceAreaMag[c0]/ds1Mag;

	  const X dFluxLeft = diffCoeffLeft*(xCell[c1]-xCell[c0]);

	  rCell[c0] += dFluxLeft;

	  offdiagC1 = diffCoeffLeft;

	  diag[c0] = -1*diffCoeffLeft;

	  //right flux (between cells c2 and c3)

	  const T_Scalar diffCoeffRight = mdCell[c3]*faceAreaMag[c0]/ds2Mag;

	  const X dFluxRight = diffCoeffRight*(xCell[c2]-xCell[c3]);

	  rCell[c0] -= dFluxRight;

	  offdiagC2 = 0;//-1*diffCoeffRight;

          offdiagC3 = diffCoeffRight;

	  //apply species jump condition to right shell cell

	  //const T_Scalar alpha = 0.1;

	  OffDiag& offdiagC0 = matrix.getCoeff(c2,  c0);
          
          diag[c2] = -1;

	  offdiagC0 = _A_coeff;

	  rCell[c2] = _A_coeff*xCell[c0] + _B_coeff - xCell[c2];

	  } */
  }
private:
  
  const GeomFields& _geomFields;
  Field& _varField;
  const Field& _mass_diffusivity;
  const T_Scalar _A_coeff;
  const T_Scalar _B_coeff;
  
};

#endif
