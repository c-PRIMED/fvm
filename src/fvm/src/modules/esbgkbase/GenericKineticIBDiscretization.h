#ifndef _GENERICKINETICIBDISCRETIZATION_H_
#define _GENERICKINETICIBDISCRETIZATION_H_

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
class GenericKineticIBDiscretization: public Discretization
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

  GenericKineticIBDiscretization (const MeshList& meshes,
				  const GeomFields& geomFields,
				  Field& varField,
				  const double cx,
				  const double cy,
                                  const double cz,
				  MacroFields& macroFields) :
  Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _cx(cx),
    _cy(cy),
    _cz(cz),
    _macroFields(macroFields)
{}

  void discretize(const Mesh& mesh, MultiFieldMatrix& matrix,
                  MultiField& xField, MultiField& rField)
  {
    const StorageSite& ibFaces = mesh.getIBFaces();
     if (ibFaces.getCount() == 0)
      return;
    VectorT3Array& v = dynamic_cast<VectorT3Array&>(_macroFields.velocity[ibFaces]);  
    const StorageSite& cells = mesh.getCells();
    const StorageSite& faces = mesh.getFaces();
    
    const CRConnectivity& faceCells = mesh.getAllFaceCells();

    const MultiField::ArrayIndex xIndex(&_varField,&cells);
    CCMatrix& dRdX = dynamic_cast<CCMatrix&>(matrix.getMatrix(xIndex,xIndex));
    CCAssembler& assembler = dRdX.getPairWiseAssembler(faceCells);
    DiagArray& dRdXDiag=dRdX.getDiag();

    XArray& xcell =
      dynamic_cast<XArray&> (xField[xIndex]);
    XArray& varcell =
      dynamic_cast<XArray&> (_varField[cells]);
    XArray& rCell = dynamic_cast<XArray&>(rField[xIndex]); 
    
    const XArray& xib =
      dynamic_cast<const XArray&>(_varField[mesh.getIBFaces()]);

    const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);

   const IntArray& ibFaceIndex = dynamic_cast<const IntArray&>(_geomFields.ibFaceIndex[faces]);

   const Field& areaMagField = _geomFields.areaMag;
   const XArray& faceAreaMag = dynamic_cast<const XArray&>(areaMagField[faces]);
   const Field& areaField = _geomFields.area;
   const VectorT3Array& faceArea=dynamic_cast<const VectorT3Array&>(areaField[faces]); 

    // used for computing the average value in IBTYPE_BOUNDARY cells.
    // we can't use the xcell storage during the loop below since a
    // cell may have more than one ib face and the current value of
    // x in the cell is required to correctly impose the dirichlet
    // condition. so we accumulate the boundary cell value and counts
    // in the following arrays and after all the faces have been
    // visited, overwrite any boundary cell values with the average of
    // all the ib faces
    const int nCells = cells.getCountLevel1();
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
            const X& xface = xib[ibFace];
	    const X uwall = v[ibFace][0];
	    const X vwall = v[ibFace][1];
	    const X wwall = v[ibFace][2];
            if (ibType[c0] == Mesh::IBTYPE_FLUID)
            {
		const VectorT3 en = faceArea[f]/faceAreaMag[f];
		const X c_dot_en = _cx*en[0]+_cy*en[1]+_cz*en[2];
		const X wallv_dot_en = uwall*en[0]+vwall*en[1]+wwall*en[2];
		if(c_dot_en - wallv_dot_en <0 ){
		  rCell[c0] += assembler.getCoeff01(f)*(xface-varcell[c1]);
		  rCell[c1] = NumTypeTraits<X>::getZero();
		  assembler.getCoeff01(f) = OffDiag(0);
		  dRdX.setDirichlet(c1);
		  xB[c1] += xface;
		  wB[c1]++;
	      }
	    }
	
            else
	      {
		const VectorT3 en = faceArea[f]/faceAreaMag[f];
		const X c_dot_en = _cx*en[0]+_cy*en[1]+_cz*en[2];
		const X wallv_dot_en = uwall*en[0]+vwall*en[1]+wwall*en[2];
		if(c_dot_en - wallv_dot_en> 0){
		  rCell[c1] += assembler.getCoeff10(f)*(xface-varcell[c0]);
		  rCell[c0] = NumTypeTraits<X>::getZero();
		  assembler.getCoeff10(f) = OffDiag(0);
		  dRdX.setDirichlet(c0);
		  xB[c0] += xface;
		  wB[c0]++;
		}	
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
	    assembler.getCoeff01(f)=NumTypeTraits<OffDiag>::getZero();
	    dRdXDiag[c0]=NumTypeTraits<Diag>::getNegativeUnity();
	    assembler.getCoeff10(f)=NumTypeTraits<OffDiag>::getZero();
	    dRdXDiag[c1]=NumTypeTraits<Diag>::getNegativeUnity();

        }
    }

    // set the phi for boundary cells as average of the ib face values
    for(int c=0; c<nCells; c++)
    {
        if (wB[c] > 0)
          varcell[c] = xB[c] / T_Scalar(wB[c]);
	  
    }
  }

private:
  const GeomFields& _geomFields;
  Field& _varField;
  const double _cx;
  const double _cy;
  const double _cz;
  MacroFields& _macroFields;
};

#endif
