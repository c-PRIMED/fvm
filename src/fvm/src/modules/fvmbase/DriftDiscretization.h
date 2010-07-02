#ifndef _DRIFTDISCRETIZATION_H_
#define _DRIFTDISCRETIZATION_H_

#include "CRMatrix.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"


template<class X, class Diag, class OffDiag>
class DriftDiscretization : public Discretization
{
public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;  
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;
  typedef Array<int> IntArray;
  typedef Array<X> XArray;
  
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

 
  DriftDiscretization(const MeshList& meshes,
                           const GeomFields& geomFields,
                           Field& varField,
                           const Field& convectingFluxField,
			   const bool useCentralDifference=false) :
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _convectingFluxField(convectingFluxField),
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
   
    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    CCMatrix& matrix =
      dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));

    const CRConnectivity& faceCells = mesh.getAllFaceCells();

    CCAssembler& assembler = matrix.getPairWiseAssembler(faceCells);

    DiagArray& diag = matrix.getDiag();

    const XArray& xCell = dynamic_cast<const XArray&>(xField[cVarIndex]);

    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);

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
	//drift only applies to conductance charge. 
          for(int f=0; f<nFaces; f++)
          {
              const int c0 = faceCells(f,0);
              const int c1 = faceCells(f,1);
              const T_Scalar faceCFlux = convectingFlux[f];

              T_Scalar varFlux;
            
              if (faceCFlux > T_Scalar(0))
              {
                  varFlux = faceCFlux*xCell[c0][1];
                  diag[c0][3] -= faceCFlux;
                  assembler.getCoeff10(f)[3] += faceCFlux;
              }
              else
              {
                  varFlux = faceCFlux*xCell[c1][1];
                  diag[c1][3] += faceCFlux;
                  assembler.getCoeff01(f)[3] -= faceCFlux;
              }

	      rCell[c0][1] -= varFlux;
              rCell[c1][1] += varFlux;
          }
    }


  }
private:
  const GeomFields& _geomFields;
  const Field& _varField;
  const Field& _convectingFluxField; 
  const bool _useCentralDifference;
};

#endif
