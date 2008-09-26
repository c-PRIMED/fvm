#ifndef _CONVECTIONDISCRETIZATION_H_
#define _CONVECTIONDISCRETIZATION_H_

#include "CRMatrix.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "UMesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "FieldSet.h"
#include "ConnectivityField.h"
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
  
  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  typedef Array<XGrad> GradArray;

  
  ConvectionDiscretization(const Args& args):
    Discretization(args),
    _varField(getChildRef<Field>("varField")),
    _convectingFluxField(getChildRef<Field>("convectingFluxField")),
    _geomFields(getChildRef<FieldSet>("geomFields")),
    _coordField(_geomFields.getChildRef<Field>("coordinates")),
    _areaField(_geomFields.getChildRef<Field>("area")),
    _areaMagField(_geomFields.getChildRef<Field>("areaMagnitude")),
    _varGradientField(getChildRef<Field>("varGradientField")),
    _continuityResidualField(getChildRef<Field>("continuityResidualField"))
  {}

  void discretize(const Mesh& gmesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {

    const UMesh& mesh = SafeCast<UMesh>(gmesh);
    const StorageSite& cells = mesh.getCells();
    const StorageSite& faces = mesh.getFaces();


    // should there be some other checks ?
    if (!_convectingFluxField.hasArray(faces))
      return;

    const TArray& convectingFlux = SafeCast<TArray>(_convectingFluxField[faces]);
    const TArray& continuityResidual = SafeCast<TArray>(_continuityResidualField[cells]);

    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    CCMatrix& matrix = SafeCast<CCMatrix>(mfmatrix.getMatrix(cVarIndex,cVarIndex));

    const CRConnectivity& faceCells = mesh.getAllFaceCells();

    CCAssembler& assembler = matrix.getPairWiseAssembler(faceCells);
    DiagArray& diag = matrix.getDiag();

    const XArray& xCell = SafeCast<XArray>(xField[cVarIndex]);
    XArray& rCell = SafeCast<XArray>(rField[cVarIndex]);

    //const GradArray& xGradCell = SafeCast<GradArray>(_varGradientField[cells]);

    
    const int nFaces = faces.getCount();
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
    }

    const int nCells = cells.getSelfCount();
    for(int c=0;c<nCells;c++)
    {
        const T_Scalar cImb = continuityResidual[c];
        diag[c] += cImb;
    }

  }
  DECLARE_HT("ConvectionDiscretization<"
             +NumTypeTraits<X>::getTypeName()+","
             +NumTypeTraits<Diag>::getTypeName()+","
             +NumTypeTraits<OffDiag>::getTypeName()
             +">");

private:
  const Field& _varField;
  const Field& _convectingFluxField; 
  const FieldSet& _geomFields;
  const Field& _coordField;
  const Field& _areaField;
  const Field& _areaMagField;
  const Field& _varGradientField;
  const Field& _continuityResidualField;
};

template<class X, class Diag, class OffDiag>
void
ConvectionDiscretization<X,Diag,OffDiag>::addMethods()
{
  INHERIT_METHODS(Discretization);
}

REGISTER_HT_TEMPLATE(MULTI_ARG3(<class X, class Diag, class OffDiag>), ConvectionDiscretization,
                     MULTI_ARG3(<X,Diag,OffDiag>));


#endif
