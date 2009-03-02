#ifndef _MOMENTUMIBDISCRETIZATION_H_
#define _MOMENTUMIBDISCRETIZATION_H_

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
class MomentumIBDiscretization: public Discretization
{
public:

  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::OffDiagArray OffDiagArray;

  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  MomentumIBDiscretization (const MeshList& meshes,
		     const GeomFields& geomFields,
		     FlowFields& flowFields) :
    Discretization(meshes),
    _geomFields(geomFields),
    _flowFields(flowFields)
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

    const MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(vIndex,vIndex));
    CCAssembler& assembler = matrix.getPairWiseAssembler(faceCells);

    VectorT3Array& cellVelocity =
      dynamic_cast<VectorT3Array&> (_flowFields.velocity[cells]);

    XArray& rCell = dynamic_cast<XArray&>(rField[vIndex]); 
    
    const VectorT3Array& ibVelocity =
      dynamic_cast<const VectorT3Array&>(_flowFields.velocity[mesh.getIBFaces()]);

    const Array<int>& ibType = mesh.getIBType();
    const int nIBFaces = ibFaces.getCount();

    // used to keep track of the current ib face index
    int ibFace =0;
    
    const int nFaces = faces.getCount();
    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);

        if (((ibType[c0] == Mesh::IBTYPE_FLUID) && (ibType[c1] == Mesh::IBTYPE_BOUNDARY)) ||
            ((ibType[c1] == Mesh::IBTYPE_FLUID) && (ibType[c0] == Mesh::IBTYPE_BOUNDARY)))
        {
            if (ibFace >= nIBFaces)
              throw CException("incorrect number of IB faces");
            // this is an iBFace, determine which cell is interior and which boundary
            const VectorT3& faceVelocity = ibVelocity[ibFace];
            if (ibType[c0] == Mesh::IBTYPE_FLUID)
            {
                rCell[c0] += assembler.getCoeff01(f)*(faceVelocity-cellVelocity[c1]);
                rCell[c1].zero();
                assembler.getCoeff01(f) = OffDiag(0);
                matrix.setDirichlet(c1);
            }
            else
            {
                rCell[c1] += assembler.getCoeff10(f)*(faceVelocity-cellVelocity[c0]);
                rCell[c0].zero();
                assembler.getCoeff10(f) = OffDiag(0);
                matrix.setDirichlet(c0);
            }
            ibFace++;
        }
        else if ((ibType[c0] == Mesh::IBTYPE_FLUID) &&
            (ibType[c1] == Mesh::IBTYPE_FLUID))
        {
            // leave as  is
        }
        else
        {
            // setup to get zero corrections
            rCell[c0].zero();
            rCell[c1].zero();
            matrix.setDirichlet(c0);
            matrix.setDirichlet(c1);
        }
    }
    if (ibFace != nIBFaces)
      throw CException("incorrect number of IB faces");
  }

private:
  const GeomFields& _geomFields;
  FlowFields& _flowFields;
};

#endif
