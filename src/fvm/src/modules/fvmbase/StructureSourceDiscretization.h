#ifndef _STRUCTURESOURCEDISCRETIZATION_H_
#define _STRUCTURESOURCEDISCRETIZATION_H_

#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "GeomFields.h"
#include "CRConnectivity.h"
#include "CRMatrixRect.h"
#include "Vector.h"

template<class T, class Diag, class OffDiag>
class StructureSourceDiscretization : public Discretization
{
public:

  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef Gradient<VectorT3> VGradType;
  typedef Array<Gradient<VectorT3> > VGradArray;

  typedef CRMatrix<Diag,OffDiag,VectorT3> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;

  StructureSourceDiscretization(const MeshList& meshes,
				const GeomFields& geomFields,
				Field& varField,
				const Field& muField,
				const Field& lambdaField,
				const Field& varGradientField)  :
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _muField(muField),
    _lambdaField(lambdaField),
    _varGradientField(varGradientField)
   {}

                          
  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();
    const StorageSite& faces = mesh.getFaces();
    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);

    const VectorT3Array& faceArea =
      dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);

    const TArray& faceAreaMag =
      dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);

    const VectorT3Array& cellCentroid =
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);

    const VectorT3Array& faceCentroid =
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[faces]);


    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);

    const CRConnectivity& faceCells = mesh.getAllFaceCells();

    VectorT3Array& rCell = dynamic_cast<VectorT3Array&>(rField[cVarIndex]);
    const VectorT3Array& xCell = dynamic_cast<const VectorT3Array&>(xField[cVarIndex]);

    const VGradArray& vGradCell =
      dynamic_cast<const VGradArray&>(_varGradientField[cells]);

    const TArray& muCell =
      dynamic_cast<const TArray&>(_muField[cells]);

    const TArray& lambdaCell =
      dynamic_cast<const TArray&>(_lambdaField[cells]);

    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,
                                                             cVarIndex));
    CCAssembler& assembler = matrix.getPairWiseAssembler(faceCells);
    DiagArray& diag = matrix.getDiag();


    const int nFaces = faces.getCount();

    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);

	const VectorT3& Af = faceArea[f];
        VectorT3 ds=cellCentroid[c1]-cellCentroid[c0];

        T vol0 = cellVolume[c0];
        T vol1 = cellVolume[c1];

        T faceMu(1.0);
	T faceLambda(1.0);

        if (vol0 == 0.)
       	{
            faceMu = muCell[c1];
	    faceLambda = lambdaCell[c1];
	}
        else if (vol1 == 0.)
	{
            faceMu = muCell[c0];
	    faceLambda = lambdaCell[c0];
	}
        else
	{
            faceMu = harmonicAverage(muCell[c0],muCell[c1]);
	    faceLambda = harmonicAverage(lambdaCell[c0],lambdaCell[c1]);
	}
	
	const VGradType gradF = (vGradCell[c0]*vol0+vGradCell[c1]*vol1)/(vol0+vol1);

	VectorT3 source(NumTypeTraits<VectorT3>::getZero());
        const T divU = (gradF[0][0] + gradF[1][1] + gradF[2][2]);
        
        // mu*grad U ^ T + lambda * div U I 
	source[0] = faceMu*(gradF[0][0]*Af[0] + gradF[0][1]*Af[1] + gradF[0][2]*Af[2])
          + faceLambda*divU*Af[0];
        
	source[1] = faceMu*(gradF[1][0]*Af[0] + gradF[1][1]*Af[1] + gradF[1][2]*Af[2])
          + faceLambda*divU*Af[1];
        
	source[2] = faceMu*(gradF[2][0]*Af[0] + gradF[2][1]*Af[1] + gradF[2][2]*Af[2])
          + faceLambda*divU*Af[2];


        // mu*gradU, primary part
        
        const T diffMetric = faceAreaMag[f]*faceAreaMag[f]/dot(faceArea[f],ds);

        source += faceMu*diffMetric*(xCell[c1]-xCell[c0]);

        // mu*gradU, secondart part

        const VectorT3 secondaryCoeff = faceMu*(faceArea[f]-ds*diffMetric);

        source += gradF*secondaryCoeff;

        // add flux to the residual of c0 and c1
        rCell[c0] += source;
	rCell[c1] -= source;

        // for Jacobian, use 2*mu + lambda as the diffusivity
        const T faceDiffusivity = T(2.0) * faceMu + faceLambda;
        const T diffCoeff = faceDiffusivity*diffMetric;

        OffDiag& a01 = assembler.getCoeff01(f);
        OffDiag& a10 = assembler.getCoeff10(f);
        
        a01 +=diffCoeff;
        a10 +=diffCoeff;

        Diag& a00 = diag[c0];
        Diag& a11 = diag[c1];
        
        a00 -= diffCoeff;
        a11 -= diffCoeff;
    }
  }
    

private:
  const GeomFields& _geomFields;
  Field& _varField;
  const Field& _muField;
  const Field& _lambdaField;
  const Field& _varGradientField;
};

#endif
