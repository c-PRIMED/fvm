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

template<class T>
class StructureSourceDiscretization : public Discretization
{
public:

  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef Gradient<VectorT3> VGradType;
  typedef Array<Gradient<VectorT3> > VGradArray;
 
  StructureSourceDiscretization(const MeshList& meshes,
				const GeomFields& geomFields,
				Field& varField,
				const Field& etaField,
				const Field& lambdaField,
				const Field& varGradientField)  :
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _etaField(etaField),
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

    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);

    const CRConnectivity& faceCells = mesh.getAllFaceCells();

    VectorT3Array& rCell = dynamic_cast<VectorT3Array&>(rField[cVarIndex]);

    const VGradArray& vGradCell =
      dynamic_cast<const VGradArray&>(_varGradientField[cells]);

    const TArray& etaCell =
      dynamic_cast<const TArray&>(_etaField[cells]);

    const TArray& lambdaCell =
      dynamic_cast<const TArray&>(_lambdaField[cells]);

    const int nFaces = faces.getCount();

    VectorT3 source(NumTypeTraits<VectorT3>::getZero());

    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);

        T vol0 = cellVolume[c0];
        T vol1 = cellVolume[c1];

        T faceEta(1.0);
	T faceLambda(1.0);

        if (vol0 == 0.)
       	{
            faceEta = etaCell[c1];
	    faceLambda = lambdaCell[c1];
	}
        else if (vol1 == 0.)
	{
            faceEta = etaCell[c0];
	    faceLambda = lambdaCell[c0];
	}
        else
	{
            faceEta = harmonicAverage(etaCell[c0],etaCell[c1]);
	    faceLambda = harmonicAverage(lambdaCell[c0],lambdaCell[c1]);
	}

        const VGradType gradF = (vGradCell[c0]*vol0+vGradCell[c1]*vol1)/(vol0+vol1);

	source[0] = faceEta*
	                     ((gradF[0])[0]*(faceArea[f])[0]
			     +(gradF[0])[1]*(faceArea[f])[1]
			     +(gradF[0])[2]*(faceArea[f])[2])
	           +faceLambda*
	                     ((gradF[0])[0]
			     +(gradF[1])[1]
			     +(gradF[2])[2])*(faceArea[f])[0];
	source[1] = faceEta*
	                     ((gradF[1])[0]*(faceArea[f])[0]
			     +(gradF[1])[1]*(faceArea[f])[1]
			     +(gradF[1])[2]*(faceArea[f])[2])
	           +faceLambda*
	                     ((gradF[0])[0]
	                     +(gradF[1])[1]
	                     +(gradF[2])[2])*(faceArea[f])[1];
	source[2] = faceEta*
	                     ((gradF[2])[0]*(faceArea[f])[0]
			     +(gradF[2])[1]*(faceArea[f])[1]
			     +(gradF[2])[2]*(faceArea[f])[2])
                   +faceLambda*
	                     ((gradF[0])[0]
	                     +(gradF[1])[1]
	                     +(gradF[2])[2])*(faceArea[f])[2];
	rCell[c0] += source;
	rCell[c1] -= source;
    }
  }
    

private:
  const GeomFields& _geomFields;
  Field& _varField;
  const Field& _etaField;
  const Field& _lambdaField;
  const Field& _varGradientField;
};

#endif
