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
				const Field& eta1Field,
				const Field& varGradientField)  :
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _etaField(etaField),
    _eta1Field(eta1Field),
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

    const TArray& eta1Cell =
      dynamic_cast<const TArray&>(_eta1Field[cells]);

    const int nFaces = faces.getCount();

    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);

	const VectorT3& Af = faceArea[f];

        T vol0 = cellVolume[c0];
        T vol1 = cellVolume[c1];

        T faceEta(1.0);
	T faceEta1(1.0);

        if (vol0 == 0.)
       	{
            faceEta = etaCell[c1];
	    faceEta1 = eta1Cell[c1];
	}
        else if (vol1 == 0.)
	{
            faceEta = etaCell[c0];
	    faceEta1 = eta1Cell[c0];
	}
        else
	{
            faceEta = harmonicAverage(etaCell[c0],etaCell[c1]);
	    faceEta1 = harmonicAverage(eta1Cell[c0],eta1Cell[c1]);
	}
	
	const VGradType gradF = (vGradCell[c0]*vol0+vGradCell[c1]*vol1)/(vol0+vol1);

	VectorT3 source(NumTypeTraits<VectorT3>::getZero());

	source[0] = faceEta*
	                     (gradF[0][0]*Af[0]
			     +gradF[0][1]*Af[1]
			     +gradF[0][2]*Af[2])
	           +faceEta1*
	                     (gradF[0][0]
			     +gradF[1][1]
			     +gradF[2][2])*Af[0];
	source[1] = faceEta*
	                     (gradF[1][0]*Af[0]
			     +gradF[1][1]*Af[1]
			     +gradF[1][2]*Af[2])
	           +faceEta1*
	                     (gradF[0][0]
	                     +gradF[1][1]
	                     +gradF[2][2])*Af[1];
	source[2] = faceEta*
	                     (gradF[2][0]*Af[0]
			     +gradF[2][1]*Af[1]
			     +gradF[2][2]*Af[2])
	           +faceEta1*
	                     (gradF[0][0]
	                     +gradF[1][1]
	                     +gradF[2][2])*Af[2];
	rCell[c0] += source;
	rCell[c1] -= source;
    }
  }
    

private:
  const GeomFields& _geomFields;
  Field& _varField;
  const Field& _etaField;
  const Field& _eta1Field;
  const Field& _varGradientField;
};

#endif
