#ifndef _SOURCEDISCRETIZATIONENE_H_
#define _SOURCEDISCRETIZATIONENE_H_

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
#include "GradientModel.h"


template<class T>
class SourceDiscretizationene : public Discretization
{
public:
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef Gradient<VectorT3> VGradType;
  typedef Array<Gradient<VectorT3> > VGradArray;
  typedef Array<T> TArray;


  SourceDiscretizationene(const MeshList& meshes,
                          const GeomFields& geomFields,
                          Field& varField,
                          const Field& velocityField,
                          Field& muField,
                          Field& dissipationField,
                          Field& densityField,
                          Field& sourcekField,
                          const Field& gradientField) :

      Discretization(meshes),
      _geomFields(geomFields),
      _varField(varField),
      _velocityField(velocityField),
      _muField(muField),
      _dissipationField(dissipationField),
      _densityField(densityField),
      _sourcekField(sourcekField),
      _gradientField(gradientField)

   {}
  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {


    const StorageSite& cells = mesh.getCells();
    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);

    const VGradArray& vGrad =
      dynamic_cast<const VGradArray&>(_gradientField[cells]);

    const TArray& muCell =
      dynamic_cast<TArray&>(_muField[cells]);

    const TArray& eCell =
      dynamic_cast<TArray&>(_dissipationField[cells]);

    const TArray & rhoCell =
      dynamic_cast<const TArray&>(_densityField[cells]);

    TArray & sourceCell =
      dynamic_cast<TArray&>(_sourcekField[cells]);

  
    TArray& rCell =
      dynamic_cast<TArray&>(rField[cVarIndex]);
   
  
    const int nCells = cells.getCount();  

    T sum;
    for(int n=0; n<nCells; n++)
    {
        const VGradType& vg = vGrad[n];
        VGradType vgSquare = vGrad[n]; 
        T x = eCell[n]*rhoCell[n];       

        for(int i=0;i<3;i++)
        {
         for(int j=0;j<3;j++)
      
         {
           vgSquare[i][j] =  vg[i][j]*vg[i][j]+vg[i][j]*vg[j][i] ;
           sum += vgSquare[i][j];
           //vgSquare[i][j] += vg[i][j]*vg[j][i];
          // sourceCell[n] = vgSquare[i][j]*muCell[n]-x;
         }
        }
        sourceCell[n] = sum*muCell[n]-x;
        rCell[n] +=sourceCell[n]*cellVolume[n];
    }
}
private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  Field& _varField;
  const  Field& _velocityField;
  Field& _muField;
  Field& _dissipationField;
  Field& _densityField;
  Field& _sourcekField;
  const Field& _gradientField;

};

#endif

