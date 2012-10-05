// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _WALLDISCRETIZATION_H_
#define _WALLDISCRETIZATION_H_

#include <math.h>
#include "CRMatrix.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "DiagonalMatrix.h"
#include "Gradient.h"
#include "DiagonalTensor.h"

template<class T>
inline T harmonicAvg(const T& x0, const T& x1)
{
  const T sum = x0+x1;
  if (x0+x1 != NumTypeTraits<T>::getZero())
    return 2.0*x0*x1/sum;
  else
    return sum;
}


template<class X, class Diag, class OffDiag>
class WallDiscretization : public Discretization
{
public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
 
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;

  typedef Gradient<X> XGrad;
  typedef Array<int> IntArray;

  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  typedef Array<XGrad> GradArray;

  
 WallDiscretization(const MeshList& meshes,
                     const GeomFields& geomFields,
                     Field& varField,
                     Field& energyField,
                     Field& densityField,
                     Field& wallstressField,
                     Field& parallelvelocityField,
                     Field& tauwallField,
                     Field& diffusivityField,
                     Field& muField,
                     const Field& varGradientField,
                     const T_Scalar thickness=0.0):

    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _energyField(energyField),
    _densityField(densityField),
    _wallstressField(wallstressField),
    _parallelvelocityField(parallelvelocityField),
    _tauwallField(tauwallField),
    _diffusivityField(diffusivityField),
    _muField(muField),
    _varGradientField(varGradientField),
    _thickness(thickness)

 {}

 FlowModelOptions<double>&   getOptions(){return _options;}

 void discretize(const Mesh& mesh,MultiFieldMatrix& mfmatrix,
                 MultiField& xField, MultiField& rField)

 {


    const StorageSite& cells = mesh.getCells();

    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,
                                                             cVarIndex));

    const VectorT3Array& cellCentroid =
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);

    const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
 const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);


  const TArray & rhoCell =
    dynamic_cast<const TArray&>(_densityField[cells]);

  const TArray & kCell =
   dynamic_cast<const TArray&>(_energyField[cells]);

  const TArray & muCell=
    dynamic_cast<const TArray&>(_muField[cells]);


 const VectorT3Array& U =
    dynamic_cast<const VectorT3Array&>(_varField[cells]);

   VectorT3Array& Up =
    dynamic_cast<VectorT3Array&>(_parallelvelocityField[cells]);

  const TArray& diffCell =
      dynamic_cast<const TArray&>(_diffusivityField[cells]);

    const GradArray& xGradCell =
      dynamic_cast<const GradArray&>(_varGradientField[cells]);

   const XArray& xCell = dynamic_cast<const XArray&>(xField[cVarIndex]);
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);



foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
   {
     const FaceGroup& fg = *fgPtr;
     const StorageSite& faces = fg.site;
     if (fg.groupType =="wall")
       {
       const TArray& faceAreaMag =
           dynamic_cast<const TArray &>(_geomFields.areaMag[faces]);

       const VectorT3Array& faceArea=
           dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
 
      const VectorT3Array& faceCentroid =
              dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[faces]);



       VectorT3Array& tauCell =
           dynamic_cast<VectorT3Array&>(_wallstressField[faces]);

 //tau parallel to the wall in direction of flow

        VectorT3Array& TauWallCell =
           dynamic_cast<VectorT3Array&>(_tauwallField[faces]);

        const CRConnectivity& faceCells = mesh.getFaceCells(faces);

        CCAssembler& assembler = matrix.getPairWiseAssembler(faceCells);

        DiagArray& diag = matrix.getDiag();

        const int nFaces = faces.getCount();

         T_Scalar vonk = _options.vk; //vonKarman
         T_Scalar Emp = _options.emp; //Empirical constant
         T_Scalar Cmu = _options.cmu; //Cmu constant of k-e model
         T_Scalar onefourth(0.25);
         T_Scalar eleven(11.225);
    
      for(int f=0; f<nFaces; f++)
           {

             T_Scalar ystar = 0; T_Scalar wallMetric =0;
             const int c0 = faceCells(f,0);
             const int c1 = faceCells(f,1);
             VectorT3 n = faceArea[f]/faceAreaMag[f];

             T_Scalar vol0 = cellVolume[c0];
             T_Scalar vol1 = cellVolume[c1];
             //Computation of Diffusion flux

             VectorT3 ds = cellCentroid[c1]-cellCentroid[c0];

  // for ib faces ignore the solid cell and use the face centroid for diff metric
                if (((ibType[c0] == Mesh::IBTYPE_FLUID)
                     && (ibType[c1] == Mesh::IBTYPE_BOUNDARY)) ||
                    ((ibType[c1] == Mesh::IBTYPE_FLUID)
                     && (ibType[c0] == Mesh::IBTYPE_BOUNDARY)))
                  {
                    if (ibType[c0] == Mesh::IBTYPE_FLUID)
                      {
                        vol1 = 0.;
                        ds = faceCentroid[f]-cellCentroid[c0];
                      }
                    else
                      {
                        vol0 = 0.;
                        ds = cellCentroid[c1]-faceCentroid[f];
                      }
                  }
                //cout << "ds" << ds << endl; 

                T_Scalar faceDiffusivity(1.0);
                if (vol0 == 0.)
                  faceDiffusivity = diffCell[c1];
                else if (vol1 == 0.)
                  faceDiffusivity = diffCell[c0];
                else
                  faceDiffusivity = harmonicAvg(diffCell[c0],diffCell[c1]);
                const T_Scalar diffMetric = faceAreaMag[f]*faceAreaMag[f]/dot(faceArea[f],ds);
                const T_Scalar diffCoeff = faceDiffusivity*diffMetric;
                const VectorT3 secondaryCoeff = faceDiffusivity*(faceArea[f]-ds*diffMetric);

                const XGrad gradF = (xGradCell[c0]*vol0+xGradCell[c1]*vol1)/(vol0+vol1);

                const X dFluxSecondary = gradF*secondaryCoeff;

                const X dFlux = diffCoeff*(xCell[c1]-xCell[c0]) + dFluxSecondary;

        const T_Scalar yp = ds[0]*n[0] + ds[1]*n[1] + ds[2]*n[2];

             ystar = (rhoCell[c0]*sqrt(kCell[c0])*pow(Cmu,onefourth)*yp)/muCell[c0];

             T_Scalar v_dot_n = U[c0][0]*n[0] + U[c0][1]*n[1] + U[c0][2]*n[2];

             for (int i=0; i<3; i++)
              {
                Up[c0][i] = U[c0][i] - v_dot_n*n[i];


                if (ystar > eleven)
                 {
                   wallMetric = (vonk*rhoCell[c0]*pow(Cmu,onefourth)*sqrt(kCell[c0]))/log(Emp*ystar);

                 }

               else
                 {

                   wallMetric = muCell[c0]/yp;

                 }

             tauCell[f][i] = Up[c0][i]*wallMetric;

             T_Scalar tau_dot_n = tauCell[f][0]*n[0]+ tauCell[f][1]*n[1]+ tauCell[f][2]*n[2];

             TauWallCell[f][i] = tauCell[f][i] - tau_dot_n*n[i];

           }
           //  T_Scalar wFlux = TauWallCell[f][0]*faceArea[f][0]+TauWallCell[f][1]*faceArea[f][1]+TauWallCell[f][2]*faceArea[f][2];
            
         // cout << "flux" << TauWallCell[f]-dFlux << endl;
          

            X flux = TauWallCell[f]-dFlux ;

            rCell[c0] += flux;
            rCell[c1] -= flux;


            rCell[c0] += dFlux;
            rCell[c1] -= dFlux;

            //diag[c0] +=diffCoeff;
            //assembler.getCoeff01(f) -=diffCoeff;

            T_Scalar wallCoeff = wallMetric- diffCoeff;

            diag[c1] -=wallCoeff;
            assembler.getCoeff10(f) +=wallCoeff;

         }
       }
     }
}
private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  Field& _varField;
  Field& _energyField;
  Field& _densityField;
  Field& _wallstressField;
  Field& _parallelvelocityField;
  Field& _tauwallField;
  Field& _diffusivityField;
  Field& _muField;
 const Field& _varGradientField;
  FlowModelOptions<double> _options;
  const T_Scalar _thickness;
};

#endif

