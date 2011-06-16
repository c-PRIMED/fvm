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
#include "NumType.h"
#include "Array.h"
#include "GeomFields.h"
#include "FlowFields.h"

template<class T, class Diag, class OffDiag>
class WallDiscretization : public Discretization
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

  WallDiscretization(const MeshList& meshes,
                     const GeomFields& geomFields,
                     Field& varField,
                     Field& energyField,
                     Field& densityField,
                     Field& parallelvelocityField,
                     Field& wallstressField,
                     Field& tauwallField,
                     Field& muField):

    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _energyField(energyField),
    _densityField(densityField),
    _parallelvelocityField(parallelvelocityField),
    _wallstressField(wallstressField),
    _tauwallField(tauwallField),
    _muField(muField)

 {}
 // KeModelOptions<T>&   getOptions() {return _options;}
  FlowModelOptions<T>&   getOptions(){return _options;}

 void discretize(const Mesh& mesh,MultiFieldMatrix& mfmatrix,
                 MultiField& xField, MultiField& rField)

 {

  const StorageSite& cells = mesh.getCells();

  const MultiField::ArrayIndex cVarIndex(&_varField,&cells);

 const VectorT3Array& cellCentroid =
    dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);

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


  TArray& rCell =
    dynamic_cast<TArray&>(rField[cVarIndex]);

 //const TArray& xCell = dynamic_cast<const TArray&>(xField[cVarIndex]);
 

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

         VectorT3Array& tau =
           dynamic_cast<VectorT3Array&>(_wallstressField[faces]);

 //tau parallel to the wall in direction of flow

        VectorT3Array& TauWall =
           dynamic_cast<VectorT3Array&>(_tauwallField[faces]);
 
        const CRConnectivity& faceCells = mesh.getFaceCells(faces);
        
        CCMatrix& matrix = 
           dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,
                                                             cVarIndex)); 
        CCAssembler& assembler = matrix.getPairWiseAssembler(faceCells);
        
        DiagArray& diag = matrix.getDiag();


        const int nFaces = faces.getCount();


         T vonk = _options.vk; //vonKarman
         T Emp = _options.emp; //Empirical constant
         T Cmu = _options.cmu; //Cmu constant of k-e model
         T onefourth(0.25);
         T eleven(11.225);
 
         for(int f=0; f<nFaces; f++)
           {

             T ystar = 0; T diffMetric =0; 

             VectorT3 n = faceArea[f]/faceAreaMag[f];

             const int c0 = faceCells(f,0);
             const int c1 = faceCells(f,1);

            // Diag& a00 = diag[c0];
             Diag& dRdX1 = diag[c1];
            // OffDiag& a01 = assembler.getCoeff01(f);
             OffDiag& dRdX0 = assembler.getCoeff10(f);


             ystar = (rhoCell[c0]*sqrt(kCell[c0])*pow(Cmu,onefourth)*cellCentroid[c0][1])/muCell[c0];

             T v_dot_n = U[c0][0]*n[0] + U[c0][1]*n[1] + U[c0][2]*n[2];
 
             for (int i=0; i<3; i++)
              { 
                Up[c0][i] = U[c0][i] - v_dot_n*n[i];
              
   
                if (ystar > eleven)
                 {
                   diffMetric = (vonk*rhoCell[c0]*pow(Cmu,onefourth)*sqrt(kCell[c0]))/log(Emp*ystar);

                 }  

               else
                 {
          
                   diffMetric = muCell[c0]/cellCentroid[c0][1];
                 
                 }

             tau[f][i] = Up[c0][i]*diffMetric;

             T tau_dot_n = tau[f][0]*n[0]+ tau[f][1]*n[1]+ tau[f][2]*n[2];
    
             TauWall[f][i] = tau[f][i] - tau_dot_n*n[i];
    
           }
             T force = TauWall[f][0]*faceArea[f][0]+TauWall[f][1]*faceArea[f][1]+TauWall[f][2]*faceArea[f][2];

            rCell[c0] += force;
            rCell[c1] -= force;

            dRdX1 +=diffMetric;
            dRdX0 +=diffMetric;
 
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
  Field& _parallelvelocityField;
  Field& _wallstressField;
  Field& _tauwallField;
  Field& _muField;
  FlowModelOptions<T> _options;
};

#endif










