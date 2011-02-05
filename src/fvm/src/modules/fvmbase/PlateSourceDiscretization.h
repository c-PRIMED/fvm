#ifndef _PLATESOURCEDISCRETIZATION_H_
#define _PLATESOURCEDISCRETIZATION_H_

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

template<class T, class Diag, class OffDiag>
class PlateSourceDiscretization : public Discretization
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

  typedef GradientModel<VectorT3> VGradModelType;
  typedef typename VGradModelType::GradMatrixType VGradMatrix;
  
  PlateSourceDiscretization(const MeshList& meshes,
				const GeomFields& geomFields,
				Field& varField,
				const Field& ymField,
				const Field& nuField,
			        const Field& varGradientField,
			        const Field& thicknessField,
			        const Field& forceField,
			        const T& scf,                                
                                bool fullLinearization=true)  :
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _ymField(ymField),
    _nuField(nuField),
    _varGradientField(varGradientField),
    _thicknessField(thicknessField),
    _forceField(forceField),
    _scf(scf),
    _fullLinearization(fullLinearization)
   {}

  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();

    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);

    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);

    const TArray& forceCell =
      dynamic_cast<const TArray&>(_forceField[cells]);

    VectorT3Array& rCell = dynamic_cast<VectorT3Array&>(rField[cVarIndex]);

    for(int c=0; c<cells.getSelfCount(); c++)
    {
        rCell[c][2]-=forceCell[c]*cellVolume[c];
	if(c==0)
	{
	    cout<<"rCell = "<<rCell[c][2]<<endl;
	}
    }

    const StorageSite& iFaces = mesh.getInteriorFaceGroup().site;
    
    discretizeFaces(mesh, iFaces, mfmatrix, xField, rField, false, false);

    /*
    foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
    {
        const FaceGroup& fg = *fgPtr;
        const StorageSite& faces = fg.site;
        discretizeFaces(mesh, faces, mfmatrix, xField, rField, false, false);
    }
    */
        
    // boundaries
    foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
    {
        const FaceGroup& fg = *fgPtr;
        const StorageSite& faces = fg.site;
	if (fg.groupType!="interior")
	{
	    discretizeFaces(mesh, faces, mfmatrix, xField, rField,
			    fg.groupType!="interface",
			    fg.groupType=="symmetry");
	}
    }
  }

                          
  void discretizeFaces(const Mesh& mesh, const StorageSite& faces,
                       MultiFieldMatrix& mfmatrix,
                       MultiField& xField, MultiField& rField,
                       const bool isBoundary, const bool isSymmetry)
  {
    const StorageSite& cells = mesh.getCells();

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

    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

    VectorT3Array& rCell = dynamic_cast<VectorT3Array&>(rField[cVarIndex]);
    const VectorT3Array& xCell = dynamic_cast<const VectorT3Array&>(xField[cVarIndex]);

    const VGradArray& vGradCell =
      dynamic_cast<const VGradArray&>(_varGradientField[cells]);

    const TArray& ymCell =
      dynamic_cast<const TArray&>(_ymField[cells]);

    const TArray& nuCell =
      dynamic_cast<const TArray&>(_nuField[cells]);

    const TArray& tCell =
      dynamic_cast<const TArray&>(_thicknessField[cells]);

    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,
                                                             cVarIndex));
    CCAssembler& assembler = matrix.getPairWiseAssembler(faceCells);
    DiagArray& diag = matrix.getDiag();

    const int nFaces = faces.getCount();

    const VGradMatrix& vgMatrix = VGradModelType::getGradientMatrix(mesh,_geomFields);
    const CRConnectivity& cellCells = mesh.getCellCells();
    const Array<int>& ccRow = cellCells.getRow();
    const Array<int>& ccCol = cellCells.getCol();

    const int nInteriorCells = cells.getSelfCount();

    const T zero(0.0);
    const T one(1.0);
    const T two(2.0);
    const T three(3.0);
    const T twelve(12.0);

    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);

	const VectorT3& Af = faceArea[f];
        const VectorT3 en = Af/faceAreaMag[f];
        
        VectorT3 ds=cellCentroid[c1]-cellCentroid[c0];
	VectorT3 dzeta = faceCentroid[f]-cellCentroid[c0];

        T vol0 = cellVolume[c0];
        T vol1 = cellVolume[c1];

        T wt0 = vol0/(vol0+vol1);
        T wt1 = vol1/(vol0+vol1);

        if (isBoundary && !isSymmetry)
        {  
            wt0 = T(1.0);
            wt1 = T(0.);
        }
        
        T faceD(1.0);
	T D0(1.0);
	T D1(1.0);

	T faceG(1.0);
	T G0(1.0);
	T G1(1.0);

	T faceB0(1.0);
	T faceB1(1.0);

	T df0(0.0);
	T df1(0.0);
	
	T faceNu(1.0);

	D0 = ymCell[c0]*pow(tCell[c0],three)/(twelve*(one - nuCell[c0]*nuCell[c0]));
	D1 = ymCell[c1]*pow(tCell[c1],three)/(twelve*(one - nuCell[c1]*nuCell[c1]));

	G0 = _scf*ymCell[c0]*tCell[c0]/(two*(one+nuCell[c0]));
	G1 = _scf*ymCell[c1]*tCell[c1]/(two*(one+nuCell[c1]));

	df0 = dzeta[0];
	df1 = dzeta[1];
	
	if(f==100)
	{
	    cout<<"df0 ="<<df0<<" df1 = "<<df1<<endl;
	    cout<<"D ="<<D0<<" and "<<D1<<endl;
	    cout<<"D ="<<G0<<" and "<<G1<<endl;
	}

        Diag& a00 = diag[c0];
        Diag& a11 = diag[c1];
        OffDiag& a01 = assembler.getCoeff01(f);
        OffDiag& a10 = assembler.getCoeff10(f);

        if (vol0 == 0.)
       	{
	    faceD = D1;
	    faceG = G1;
	    faceNu = nuCell[c1];
	}
        else if (vol1 == 0.)
	{
	    faceD = D0;
	    faceG = G0;
	    faceNu = nuCell[c0];
	}
        else
	{
            faceD = harmonicAverage(D0,D1);
	    faceG = harmonicAverage(G0,G1);
	    faceNu = harmonicAverage(nuCell[c0],nuCell[c1]);
	}

        faceD = D0*wt0 + D1*wt1;
	faceG = G0*wt0 + G1*wt1;
	faceNu = nuCell[c0]*wt0 + nuCell[c1]*wt1;

	faceB0 = xCell[c0][0]*wt0 + xCell[c1][0]*wt1;
	faceB1 = xCell[c0][1]*wt0 + xCell[c1][1]*wt1;

        if(f==100)
	{
            cout<<faceD<<" "<<faceG<<" "<<faceNu<<endl;
	}

	const VGradType gradF = (vGradCell[c0]*wt0 + vGradCell[c1]*wt1);

	VectorT3 source(NumTypeTraits<VectorT3>::getZero());

        // mu*grad U ^ T + lambda * div U I
	source[0] = -faceD*((gradF[0][0]+faceNu*gradF[1][1])*Af[0]+
			    ((one-faceNu)/two)*(gradF[0][1]+gradF[1][0])*Af[1])+
	  faceG*df0*((gradF[2][0]+faceB0)*Af[0]+(gradF[2][1]+faceB1)*Af[1]);
        
        source[1] = -faceD*(((one-faceNu)/two)*(gradF[0][1]+gradF[1][0])*Af[0]+
			    (gradF[1][1]+faceNu*gradF[0][0])*Af[1])+
	  faceG*df1*((gradF[2][0]+faceB0)*Af[0]+(gradF[2][1]+faceB1)*Af[1]);
        
	source[2] = -faceG*((gradF[2][0]+faceB0)*Af[0]+(gradF[2][1]+faceB1)*Af[1]);

	Diag coeffPair;

	coeffPair(0,0)=faceG*df0*Af[0];
	coeffPair(0,1)=faceG*df0*Af[1];
	coeffPair(0,2)=zero;

	coeffPair(1,0)=faceG*df1*Af[0];
	coeffPair(1,1)=faceG*df1*Af[1];
	coeffPair(1,2)=zero;

	coeffPair(2,0)=-faceG*Af[0];
	coeffPair(2,1)=-faceG*Af[1];
	coeffPair(2,2)=zero;
			
	a00 += wt0*coeffPair;
	a01 += wt1*coeffPair;
	a10 -= wt0*coeffPair;
	a11 -= wt1*coeffPair;
	
        VectorT3 s0(NumTypeTraits<VectorT3>::getZero());
        VectorT3 s1(NumTypeTraits<VectorT3>::getZero());
        if (_fullLinearization)
        {
            // loop over level 0 neighbors
            for(int nnb = ccRow[c0]; nnb<ccRow[c0+1]; nnb++)
            {
                const int nb = ccCol[nnb];

                // get the coefficient from the gradient matrix that uses level 2 connectivity
                // getting a copy rather than reference since we might modify it
                VectorT3 g_nb = vgMatrix.getCoeff(c0,nb);
#if 1
                if (isSymmetry)
                {
                    const T gnb_dot_nx2 = T(2.0)*dot(en,g_nb);
                    g_nb = gnb_dot_nx2*en;
                }
#endif
                Diag coeff;
		
		coeff(0,0)=-wt0*faceD*(Af[0]*g_nb[0]+((one-faceNu)/two)*Af[1]*g_nb[1]);
		coeff(0,1)=-wt0*faceD*(((one-faceNu)/two)*Af[1]*g_nb[0]+faceNu*Af[0]*g_nb[1]);
		coeff(0,2)=wt0*faceG*df0*(g_nb[0]*Af[0]+g_nb[1]*Af[1]);

		coeff(1,0)=-wt0*faceD*(((one-faceNu)/two)*Af[0]*g_nb[1]+faceNu*Af[1]*g_nb[0]);
		coeff(1,1)=-wt0*faceD*(Af[1]*g_nb[1]+((one-faceNu)/two)*Af[0]*g_nb[0]);
		coeff(1,2)=wt0*faceG*df1*(g_nb[0]*Af[0]+g_nb[1]*Af[1]);

		coeff(2,0)=zero;
		coeff(2,1)=zero;
		coeff(2,2)=-wt0*faceG*(g_nb[0]*Af[0]+g_nb[1]*Af[1]);
		
#if 0
                if (isSymmetry)
                {
                    SquareTensor<T,3> R;
                    
                    for(int i=0;i<3;i++)
                      for(int j=0; j<3; j++)
                      {
                          if (i==j)
                            R(i,j) = 1.0 - 2*en[i]*en[j];
                          else
                            R(i,j) = - 2*en[i]*en[j];
                      }

                    Diag coeff1(R*coeff*R);
                    coeff += coeff1;
                      
                }
#endif
                OffDiag& a0_nb = matrix.getCoeff(c0,nb);

                a0_nb += coeff;
                a00 -= coeff;
                a10 += coeff;

                if (c1 != nb)
                {
                    OffDiag& a1_nb = matrix.getCoeff(c1,nb);
                    a1_nb -= coeff;
                }
                else
		{
                    a11 -= coeff;
		}
		if(f==100)
		  {
		    cout<<"coeff "<<coeff<<endl;
		  }
            }


            if (!isBoundary)
            {
                for(int nnb = ccRow[c1]; nnb<ccRow[c1+1]; nnb++)
                {
                    const int nb = ccCol[nnb];
                    const VectorT3& g_nb = vgMatrix.getCoeff(c1,nb);

		    Diag coeff;
		    
		    coeff(0,0)=-wt1*faceD*(Af[0]*g_nb[0]+((one-faceNu)/two)*Af[1]*g_nb[1]);
		    coeff(0,1)=-wt1*faceD*(((one-faceNu)/two)*Af[1]*g_nb[0]+faceNu*Af[0]*g_nb[1]);
		    coeff(0,2)=wt1*faceG*df0*(g_nb[0]*Af[0]+g_nb[1]*Af[1]);

		    coeff(1,0)=-wt1*faceD*(((one-faceNu)/two)*Af[0]*g_nb[1]+faceNu*Af[1]*g_nb[0]);
		    coeff(1,1)=-wt1*faceD*(Af[1]*g_nb[1]+((one-faceNu)/two)*Af[0]*g_nb[0]);
		    coeff(1,2)=wt1*faceG*df1*(g_nb[0]*Af[0]+g_nb[1]*Af[1]);

		    coeff(2,0)=zero;
		    coeff(2,1)=zero;
		    coeff(2,2)=-wt1*faceG*(g_nb[0]*Af[0]+g_nb[1]*Af[1]);

                    OffDiag& a1_nb = matrix.getCoeff(c1,nb);
                    a1_nb -= coeff;
                    a11 += coeff;
                    a01 -= coeff;

                    if (c0 != nb)
                    {
                        OffDiag& a0_nb = matrix.getCoeff(c0,nb);
                        
                        a0_nb += coeff;
                    }
                    else
                      a00 += coeff;
                }
            }
        }
	if(c0==0)
	{

	    cout<<"a00 = "<<a00<<endl;
	    cout<<"a01 = "<<a01<<endl;
	    cout<<"a10 = "<<a10<<endl;
	    cout<<"a11 = "<<a11<<endl;
	    cout<<endl;
	}

        // add flux to the residual of c0 and c1
        rCell[c0] += source;
	rCell[c1] -= source;

    }
  }
    

private:
  const GeomFields& _geomFields;
  Field& _varField;
  const Field& _ymField;
  const Field& _nuField;
  const Field& _varGradientField;
  const Field& _thicknessField;
  const Field& _forceField;
  const T& _scf;
  const bool _fullLinearization;
};

#endif
