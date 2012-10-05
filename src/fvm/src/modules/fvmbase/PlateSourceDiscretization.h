// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

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
  typedef Vector<T,4> VectorT4;
  typedef Array<VectorT4> VectorT4Array;
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
			    const Field& residualStressField,
			    const Field& thicknessField,
			    const Field& forceField,
			    const T& scf,
			    const Field& devStressField,
			    const Field& VMStressField,
			    Field& plasticStrainField,
			    Field& plasticStrainOutField,
			    Field& plasticStrainN1Field,
			    Field& plasticMomentField,
			    const T& A,
			    const T& B,
			    const T& mm,
			    const T& nn,
			    const T& Sy0,
			    const int& nz,
			    const T& timeStep,
			    const int& creepModel,
			    const bool& creep,
			    bool fullLinearization=true)  :
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _ymField(ymField),
    _nuField(nuField),
    _varGradientField(varGradientField),
    _residualStressField(residualStressField),
    _thicknessField(thicknessField),
    _forceField(forceField),
    _scf(scf),
    _devStressField(devStressField),
    _VMStressField(VMStressField),
    _plasticStrainField(plasticStrainField),
    _plasticStrainOutField(plasticStrainOutField),
    _plasticStrainN1Field(plasticStrainN1Field),
    _plasticMomentField(plasticMomentField),
    _A(A),
    _B(B),
    _m(mm),
    _n(nn),
    _Sy0(Sy0),
    _nz(nz),
    _timeStep(timeStep),
    _creepModel(creepModel),
    _creep(creep),
    _fullLinearization(fullLinearization)
   {}

  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();

    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    CCMatrix& matrix =
      dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));

    DiagArray& diag = matrix.getDiag();

    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);

    const TArray& forceCell =
      dynamic_cast<const TArray&>(_forceField[cells]);

    VectorT3Array& rCell = dynamic_cast<VectorT3Array&>(rField[cVarIndex]);
    const VectorT3Array& xCell = dynamic_cast<const VectorT3Array&>(xField[cVarIndex]);

    const VGradArray& residualCell =
      dynamic_cast<const VGradArray&>(_residualStressField[cells]);

    const TArray& ymCell =
      dynamic_cast<const TArray&>(_ymField[cells]);

    const TArray& nuCell =
      dynamic_cast<const TArray&>(_nuField[cells]);

    const TArray& tCell =
      dynamic_cast<const TArray&>(_thicknessField[cells]);
    
    const VectorT4Array& devStress =
      dynamic_cast<const VectorT4Array&>(_devStressField[cells]);

    const TArray& VMStress =
      dynamic_cast<const TArray&>(_VMStressField[cells]);

    VectorT4Array& plasticStrain =
      dynamic_cast<VectorT4Array&>(_plasticStrainField[cells]);

    VectorT3Array& plasticStrainOut =
      dynamic_cast<VectorT3Array&>(_plasticStrainOutField[cells]);

    VectorT4Array& plasticStrainN1 =
      dynamic_cast<VectorT4Array&>(_plasticStrainN1Field[cells]);
    
    VectorT3Array& plasticMoment =
      dynamic_cast<VectorT3Array&>(_plasticMomentField[cells]);

    for(int c=0; c<cells.getSelfCount(); c++)
    {
        rCell[c][2]-=forceCell[c]*cellVolume[c];
    }

    //residual stress
    for(int c=0; c<cells.getSelfCount(); c++)
    {
        const T HSigmazzV = tCell[c]*(residualCell[c])[2][2]*cellVolume[c];
        rCell[c][0] += HSigmazzV*xCell[c][0];
	(diag[c])(0,0) += HSigmazzV;
	rCell[c][1] += HSigmazzV*xCell[c][1];
	(diag[c])(1,1) += HSigmazzV;
    }
    //

    const int nCells = cells.getCountLevel1();
    const T zero(0.0);
    const T half(0.5);
    const T one(1.0);
    const T two(2.0);
    const T three(3.0);
    const T four(4.0);
    const T six(6.0);
    
    if(_creep)
    {
	if (_creepModel==1) 
	{
	    for(int n=0; n<nCells; n++)
	    {
	        int nn = n*(_nz+1);
	        for(int k=0; k<=_nz; k++)
		{
	   
		    T VMPlasticStrain = sqrt(half*
					     (pow(plasticStrain[nn+k][0]-plasticStrain[nn+k][1],2.0) +
					      pow(plasticStrain[nn+k][1]-plasticStrain[nn+k][2],2.0) +
					      pow(plasticStrain[nn+k][2]-plasticStrain[nn+k][0],2.0) +
					      six*(pow(plasticStrain[nn+k][3],2.0))));
		    
		    //T VMPlasticStrain = sqrt(pow(plasticStrain[nn+k][0],2.0));
		    T Sy = _Sy0*(one + _B*pow(VMPlasticStrain,_n));
		    T mult = _A*(pow((VMStress[nn+k]/Sy),_m))/VMStress[nn+k];
		    if(k==_nz/2)
		      mult = zero;
		    for(int i=0;i<4;i++)
		      plasticStrain[nn+k][i] = plasticStrainN1[nn+k][i]+mult*devStress[nn+k][i]*_timeStep;
		}
		plasticStrainOut[n][0] = plasticStrain[nn+_nz][0];
		plasticStrainOut[n][1] = plasticStrain[nn+_nz][1];
		plasticStrainOut[n][2] = plasticStrain[nn+_nz][3];
	    }
	}
	
	for(int n=0; n<nCells; n++)
	{
   	    T var1 = ymCell[n]/(one - pow(nuCell[n],two));
	    T var2 = one - nuCell[n];
	    T var3 = (tCell[n]/T(_nz))/three;
	    T tHalf(0.0);
	    T tempxx(0.0);
	    T tempyy(0.0);
	    T tempxy(0.0);
	    tHalf = tCell[n]/two;
	    int nn = n*(_nz+1);

	    tempxx += (-tHalf)*(plasticStrain[nn][0] + nuCell[n]*plasticStrain[nn][1]);
	    tempyy += (-tHalf)*(plasticStrain[nn][1] + nuCell[n]*plasticStrain[nn][0]);
	    tempxy += (-tHalf)*var2*plasticStrain[nn][3];
	    for(int k=1; k<_nz; k++)
	    {
	        T n1 = four;
		if((k%2)==0)n1 = two;
		T zz(tCell[n]*(T(k)-(T(_nz)/T(2)))/T(_nz));
	        tempxx += T(n1)*zz*(plasticStrain[nn+k][0] + nuCell[n]*plasticStrain[nn+k][1]);
		tempyy += T(n1)*zz*(plasticStrain[nn+k][1] + nuCell[n]*plasticStrain[nn+k][0]); 
		tempxy += T(n1)*zz*var2*plasticStrain[nn+k][3];
	    }
	    tempxx += (tHalf)*(plasticStrain[nn+_nz][0] + nuCell[n]*plasticStrain[nn+_nz][1]);
            tempyy += (tHalf)*(plasticStrain[nn+_nz][1] + nuCell[n]*plasticStrain[nn+_nz][0]);
            tempxy += (tHalf)*var2*plasticStrain[nn+_nz][3];

	    plasticMoment[n][0] = var1*var3*tempxx;
	    plasticMoment[n][1] = var1*var3*tempyy;
	    plasticMoment[n][2] = var1*var3*tempxy;

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

    const VGradArray& residualCell =
      dynamic_cast<const VGradArray&>(_residualStressField[cells]);

    const TArray& ymCell =
      dynamic_cast<const TArray&>(_ymField[cells]);

    const TArray& nuCell =
      dynamic_cast<const TArray&>(_nuField[cells]);

    const TArray& tCell =
      dynamic_cast<const TArray&>(_thicknessField[cells]);

    VectorT3Array& plasticMoment =
      dynamic_cast<VectorT3Array&>(_plasticMomentField[cells]);

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
	VectorT3 dzeta0 = faceCentroid[f]-cellCentroid[c0];
	VectorT3 dzeta1 = faceCentroid[f]-cellCentroid[c1];

        const T diffMetric = faceAreaMag[f]*faceAreaMag[f]/dot(faceArea[f],ds);
        const VectorT3 secondaryCoeff = (faceArea[f]-ds*diffMetric);

        const T dfx0(dzeta0[0]);
        const T dfy0(dzeta0[1]);
        
        const T dfx1(dzeta1[0]);
        const T dfy1(dzeta1[1]);
        
        T vol0 = cellVolume[c0];
        T vol1 = cellVolume[c1];

        T wt0 = vol0/(vol0+vol1);
        T wt1 = vol1/(vol0+vol1);

        T bwt0(wt0);
        T bwt1(wt1);
        if (isBoundary && !isSymmetry)
        {  
            wt0 = T(1.0);
            wt1 = T(0.);
            //bwt0 = T(0.);
            // bwt1 = T(1.);
        }
        
        T faceD(1.0);
	T D0(1.0);
	T D1(1.0);

	T faceG(1.0);
	T G0(1.0);
	T G1(1.0);

	T faceB0(1.0);
	T faceB1(1.0);

	T faceNu(1.0);

	T faceM0(1.0);
	T faceM1(1.0);
	T faceM2(1.0);

	T faceT(1.0);

	D0 = ymCell[c0]*pow(tCell[c0],three)/(twelve*(one - nuCell[c0]*nuCell[c0]));
	D1 = ymCell[c1]*pow(tCell[c1],three)/(twelve*(one - nuCell[c1]*nuCell[c1]));

	G0 = _scf*ymCell[c0]*tCell[c0]/(two*(one+nuCell[c0]));
	G1 = _scf*ymCell[c1]*tCell[c1]/(two*(one+nuCell[c1]));

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

	faceB0 = xCell[c0][0]*bwt0 + xCell[c1][0]*bwt1;
	faceB1 = xCell[c0][1]*bwt0 + xCell[c1][1]*bwt1;

	faceT = tCell[c0]*bwt0 + tCell[c1]*bwt1;

	if(_creep)
	{
	    faceM0 = plasticMoment[c0][0]*bwt0 + plasticMoment[c1][0]*bwt1;
	    faceM1 = plasticMoment[c0][1]*bwt0 + plasticMoment[c1][1]*bwt1;
	    faceM2 = plasticMoment[c0][2]*bwt0 + plasticMoment[c1][2]*bwt1;
	}

        //residual stress contribution//
        const VGradType residualStressF = (residualCell[c0]*wt0 + residualCell[c1]*wt1);
        const T diffRMetric = dot(faceArea[f],(residualStressF*faceArea[f]))/dot(faceArea[f],ds);
        VectorT3 secondaryRCoeff = (residualStressF*faceArea[f]-ds*diffRMetric);
	secondaryRCoeff[0] *= pow(faceT,three)/twelve;
	secondaryRCoeff[1] *= pow(faceT,three)/twelve;
	secondaryRCoeff[2] *= faceT;
        //

	const VGradType gradF = (vGradCell[c0]*wt0 + vGradCell[c1]*wt1);

	VectorT3 source(NumTypeTraits<VectorT3>::getZero());

        // primary w flux, contribution from grad w
        
        T wFlux = faceG*diffMetric*(xCell[c1][2]-xCell[c0][2]);

        // secondary w flux
        wFlux += faceG*(gradF*secondaryCoeff)[2];

        // contribution from beta term

        wFlux += faceG*(faceB0*Af[0]+faceB1*Af[1]);
        
	// primary MxFlux
	T MxFlux = -faceD*diffMetric*(xCell[c1][0]-xCell[c0][0]);

        // secondary MxFlux
        MxFlux -= faceD*(gradF*secondaryCoeff)[0];

	// primary MyFlux
	T MyFlux = -faceD*diffMetric*(xCell[c1][1]-xCell[c0][1]);

        // secondary MyFlux
        MyFlux -= faceD*(gradF*secondaryCoeff)[1];

        // source for cell 0, uses dfx0 and dfy0
        source[0] = -faceD*((faceNu*gradF[1][1])*Af[0]+
                            ((one-faceNu)/two)*(gradF[0][1])*Af[1]-
                            ((one+faceNu)/two)*(gradF[1][0])*Af[1])+
          + dfx0*wFlux + MxFlux;

        source[1] = -faceD*(((one-faceNu)/two)*(gradF[1][0])*Af[0]-
			    ((one+faceNu)/two)*(gradF[0][1])*Af[0]+
			    (faceNu*gradF[0][0])*Af[1])+
	  + dfy0*wFlux + MyFlux ;
    
	source[2] = -wFlux;

	if(_creep)
	{
	    source[0] += faceM0*Af[0] + faceM2*Af[1];
	    source[1] += faceM2*Af[0] + faceM1*Af[1];
	}

        // add flux to the residual of c0 
        rCell[c0] += source;


        // source for cell 1, uses dfx1 and dfy1
	source[0] = -faceD*((faceNu*gradF[1][1])*Af[0]+
			    ((one-faceNu)/two)*(gradF[0][1])*Af[1]-
			    ((one+faceNu)/two)*(gradF[1][0])*Af[1])+
	  + dfx1*wFlux + MxFlux;
        
        source[1] = -faceD*(((one-faceNu)/two)*(gradF[1][0])*Af[0]-
			    ((one+faceNu)/two)*(gradF[0][1])*Af[0]+
			    (faceNu*gradF[0][0])*Af[1])+
	  + dfy1*wFlux + MyFlux ;
        
	source[2] = -wFlux;

        if(_creep)
	{
            source[0] += faceM0*Af[0] + faceM2*Af[1];
            source[1] += faceM2*Af[0] + faceM1*Af[1];
	}

        rCell[c1] -= source;


        // linearization of wflux, primary part

        a00(0,2) += -diffMetric*faceG*dfx0;
        a00(1,2) += -diffMetric*faceG*dfy0;
        a00(2,2) += diffMetric*faceG;

        a01(0,2) += diffMetric*faceG*dfx0;
        a01(1,2) += diffMetric*faceG*dfy0;
        a01(2,2) += -diffMetric*faceG;
        
        a11(0,2) += -diffMetric*faceG*dfx1;
        a11(1,2) += -diffMetric*faceG*dfy1;
        a11(2,2) += diffMetric*faceG;

        a10(0,2) += diffMetric*faceG*dfx1;
        a10(1,2) += diffMetric*faceG*dfy1;
        a10(2,2) += -diffMetric*faceG;
        
        // linearization of wflux, primary part

        a00(0,0) += diffMetric*faceD;
        a00(1,1) += diffMetric*faceD;

        a01(0,0) += -diffMetric*faceD;
        a01(1,1) += -diffMetric*faceD;

        a11(0,0) += diffMetric*faceD;
        a11(1,1) += diffMetric*faceD;

        a10(0,0) += -diffMetric*faceD;
        a10(1,1) += -diffMetric*faceD;   

        // linearization of beta term
	Diag coeffPair;

	coeffPair(0,0)=faceG*dfx0*Af[0];
	coeffPair(0,1)=faceG*dfx0*Af[1];
	coeffPair(0,2)=zero;

	coeffPair(1,0)=faceG*dfy0*Af[0];
	coeffPair(1,1)=faceG*dfy0*Af[1];
	coeffPair(1,2)=zero;

	coeffPair(2,0)=-faceG*Af[0];
	coeffPair(2,1)=-faceG*Af[1];
	coeffPair(2,2)=zero;
			
	a00 += bwt0*coeffPair;
	a01 += bwt1*coeffPair;

	coeffPair(0,0)=faceG*dfx1*Af[0];
	coeffPair(0,1)=faceG*dfx1*Af[1];
	coeffPair(0,2)=zero;

	coeffPair(1,0)=faceG*dfy1*Af[0];
	coeffPair(1,1)=faceG*dfy1*Af[1];
	coeffPair(1,2)=zero;

	coeffPair(2,0)=-faceG*Af[0];
        coeffPair(2,1)=-faceG*Af[1];
	coeffPair(2,2)=zero;
			
        
	a10 -= bwt0*coeffPair;
	a11 -= bwt1*coeffPair;
	
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
		
		coeff(0,0)=-wt0*faceD*(-((one+faceNu)/two)*Af[1]*g_nb[1]);
		coeff(0,1)=-wt0*faceD*(((one-faceNu)/two)*Af[1]*g_nb[0]+faceNu*Af[0]*g_nb[1]);
		coeff(0,2)=wt0*faceG*dfx0*(g_nb[0]*secondaryCoeff[0]+g_nb[1]*secondaryCoeff[1]);

		coeff(1,0)=-wt0*faceD*(((one-faceNu)/two)*Af[0]*g_nb[1]+faceNu*Af[1]*g_nb[0]);
		coeff(1,1)=-wt0*faceD*(-((one+faceNu)/two)*Af[0]*g_nb[0]);
		coeff(1,2)=wt0*faceG*dfy0*(g_nb[0]*secondaryCoeff[0]+g_nb[1]*secondaryCoeff[1]);

		coeff(2,0)=zero;
		coeff(2,1)=zero;
		coeff(2,2)=-wt0*faceG*(g_nb[0]*secondaryCoeff[0]+g_nb[1]*secondaryCoeff[1]);

                //secondary coefficient
                for(int k=0; k<3; k++)
		{
                    coeff(0,0) -= wt0*faceD*secondaryCoeff[k]*g_nb[k];
                    coeff(1,1) -= wt0*faceD*secondaryCoeff[k]*g_nb[k];
                }

		//residual stress
		for(int i=0; i<3; i++)
		{
		  for(int k=0; k<3; k++)
		    coeff(i,i) -= wt0*secondaryRCoeff[k]*g_nb[k];
                }
		//
		
                OffDiag& a0_nb = matrix.getCoeff(c0,nb);

                a0_nb += coeff;
                a00 -= coeff;

                
		coeff(0,0)=-wt0*faceD*(-((one+faceNu)/two)*Af[1]*g_nb[1]);
		coeff(0,1)=-wt0*faceD*(((one-faceNu)/two)*Af[1]*g_nb[0]+faceNu*Af[0]*g_nb[1]);
		coeff(0,2)=wt0*faceG*dfx1*(g_nb[0]*secondaryCoeff[0]+g_nb[1]*secondaryCoeff[1]);

		coeff(1,0)=-wt0*faceD*(((one-faceNu)/two)*Af[0]*g_nb[1]+faceNu*Af[1]*g_nb[0]);
		coeff(1,1)=-wt0*faceD*(-((one+faceNu)/two)*Af[0]*g_nb[0]);
		coeff(1,2)=wt0*faceG*dfy1*(g_nb[0]*secondaryCoeff[0]+g_nb[1]*secondaryCoeff[1]);

		coeff(2,0)=zero;
		coeff(2,1)=zero;
		coeff(2,2)=-wt0*faceG*(g_nb[0]*secondaryCoeff[0]+g_nb[1]*secondaryCoeff[1]);

                //secondary coefficient
                for(int k=0; k<3; k++)
		{
                    coeff(0,0) -= wt0*faceD*secondaryCoeff[k]*g_nb[k];
                    coeff(1,1) -= wt0*faceD*secondaryCoeff[k]*g_nb[k];
                }

		//residual stress
		for(int i=0; i<3; i++)
		{
		    for(int k=0; k<3; k++)
		      coeff(i,i) -= wt0*secondaryRCoeff[k]*g_nb[k];
		}
		//

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
            }


            if (!isBoundary)
            {
                for(int nnb = ccRow[c1]; nnb<ccRow[c1+1]; nnb++)
                {
                    const int nb = ccCol[nnb];
                    const VectorT3& g_nb = vgMatrix.getCoeff(c1,nb);

		    Diag coeff;
		    
		    coeff(0,0)=-wt1*faceD*(-((one+faceNu)/two)*Af[1]*g_nb[1]);
		    coeff(0,1)=-wt1*faceD*(((one-faceNu)/two)*Af[1]*g_nb[0]+faceNu*Af[0]*g_nb[1]);
		    coeff(0,2)=wt1*faceG*dfx1*(g_nb[0]*secondaryCoeff[0]+g_nb[1]*secondaryCoeff[1]);

		    coeff(1,0)=-wt1*faceD*(((one-faceNu)/two)*Af[0]*g_nb[1]+faceNu*Af[1]*g_nb[0]);
		    coeff(1,1)=-wt1*faceD*(-((one+faceNu)/two)*Af[0]*g_nb[0]);
		    coeff(1,2)=wt1*faceG*dfy1*(g_nb[0]*secondaryCoeff[0]+g_nb[1]*secondaryCoeff[1]);

		    coeff(2,0)=zero;
		    coeff(2,1)=zero;
		    coeff(2,2)=-wt1*faceG*(g_nb[0]*secondaryCoeff[0]+g_nb[1]*secondaryCoeff[1]);

                    //secondary coefficient
                    for(int k=0; k<3; k++)
		    {
                        coeff(0,0) -= wt1*faceD*secondaryCoeff[k]*g_nb[k];
                        coeff(1,1) -= wt1*faceD*secondaryCoeff[k]*g_nb[k];
		    }

		    //residual stress
		    for(int i=0; i<3; i++)
		    {
			for(int k=0; k<3; k++)
			  coeff(i,i) -= wt1*secondaryRCoeff[k]*g_nb[k];
		    }
		    //

                    OffDiag& a1_nb = matrix.getCoeff(c1,nb);
                    a1_nb -= coeff;
                    a11 += coeff;

                    
		    coeff(0,0)=-wt1*faceD*(-((one+faceNu)/two)*Af[1]*g_nb[1]);
		    coeff(0,1)=-wt1*faceD*(((one-faceNu)/two)*Af[1]*g_nb[0]+faceNu*Af[0]*g_nb[1]);
		    coeff(0,2)=wt1*faceG*dfx0*(g_nb[0]*secondaryCoeff[0]+g_nb[1]*secondaryCoeff[1]);

		    coeff(1,0)=-wt1*faceD*(((one-faceNu)/two)*Af[0]*g_nb[1]+faceNu*Af[1]*g_nb[0]);
		    coeff(1,1)=-wt1*faceD*(-((one+faceNu)/two)*Af[0]*g_nb[0]);
		    coeff(1,2)=wt1*faceG*dfy0*(g_nb[0]*secondaryCoeff[0]+g_nb[1]*secondaryCoeff[1]);

		    coeff(2,0)=zero;
		    coeff(2,1)=zero;
		    coeff(2,2)=-wt1*faceG*(g_nb[0]*secondaryCoeff[0]+g_nb[1]*secondaryCoeff[1]);

                    //secondary coefficient
                    for(int k=0; k<3; k++)
		    {
                        coeff(0,0) -= wt1*faceD*secondaryCoeff[k]*g_nb[k];
                        coeff(1,1) -= wt1*faceD*secondaryCoeff[k]*g_nb[k];
                    }

		    //residual stress
                    for(int i=0; i<3; i++)
		    {
                        for(int k=0; k<3; k++)
                          coeff(i,i) -= wt1*secondaryRCoeff[k]*g_nb[k];
		    }
		    //

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
	//contribution due to residual stress

        //primary part
        
	VectorT3 residualSource(NumTypeTraits<VectorT3>::getZero());
        residualSource -= diffRMetric*(xCell[c1]-xCell[c0]);
	residualSource[0] *=  pow(faceT,three)/twelve;
	residualSource[1] *=  pow(faceT,three)/twelve;
	residualSource[2] *=  faceT;

        // secondary part

        residualSource -= gradF*secondaryRCoeff;

        // add flux due to residual Stress to the residual of c0 and c1
        rCell[c0] += residualSource;
        rCell[c1] -= residualSource;
     
        const T diffRCoeff = -diffRMetric;

        
        a01(0,0) +=(pow(faceT,three)/twelve)*diffRCoeff;
        a10(0,0) +=(pow(faceT,three)/twelve)*diffRCoeff;

        
        a00(0,0) -=(pow(faceT,three)/twelve)*diffRCoeff;
        a11(0,0) -=(pow(faceT,three)/twelve)*diffRCoeff;

        a01(1,1) +=(pow(faceT,three)/twelve)*diffRCoeff;
        a10(1,1) +=(pow(faceT,three)/twelve)*diffRCoeff;


        a00(1,1) -=(pow(faceT,three)/twelve)*diffRCoeff;
        a11(1,1) -=(pow(faceT,three)/twelve)*diffRCoeff;


        a01(2,2) +=faceT*diffRCoeff;
        a10(2,2) +=faceT*diffRCoeff;


        a00(2,2) -=faceT*diffRCoeff;
        a11(2,2) -=faceT*diffRCoeff;

	//residual stress

    }
  }
    

private:
  const GeomFields& _geomFields;
  Field& _varField;
  const Field& _ymField;
  const Field& _nuField;
  const Field& _varGradientField;
  const Field& _residualStressField;
  const Field& _thicknessField;
  const Field& _forceField;
  const T& _scf;
  const Field& _devStressField;
  const Field& _VMStressField;
  Field& _plasticStrainField;
  Field& _plasticStrainOutField;
  Field& _plasticStrainN1Field;
  Field& _plasticMomentField;
  const T _A;
  const T _B;
  const T _m;
  const T _n;
  const T _Sy0;
  const int _nz;
  const T _timeStep;
  const int _creepModel;
  const bool _creep;
  const bool _fullLinearization;
};

#endif
