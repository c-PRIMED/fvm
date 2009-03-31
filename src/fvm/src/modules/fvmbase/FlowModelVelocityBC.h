#ifndef _FLOWMODELVELOCITYBC_H_
#define _FLOWMODELVELOCITYBC_H_

// this file is meant to be included inside the FlowModel::Impl class
// the code is here just because FlowModel_impl.h has grown too big

  T fixedFluxContinuityBC(const StorageSite& faces,
                          const Mesh& mesh,
                          MultiFieldMatrix& matrix,
                          MultiField& xField,
                          MultiField& rField,
                          const FlowBC<T>& bc)
  {
    const StorageSite& cells = mesh.getCells();

    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

    MultiField::ArrayIndex mfIndex(&_flowFields.massFlux,&faces);
    MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
    MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);

    PPMatrix& ppMatrix =
      dynamic_cast<PPMatrix&>(matrix.getMatrix(pIndex,pIndex));

    FMatrix& dFluxdP = dynamic_cast<FMatrix&>(matrix.getMatrix(mfIndex,pIndex));

    PPAssembler& ppAssembler = ppMatrix.getPairWiseAssembler(faceCells);
    PPDiagArray& ppDiag = ppMatrix.getDiag();

    const VectorT3Array& faceArea =
      dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);

    TArray& rCell = dynamic_cast<TArray&>(rField[pIndex]);
    TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);
    const TArray& density = dynamic_cast<const TArray&>(_flowFields.density[cells]);

    VectorT3 bVelocity;
    bVelocity[0] = bc["specifiedXVelocity"];
    bVelocity[1] = bc["specifiedYVelocity"];
    bVelocity[2] = bc["specifiedZVelocity"];

    const int nFaces = faces.getCount();

    T netFlux(0.);
    
    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);

        massFlux[f] = density[c0]*dot(bVelocity,faceArea[f]);

        rCell[c0] -= massFlux[f];

        netFlux += massFlux[f];
        ppAssembler.getCoeff01(f) =0;
        ppAssembler.getCoeff10(f) =1;
        ppDiag[c1] = -1;
        rCell[c1] = 0.;
        ppMatrix.setBoundary(c1);

        dFluxdP.setCoeffL(f,T(0.));
        dFluxdP.setCoeffR(f,T(0.));
        
    }
#if 1
    if (matrix.hasMatrix(vIndex,pIndex))
    {
        VPMatrix& vpMatrix =
          dynamic_cast<VPMatrix&>(matrix.getMatrix(vIndex,pIndex));
        
        VPAssembler& vpAssembler = vpMatrix.getPairWiseAssembler(faceCells);
        VPDiagArray& vpDiag = vpMatrix.getDiag();

        for(int f=0; f<nFaces; f++)
        {
            const int c0 = faceCells(f,0);
            vpDiag[c0] += vpAssembler.getCoeff01(f);
            vpAssembler.getCoeff01(f) = 0;
        }
    }
#endif
    return netFlux;
  }

#endif
