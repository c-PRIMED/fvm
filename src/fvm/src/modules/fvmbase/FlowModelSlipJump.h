// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _FLOWMODELSLIPJUMP_H_
#define _FLOWMODELSLIPJUMP_H_

// this file is meant to be included inside the FlowModel::Impl class
// the code is here just because FlowModel_impl.h has grown too big

void slipJumpMomentumBC(const StorageSite& faces,
                        const Mesh& mesh,
                        GenericBCS<VectorT3,DiagTensorT3,T>& gbc,
                        const T accomodationCoefficient,
                        const FloatValEvaluator<VectorT3>& bVelocity)
{
  const StorageSite& cells = mesh.getCells();

  const CRConnectivity& faceCells = mesh.getFaceCells(faces);

  const VectorT3Array& faceArea =
    dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);

  const VectorT3Array& faceCentroid =
    dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[faces]);
    
  const VectorT3Array& cellCentroid =
    dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);
    
  const TArray& faceAreaMag =
    dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);

  const VectorT3Array& V = dynamic_cast<const VectorT3Array&>(_flowFields.velocity[cells]);

  const TArray& p = dynamic_cast<const TArray&>(_flowFields.pressure[cells]);
  const TArray& mu = dynamic_cast<const TArray&>(_flowFields.viscosity[cells]);

  const T opPressure( _options["operatingPressure"]);
  const T opTemperature(_options["operatingTemperature"]);
  const T molWt(_options["molecularWeight"]);
  const T Rgas = 8314.472/molWt;

  const bool incompressible = _options.incompressible;
  
  const int nFaces = faces.getCount();

  for(int f=0; f<nFaces; f++)
  {
      const int c0 = faceCells(f,0);

      // face normal

      const VectorT3 en(faceArea[f]/faceAreaMag[f]);
      
      // normal velocity
      const T Vn = dot(V[c0],en);

      // wall parallel component of cell velocity is Vc - Vn * en
      const VectorT3 Vp(V[c0] - Vn*en);

      const VectorT3 ds(faceCentroid[f]-cellCentroid[c0]);

      // normal distance between face and cell centroid
      const T dn = dot(ds,en);

      
      const T pAbs = incompressible ? opPressure : (p[c0]+opPressure);
      
      
      const T muCell = mu[c0];
      const T MeanFreePath = muCell/pAbs*sqrt(0.5*M_PI*Rgas*opTemperature);

      const T coeff = accomodationCoefficient*MeanFreePath/ (dn+(accomodationCoefficient*MeanFreePath));
      const VectorT3 Vwp(Vp*coeff);

      // specified velocity at the boundary
      const VectorT3 bv = bVelocity[f];
      
      // normal component of the specified boundary velocity
      const VectorT3 Vwn(bv-dot(bv,en)*bv);

      // the velocity at the boundary is the parallel component
      // computed above + the normal from the bc specfication
      const VectorT3 Vw(Vwn + Vwp);
      gbc.applyDirichletBC(f,Vw);
     
  }
}

#endif
