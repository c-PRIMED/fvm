#ifndef _FLOWMODELSLIPJUMP_H_
#define _FLOWMODELSLIPJUMP_H_

// this file is meant to be included inside the FlowModel::Impl class
// the code is here just because FlowModel_impl.h has grown too big

void slipJumpMomentumBC(const StorageSite& faces,
                        const Mesh& mesh,
                        GenericBCS<VectorT3,DiagTensorT3,T>& gbc,
                        const T accomodationCoefficient)
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

  const T Kn(_options["KnudsenNumber"]);
    
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

      const T coeff = 1 + dn/(accomodationCoefficient*Kn);
      const VectorT3 Vb(Vp/coeff);
      gbc.applyDirichletBC(f,Vb);
  }
}

#endif
