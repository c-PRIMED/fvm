#ifndef _COMETBOUNDARY_H_
#define _COMETBOUNDARY_H_

#include "Mesh.h"
#include <math.h>
#include "NumType.h"
#include "Array.h"
#include "Vector.h"
#include "Field.h"
#include "StorageSite.h"
#include "CRConnectivity.h"
#include "GeomFields.h"
#include "pmode.h"
#include "Kspace.h"
#include "kvol.h"

template<class T>
class COMETBoundary
{
 public :
  
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef Array<int> IntArray;
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;  
  typedef Array<VectorT3> VectorT3Array;
  typedef Kspace<T> Tkspace;
  typedef kvol<T> Tkvol;
  typedef pmode<T> Tmode;
  typedef typename Tmode::Refl_pair Refl_pair;

  COMETBoundary(const StorageSite& faces,
		const Mesh& mesh,
		const GeomFields& geomFields,
		Tkspace& kspace,
		COMETModelOptions<T>& opts,
		const int fg_id):
  _faces(faces),
    _cells(mesh.getCells()),
    _faceCells(mesh.getFaceCells(_faces)),
    _areaMagField(geomFields.areaMag),
    _faceAreaMag(dynamic_cast<const TArray&>(_areaMagField[_faces])),
    _areaField(geomFields.area),
    _faceArea(dynamic_cast<const VectorT3Array&>(_areaField[_faces])),
    _options(opts),
    _kspace(kspace), 
    _fg_id(fg_id),
    _eArray(kspace.geteArray())
    {}

  void applyTemperatureWall(FloatValEvaluator<T>& bTemp) const
  {
    for (int j=0; j<_faces.getCount();j++)
      {
	applyTemperatureWall(j,bTemp[j]);
      }
  }

  void applyTemperatureWall(int f,const T Twall) const
  {
    
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1); 
    
    int numK=_kspace.getlength();
    
    for (int k=0;k<numK;k++)
      {
	Tkvol& kv=_kspace.getkvol(k);
	int numM=kv.getmodenum();
	
	for (int m=0;m<numM;m++)
	  {
	    Tmode& mode=kv.getmode(m);
	    VectorT3 vg = mode.getv();
	    const int index=mode.getIndex()-1;
	    const VectorT3 en = _faceArea[f];
	    const T vg_dot_en = vg[0]*en[0]+vg[1]*en[1]+vg[2]*en[2];
	    const int c0ind=_kspace.getGlobalIndex(c0,index);
	    const int c1ind=_kspace.getGlobalIndex(c1,index);
	    
	    if (vg_dot_en > T_Scalar(0.0))
	      _eArray[c1ind]=_eArray[c0ind];
	    else
	      _eArray[c1ind]=mode.calce0(Twall);
	  }	
      }
  }
  
 protected:
  
  const StorageSite& _faces;
  const StorageSite& _cells;
  const CRConnectivity& _faceCells;
  const Field& _areaMagField;
  const TArray& _faceAreaMag;
  const Field& _areaField;
  const VectorT3Array& _faceArea;
  COMETModelOptions<T>& _options;
  Tkspace& _kspace;
  const int _fg_id;
  TArray& _eArray;
};

#endif
