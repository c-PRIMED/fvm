#ifndef _PHONONBOUNDARY_H_
#define _PHONONBOUNDARY_H_

#include "Mesh.h"

#include "NumType.h"
#include "Array.h"
#include "Vector.h"
#include "Field.h"
#include "StorageSite.h"
#include "CRConnectivity.h"
#include "MultiFieldMatrix.h"
#include "CRMatrix.h"
#include "FluxJacobianMatrix.h"
#include "DiagonalMatrix.h"
#include "GeomFields.h"
#include "pmode.h"
#include "Kspace.h"
#include "kvol.h"

template<class X>
class PhononBoundary
{
 public :
  
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;

  typedef Array<int> IntArray;
  
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;  
  typedef Array<VectorT3> VectorT3Array;

  typedef Array<X> XArray;
  typedef Vector<X,3> VectorX3; //also needed for phonons
  typedef Array<VectorX3> VectorX3Array;
  typedef Kspace<X> Xkspace;
  typedef kvol<X> Xkvol;
  typedef pmode<X> Xmode;
 

 PhononBoundary(const StorageSite& faces,
		const Mesh& mesh,
		const GeomFields& geomFields,
		const Xkspace& kspace,
		PhononModelOptions<X>& opts):
  
  _faces(faces),
    _cells(mesh.getCells()),
    _ibType(dynamic_cast<const IntArray&>(geomFields.ibType[_cells])), 
    _faceCells(mesh.getFaceCells(_faces)),
    _areaMagField(geomFields.areaMag),
    _faceAreaMag(dynamic_cast<const TArray&>(_areaMagField[_faces])),
    _areaField(geomFields.area),
    _faceArea(dynamic_cast<const VectorT3Array&>(_areaField[_faces])),
    _options(opts),
    _kspace(kspace)  
  {}

  PhononModelOptions<X>&   getOptions() {return _options;}  

  void applyReflectingWall(int f,const X refl) const
  {
    
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1); ///changed
    X tot_in = 0.;  //total incoming energy (to wall)
    X tot_dk3 = 0.; //total weight of incoming energy
    
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;

    // sum up energy incoming to boundary (from domain)
    int numK=_kspace.getlength();
    for (int k=0;k<numK;k++)
      {

	Xkvol& kv=_kspace.getkvol(k);
	X dk3=kv.getdk3();
	int numM=kv.getmodenum();

	for (int m=0;m<numM;m++) //mode loop beg
	  {
	    
	    Xmode& mode=kv.getmode(m);
	    Field& efield=mode.getfield();
	    VectorT3 vg = mode.getv();     // phonon group velocity
	    XArray& e_val = dynamic_cast< XArray&>(efield[_cells]);  // e"
	    const VectorT3 en = _faceArea[f]/_faceAreaMag[f];  //normal unit vector to face
	    const X vg_dot_en = vg[0]*en[0]+vg[1]*en[1]+vg[2]*en[2];
  
	    if (vg_dot_en > T_Scalar(0.0))
	      {
		tot_in = tot_in + e_val[c0]*dk3;
		tot_dk3=tot_dk3 + dk3;
	      }

	  } //mode loop end	
      }

    const X diff_refl = tot_in/tot_dk3; // value for e" leaving the wall

    // assign values for incoming (upwinded) and outgoing (reflected) to/from wall
    for (int k=0;k<numK;k++)
      {
	
	Xkvol& kv=_kspace.getkvol(k);
	int numM=kv.getmodenum();

	for (int m=0;m<numM;m++) //mode loop beg
	  {
	    
	    Xmode& mode=kv.getmode(m);
	    Field& efield=mode.getfield();
	    VectorT3 vg = mode.getv();     // phonon group velocity
	    XArray& e_val = dynamic_cast<XArray&>(efield[_cells]);  // e"
	    const VectorT3 en = _faceArea[f]/_faceAreaMag[f];  //normal unit vector to face
	    const X vg_dot_en = vg[0]*en[0]+vg[1]*en[1]+vg[2]*en[2];
	    
	    if (vg_dot_en > T_Scalar(0.0))
	      {
		e_val[c1]=e_val[c0];    // upwinded value
	      }
	    else
	      {
		e_val[c1]=diff_refl;  // diffusely reflected value
		// still need to add specularly reflected portion
	      }
	  } //mode loop end	
      }
  }

 void applyReflectingWall(FloatValEvaluator<X>& bRefl) const
  {
    for (int i=0; i<_faces.getCount();i++)
      applyReflectingWall(i,bRefl[i]);
  }


 void applyTemperatureWall(int f,const X Twall) const
  {
    
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1); ///changed
    
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;

    // sum up energy incoming to boundary (from domain)
    int numK=_kspace.getlength();
    X DK3=_kspace.getDK3();
 
    for (int k=0;k<numK;k++)
      {
	Xkvol& kv=_kspace.getkvol(k);
	int numM=kv.getmodenum();
     
	for (int m=0;m<numM;m++) //mode loop beg
	  {
	    X Tref=_options["Tref"];
	    Xmode& mode=kv.getmode(m);
	    Field& efield=mode.getfield();
	    VectorX3 vg = mode.getv();
	    X cp=mode.getcp();
	    XArray& e_val = dynamic_cast< XArray&>(efield[_cells]);  // e"
	    const VectorX3 en = _faceArea[f]/_faceAreaMag[f];  //normal unit vector to face
	    const X vg_dot_en = vg[0]*en[0]+vg[1]*en[1]+vg[2]*en[2];
	    	          
	    if (vg_dot_en > T_Scalar(0.0))
	      {
		e_val[c1]=e_val[c0];
	      }
	    else
	      {
		e_val[c1]=cp*(Twall-Tref)/DK3;
	      }
	    
	  } //mode loop end	
      }
  }

 void applyTemperatureWall(FloatValEvaluator<X>& bTemp) const
  {
    for (int j=0; j<_faces.getCount();j++)
      {
	applyTemperatureWall(j,bTemp[j]);
      }
  }

 /*
 
  void applyZeroGradientBC(int f) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);
    
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    
    const int numDirections = _quadrature.getDirCount();
    
    for (int j=0; j<numDirections; j++)
      {
	Field& fnd = *_dsfPtr.dsf[j];
	XArray& dsf = dynamic_cast< XArray&>(fnd[_cells]);
	dsf[c1]=dsf[c0];
      }
  }
  void applyZeroGradientBC() const
  {
    for (int i=0; i<_faces.getCount();i++)
      applyZeroGradientBC(i);
  } 
    */

  
 protected:

  const StorageSite& _faces;
  const StorageSite& _cells;
  const Array<int>& _ibType;
  const CRConnectivity& _faceCells;
  const Field& _areaMagField;
  const TArray& _faceAreaMag;
  const Field& _areaField;
  const VectorT3Array& _faceArea;
  PhononModelOptions<X>& _options;
  const Xkspace& _kspace;
};
  
#endif
