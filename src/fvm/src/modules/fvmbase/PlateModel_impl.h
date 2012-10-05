// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

#include "Mesh.h"

#include "NumType.h"
#include "Array.h"
#include "Field.h"
#include "CRConnectivity.h"
#include "LinearSystem.h"
#include "StorageSite.h"
#include "MultiFieldMatrix.h"
#include "CRMatrix.h"
#include "FluxJacobianMatrix.h"
#include "DiagonalMatrix.h"
#include "GenericBCS.h"
#include "Vector.h"
#include "VectorTranspose.h"
#include "DiffusionDiscretization.h"
#include "TimeDerivativePlateDiscretization.h"
#include "PlateSourceDiscretization.h"
#include "SourceDiscretization.h"
#include "Underrelaxer.h"
#include "AMG.h"
#include "Linearizer.h"
#include "GradientModel.h"
#include "GenericIBDiscretization.h"
#include "StressTensor.h"
#include "SquareTensor.h"


template<class X, class Diag, class OffDiag>
class PlateBCS
{
public:

  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;

  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;
  
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;

  typedef FluxJacobianMatrix<Diag,X> FMatrix;
  typedef DiagonalMatrix<Diag,X> BBMatrix;

  typedef Array<Diag> DiagArray;
  typedef Array<OffDiag> OffDiagArray;
  
  typedef Array<X> XArray;
  typedef Array<VectorT3> VectorT3Array;

  
  PlateBCS(const StorageSite& faces,
             const Mesh& mesh,
             const GeomFields& geomFields,
             Field& varField,
             MultiFieldMatrix& matrix,
             MultiField& xField, MultiField& rField) :
    _faces(faces),
    _cells(mesh.getCells()),
    _faceCells(mesh.getFaceCells(_faces)),
    _varField(varField),
    _xIndex(&_varField,&_cells),
    _dRdX(dynamic_cast<CCMatrix&>(matrix.getMatrix(_xIndex,_xIndex))),
    _assembler(_dRdX.getPairWiseAssembler(_faceCells)),
    _dRdXDiag(_dRdX.getDiag()),
    _x(dynamic_cast<XArray&>(xField[_xIndex])),
    _r(dynamic_cast<XArray&>(rField[_xIndex])),
    _areaMagField(geomFields.areaMag),
    _faceAreaMag(dynamic_cast<const TArray&>(_areaMagField[_faces])),
    _areaField(geomFields.area),
    _faceArea(dynamic_cast<const VectorT3Array&>(_areaField[_faces])),
    _faceCoord(dynamic_cast<const VectorT3Array&>(geomFields.coordinate[_faces])),
    _cellCoord(dynamic_cast<const VectorT3Array&>(geomFields.coordinate[_cells]))
  {}
  
  X applyDirichletBC(int f, const X& bValue) const
  {
    const int c1 = _faceCells(f,1);

    const X fluxB = -_r[c1];
    
    const X dXC1 = bValue - _x[c1];

    _dRdX.eliminateDirichlet(c1,_r,dXC1);
    _x[c1] = bValue;
    _r[c1] = NumTypeTraits<X>::getZero();
   
    return fluxB;
  }

  X applyDirichletBC(const X& bValue) const
  {
    X fluxB(NumTypeTraits<X>::getZero());
    for(int i=0; i<_faces.getCount(); i++)
      fluxB += applyDirichletBC(i,bValue);
    return fluxB;
  }
  
  X applyDirichletBC(const FloatValEvaluator<X>& bValue) const
  {
    X fluxB(NumTypeTraits<X>::getZero());
    for(int i=0; i<_faces.getCount(); i++)
      fluxB += applyDirichletBC(i,bValue[i]);
    return fluxB;
  }
  
  X applyCantileverBC(const int f,
		      const X& specifiedFlux) const
  {
    //const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);
    const X fluxB = -_r[c1];

    VectorT3 dzeta1 = _faceCoord[f]-_cellCoord[c1];

    // the current value of flux and its Jacobians
    X dFlux;
    dFlux[0]= (-(specifiedFlux[0]*_faceArea[f][0] + 
		 specifiedFlux[1]*_faceArea[f][1])*dzeta1[0]) 
      - fluxB[0];
    dFlux[1]= (-(specifiedFlux[0]*_faceArea[f][0] +
                 specifiedFlux[1]*_faceArea[f][1])*dzeta1[1])
      - fluxB[1];
    dFlux[2]= (specifiedFlux[0]*_faceArea[f][0] +
	       specifiedFlux[1]*_faceArea[f][1])
      - fluxB[2];

    // setup the equation for the boundary value; the coefficients
    // are already computed so just need to set the rhs
    _r[c1] = dFlux;

    // _dRdX.eliminate(c1,_r);
    // mark this row as a "boundary" row so that we will update it
    // after the overall system is solved
    _dRdX.setBoundary(c1);
    return fluxB;
  }

  X applyNeumannBC(const int f,
                      const X& specifiedFlux) const
  {
    //const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);
    const X fluxB = -_r[c1];

    // the current value of flux and its Jacobians
    const X dFlux = specifiedFlux*_faceAreaMag[f] - fluxB;

    // setup the equation for the boundary value; the coefficients
    // are already computed so just need to set the rhs
    _r[c1] = dFlux;

    // _dRdX.eliminate(c1,_r);
    // mark this row as a "boundary" row so that we will update it
    // after the overall system is solved
    _dRdX.setBoundary(c1);
    return fluxB;
  }


  X applyNeumannBC(const X& bFlux) const
  {
    X fluxB(NumTypeTraits<X>::getZero());
    
    for(int i=0; i<_faces.getCount(); i++)
      fluxB += applyNeumannBC(i,bFlux);
    return fluxB;
  }
  
  X applyZeroDerivativeBC() const
  {
    X fluxTotal(NumTypeTraits<X>::getZero());
    for(int f=0; f<_faces.getCount(); f++)
    {
        const int c0 = _faceCells(f,0);
        const int c1 = _faceCells(f,1);
        const X fluxB = -_r[c1];
        
        // setup the equation for the boundary value; we want xb to be equal to xc
        
        _r[c1] = _x[c0] - _x[c1];
        _dRdXDiag[c1]  = -NumTypeTraits<Diag>::getUnity();
        this->_assembler.getCoeff10(f)  = NumTypeTraits<OffDiag>::getUnity();

        
        // mark this row as a "boundary" row so that we will update it
        // after the overall system is solved
        _dRdX.setBoundary(c1);
        fluxTotal += fluxB;
    }
    return fluxTotal;
  }

  X applyNeumannBC(const FloatValEvaluator<X>& bFlux) const
  {
    X fluxB(NumTypeTraits<X>::getZero());
    for(int i=0; i<_faces.getCount(); i++)
      fluxB += applyNeumannBC(i,bFlux[i]);
    return fluxB;
  }

  void applyInterfaceBC(const int f) const
  {
    // do nothing
    // the boundary cell could be either c0 or c1 at an interface
    int cb = _faceCells(f,1);
    T_Scalar sign(NumTypeTraits<T_Scalar>::getUnity());
    if (cb < _cells.getSelfCount())
    {
        cb = _faceCells(f,0);
        sign *= -1.0;
    }
    

    _r[cb] = T_Scalar(0);

    if (sign>0)
      _assembler.getCoeff10(f) = NumTypeTraits<OffDiag>::getZero();
    else
      _assembler.getCoeff01(f) = NumTypeTraits<OffDiag>::getZero();
    
  }

  void applyInterfaceBC() const
  {
    // do nothing
  }

  void applySymmetryBC() const
  {
    for(int f=0; f<this->_faces.getCount(); f++)
    {
        const int c0 = this->_faceCells(f,0);
        const int c1 = this->_faceCells(f,1);
        

        const VectorT3 en = this->_faceArea[f]/this->_faceAreaMag[f];
        const T_Scalar xC0_dotn = dot(this->_x[c0],en);
        const X xB = this->_x[c0] - 2.*xC0_dotn * en;

        Diag dxBdxC0(Diag::getZero());
        dxBdxC0(0,0) =  1.0 - 2.*en[0]*en[0];
        dxBdxC0(0,1) =  - 2.*en[0]*en[1];
        dxBdxC0(0,2) =  - 2.*en[0]*en[2];

        dxBdxC0(1,0) =  - 2.*en[1]*en[0];
        dxBdxC0(1,1) =  1.0 - 2.*en[1]*en[1];
        dxBdxC0(1,2) =  - 2.*en[1]*en[2];

        dxBdxC0(2,0) =  - 2.*en[2]*en[0];
        dxBdxC0(2,1) =  - 2.*en[2]*en[1];
        dxBdxC0(2,2) =  1.0 - 2.*en[2]*en[2];
        
        
        const X xc1mxB = xB-this->_x[c1];
        
        // boundary value equation
        // set all neighbour coeffs to zero first and ap to  -1
        this->_dRdX.setDirichlet(c1);

        // dependance on c0
        this->_assembler.getCoeff10(f) = dxBdxC0;
        this->_r[c1] = xc1mxB;

    }
  }


  
protected:
  const StorageSite& _faces;
  const StorageSite& _cells;
  const CRConnectivity& _faceCells;
  const Field& _varField;
  const MultiField::ArrayIndex _xIndex;
  CCMatrix& _dRdX;
  CCAssembler& _assembler;
  DiagArray& _dRdXDiag;
  XArray& _x;
  XArray& _r;
  const Field& _areaMagField;
  const TArray& _faceAreaMag;
  const Field& _areaField;
  const VectorT3Array& _faceArea;
  const VectorT3Array& _faceCoord;
  const VectorT3Array& _cellCoord;
};

template<class T>
class PlateModel<T>::Impl
{
public:
  typedef T T_Scalar;
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef VectorTranspose<T,3> VectorT3T;

  typedef Array<VectorT3> VectorT3Array;
  typedef Vector<T,4> VectorT4;
  typedef Array<VectorT4> VectorT4Array;
  //typedef DiagonalTensor<T,3> DiagTensorT3;
  //typedef SquareTensor<T,3>  SquareTensorT3;
  typedef SquareTensor<T,3>  DiagTensorT3;

  typedef CRMatrix<DiagTensorT3,DiagTensorT3,VectorT3> VVMatrix;
  typedef typename VVMatrix::DiagArray VVDiagArray;


  typedef Gradient<VectorT3> VGradType;
  typedef Array<Gradient<VectorT3> > VGradArray;
  
  Impl(const GeomFields& geomFields,
       PlateFields& plateFields,
       const MeshList& meshes) :
    _meshes(meshes),
    _geomFields(geomFields),
    _plateFields(plateFields),
    _deformationGradientModel(_meshes,_plateFields.deformation,
                           _plateFields.deformationGradient,_geomFields),
    _initialDeformationNorm(),
    _niters(0)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        PlateVC<T> *vc(new PlateVC<T>());
        vc->vcType = "plate";
        _vcMap[mesh.getID()] = vc;
        foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            if ((_bcMap.find(fg.id) == _bcMap.end())&&(fg.groupType != "interior"))
            {
                PlateBC<T> *bc(new PlateBC<T>());
                
                _bcMap[fg.id] = bc;
                if (fg.groupType == "wall") 
                {
                    bc->bcType = "SpecifiedTraction";
                }
                else if (fg.groupType == "interface")
		  {
                    bc->bcType = "Interface";
		  }
		else if (fg.groupType == "symmetry")
		{
		    bc->bcType = "Symmetry";
		}
                else if ((fg.groupType == "velocity-inlet") ||
                         (fg.groupType == "pressure-inlet")) 
                {
                    bc->bcType = "SpecifiedDeformation";
                }
                else if (fg.groupType == "pressure-outlet") 
		{
                    bc->bcType = "SpecifiedForce";
		}
                else
                  throw CException("PlateModel: unknown face group type "
                                   + fg.groupType);
            }
        }
    }
  }

  void init()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const PlateVC<T>& vc = *_vcMap[mesh.getID()];
            
        const StorageSite& cells = mesh.getCells();
	//        const StorageSite& faces = mesh.getFaces();

        shared_ptr<VectorT3Array> sCell(new VectorT3Array(cells.getCountLevel1()));

        VectorT3 initialDeformation;
        initialDeformation[0] = _options["initialXRotation"];
        initialDeformation[1] = _options["initialYRotation"];
        initialDeformation[2] = _options["initialZDeformation"];
        *sCell = initialDeformation;

        _plateFields.deformation.addArray(cells,sCell);
	_plateFields.deformation.syncLocal();

        if (_options.transient)
        {
	    _plateFields.volume0.addArray(cells,
					 dynamic_pointer_cast<ArrayBase>
					 (_geomFields.volume[cells].newCopy()));
            _plateFields.volume0.syncLocal();	
	    				 
	    _plateFields.deformationN1.addArray(cells,
                                            dynamic_pointer_cast<ArrayBase>(sCell->newCopy()));
            _plateFields.deformationN1.syncLocal();
	    					    
            _plateFields.deformationN2.addArray(cells,
					    dynamic_pointer_cast<ArrayBase>(sCell->newCopy()));
            _plateFields.deformationN2.syncLocal();
	    					                
            if (_options.timeDiscretizationOrder > 1)
	      _plateFields.deformationN3.addArray(cells,
					    dynamic_pointer_cast<ArrayBase>(sCell->newCopy()));
              _plateFields.deformationN3.syncLocal();					    
					    
	    if(_options.variableTimeStep)
	    {
		_options.timeStepN1 = _options["timeStep"];
		_options.timeStepN2 = _options["timeStep"];
	    }
        }
        
        shared_ptr<VectorT3Array> stressField(new VectorT3Array(cells.getCountLevel1()));
        stressField->zero();
        _plateFields.stress.addArray(cells,stressField);
	_plateFields.stress.syncLocal();
                
        shared_ptr<VectorT4Array> devStressField(new VectorT4Array((cells.getCountLevel1())*(_options.nz+1)));
        devStressField->zero();
        _plateFields.devStress.addArray(cells,devStressField);
	_plateFields.devStress.syncLocal();

        shared_ptr<TArray> VMStressField(new TArray((cells.getCountLevel1())*(_options.nz+1)));
        VMStressField->zero();
        _plateFields.VMStress.addArray(cells,VMStressField);
	_plateFields.VMStress.syncLocal();

        shared_ptr<TArray> VMStressOutField(new TArray(cells.getCountLevel1()));
        VMStressOutField->zero();
        _plateFields.VMStressOut.addArray(cells,VMStressOutField);
	_plateFields.VMStressOut.syncLocal();

        shared_ptr<VectorT3Array> strainField(new VectorT3Array(cells.getCountLevel1()));
        strainField->zero();
        _plateFields.strain.addArray(cells,strainField);
        _plateFields.strain.syncLocal();        

        shared_ptr<VectorT4Array> plasticStrainField(new VectorT4Array((cells.getCountLevel1())*(_options.nz+1)));
        plasticStrainField->zero();
        _plateFields.plasticStrain.addArray(cells,plasticStrainField);
	_plateFields.plasticStrain.syncLocal();

        shared_ptr<VectorT3Array> plasticStrainOutField(new VectorT3Array(cells.getCountLevel1()));
        plasticStrainOutField->zero();
        _plateFields.plasticStrainOut.addArray(cells,plasticStrainOutField);
	_plateFields.plasticStrainOut.syncLocal();

        shared_ptr<VectorT4Array> plasticStrainN1Field(new VectorT4Array((cells.getCountLevel1())*(_options.nz+1)));
        plasticStrainN1Field->zero();
        _plateFields.plasticStrainN1.addArray(cells,plasticStrainN1Field);
	_plateFields.plasticStrainN1.syncLocal();

        shared_ptr<VectorT3Array> plasticMomentField(new VectorT3Array(cells.getCountLevel1()));
        plasticMomentField->zero();
        _plateFields.plasticMoment.addArray(cells,plasticMomentField); 
	_plateFields.plasticMoment.syncLocal();

        shared_ptr<TArray> rhoCell(new TArray(cells.getCountLevel1()));
        *rhoCell = vc["density"];
        _plateFields.density.addArray(cells,rhoCell);
	_plateFields.density.syncLocal();

        shared_ptr<TArray> ymCell(new TArray(cells.getCountLevel1()));
        *ymCell = vc["ym"];
        _plateFields.ym.addArray(cells,ymCell);
	_plateFields.ym.syncLocal();

        shared_ptr<TArray> nuCell(new TArray(cells.getCountLevel1()));
        *nuCell = vc["nu"];
        _plateFields.nu.addArray(cells,nuCell);
	_plateFields.nu.syncLocal();

        shared_ptr<TArray> forceCell(new TArray(cells.getCountLevel1()));
        forceCell->zero();
        _plateFields.force.addArray(cells,forceCell);
	_plateFields.force.syncLocal();

        shared_ptr<TArray> thicknessCell(new TArray(cells.getCountLevel1()));
        thicknessCell->zero();
        _plateFields.thickness.addArray(cells,thicknessCell);
	_plateFields.thickness.syncLocal();

        shared_ptr<TArray> accelerationCell(new TArray(cells.getCountLevel1()));
        accelerationCell->zero();
        _plateFields.acceleration.addArray(cells,accelerationCell);
	_plateFields.acceleration.syncLocal();
        
        shared_ptr<VectorT3Array> velCell(new VectorT3Array(cells.getCountLevel1()));
        velCell->zero();
        _plateFields.velocity.addArray(cells,velCell);
	_plateFields.velocity.syncLocal();

	//initial temparature gradient array
	shared_ptr<VGradArray> rCell(new VGradArray(cells.getCountLevel1()));
	VGradType residualStress;
	residualStress[0][0] = _options["residualStressXX"];
	residualStress[0][1] = _options["residualStressXY"];
	residualStress[0][2] = _options["residualStressXZ"];
	residualStress[1][0] = _options["residualStressXY"];
	residualStress[1][1] = _options["residualStressYY"];
	residualStress[1][2] = _options["residualStressYZ"];
	residualStress[2][0] = _options["residualStressXZ"];
	residualStress[2][1] = _options["residualStressYZ"];
	residualStress[2][2] = _options["residualStressZZ"];
	*rCell = residualStress;
	_plateFields.residualStress.addArray(cells,rCell);
	_plateFields.residualStress.syncLocal();

        // compute values of deformation flux

        // store deformation flux at interfaces
	/*
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            shared_ptr<VectorT3Array> deformationFlux(new VectorT3Array(faces.getCount()));
            deformationFlux->zero();
            _plateFields.deformationFlux.addArray(faces,deformationFlux);
	}
	*/

	// store deformation flux at boundary faces
        foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
	{
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            shared_ptr<VectorT3Array> deformationFlux(new VectorT3Array(faces.getCount()));
            deformationFlux->zero();
	    if (fg.groupType != "interior")
	    {
	        _plateFields.deformationFlux.addArray(faces,deformationFlux);
	    }
	}

    }

    _niters  =0;
    _initialDeformationNorm = MFRPtr();
  }
  
  PlateBCMap& getBCMap() {return _bcMap;}
  PlateVCMap& getVCMap() {return _vcMap;}
  PlateModelOptions<T>& getOptions() {return _options;}

  void updateTime()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        VectorT3Array& w =
          dynamic_cast<VectorT3Array&>(_plateFields.deformation[cells]);
        VectorT3Array& wN1 =
          dynamic_cast<VectorT3Array&>(_plateFields.deformationN1[cells]);
        VectorT3Array& wN2 =
          dynamic_cast<VectorT3Array&>(_plateFields.deformationN2[cells]);
        if (_options.timeDiscretizationOrder > 1)
	{
            VectorT3Array& wN3 =
              dynamic_cast<VectorT3Array&>(_plateFields.deformationN3[cells]);
            wN3 = wN2;
	}
	wN2 = wN1;
	wN1 = w;
	if(_options.variableTimeStep)
	{
	    if (_options.timeDiscretizationOrder > 1)
	    {
		_options.timeStepN2 = _options.timeStepN1;
	    }
	    _options.timeStepN1 = _options["timeStep"];
	}
        VectorT4Array& pS =
          dynamic_cast<VectorT4Array&>(_plateFields.plasticStrain[cells]);
        VectorT4Array& pSN1 =
          dynamic_cast<VectorT4Array&>(_plateFields.plasticStrainN1[cells]);
        pSN1 = pS;
    }
  }

  
  void initDeformationLinearization(LinearSystem& ls)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex wIndex(&_plateFields.deformation,&cells);

        ls.getX().addArray(wIndex,_plateFields.deformation.getArrayPtr(cells));

        const CRConnectivity& cellCells2 = mesh.getCellCells2();
        
        shared_ptr<Matrix> m(new CRMatrix<DiagTensorT3,DiagTensorT3,VectorT3>(cellCells2));

        ls.getMatrix().addMatrix(wIndex,wIndex,m);

    }
  }

  void linearizeDeformation(LinearSystem& ls)
  {
    _deformationGradientModel.compute();
    DiscrList discretizations;

    const Mesh& mesh = *_meshes[0];
    if(!_options.constForce)
      getMoment(mesh);
    shared_ptr<Discretization>
      sd(new PlateSourceDiscretization<T,DiagTensorT3,DiagTensorT3>
         (_meshes,_geomFields,
          _plateFields.deformation,
          _plateFields.ym,
          _plateFields.nu,
          _plateFields.deformationGradient,
	  _plateFields.residualStress,
          _plateFields.thickness,
          _plateFields.force,
          _options.scf,
          _plateFields.devStress,
          _plateFields.VMStress,
          _plateFields.plasticStrain,
	  _plateFields.plasticStrainOut,
          _plateFields.plasticStrainN1,
          _plateFields.plasticMoment,
          _options.A,
          _options.B,
          _options.m,
          _options.n,
          _options.Sy0,
          _options.nz,
          _options["timeStep"],
          _options.creepModel,
          _options.creep));
      
    //    shared_ptr<Discretization>
    //  bfd(new SourceDiscretization<VectorT3>
    //     (_meshes,_geomFields,
    //    _structureFields.deformation,
    //    _structureFields.bodyForce));

    discretizations.push_back(sd);
    //discretizations.push_back(bfd);
        
    if (_options.transient)
    {
        shared_ptr<Discretization>
          td(new TimeDerivativePlateDiscretization
	     <VectorT3,DiagTensorT3,DiagTensorT3>
             (_meshes,_geomFields,
              _plateFields.deformation,
              _plateFields.deformationN1,
              _plateFields.deformationN2,
              _plateFields.deformationN3,
              _plateFields.density,
	      _plateFields.thickness,
	      _plateFields.acceleration,
	      _plateFields.volume0,
	      _options.variableTimeStep,
              _options["timeStep"],
	      _options.timeStepN1,
	      _options.timeStepN2));
        
        discretizations.push_back(td);
    }
   
    /*
    shared_ptr<Discretization>
      ibm(new GenericIBDiscretization<VectorT3,DiagTensorT3,T>
	  (_meshes,_geomFields,_plateFields.deformation));
      
    discretizations.push_back(ibm);
    */
    Linearizer linearizer;

    linearizer.linearize(discretizations,_meshes,ls.getMatrix(),
                         ls.getX(), ls.getB());
    bool allNeumann = true;
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
	    if (fg.groupType != "interior")
	    {
            const StorageSite& faces = fg.site;
            const int nFaces = faces.getCount();
	    const TArray& faceAreaMag = 
	      dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);
            const PlateBC<T>& bc = *_bcMap[fg.id];
            

            PlateBCS<VectorT3,DiagTensorT3,DiagTensorT3> gbc(faces,mesh,
                                                    _geomFields,
                                                    _plateFields.deformation,
                                                    ls.getMatrix(), ls.getX(), ls.getB());
	    
            VectorT3 fluxB(NumTypeTraits<VectorT3>::getZero());
	    VectorT3 zeroFlux(NumTypeTraits<VectorT3>::getZero());
	    if (bc.bcType == "Clamped")
            {
	        FloatValEvaluator<VectorT3>
		  bDeformation(bc.getVal("specifiedXRotation"),
			       bc.getVal("specifiedYRotation"),
			       bc.getVal("specifiedZDeformation"),
			       faces);
                fluxB = gbc.applyDirichletBC(bDeformation);
                
		allNeumann = false;
            }
            else if (bc.bcType == "Symmetry")
            {
                gbc.applySymmetryBC();
                allNeumann = false;
	    }
            else if (bc.bcType == "ZeroDerivative")
            {
                gbc.applyZeroDerivativeBC();
	    }
            else if (bc.bcType == "SpecifiedTraction")
            {
	        for(int f=0; f<nFaces; f++)
                {
                    
		    fluxB += gbc.applyNeumannBC(f,zeroFlux);
        
		}
	    }
            else if (bc.bcType == "SpecifiedForce")
	    {
	        FloatValEvaluator<VectorT3>
		  bForce(bc.getVal("specifiedXForce"),
			 bc.getVal("specifiedYForce"),
			 bc.getVal("specifiedZForce"),
			 faces);
                for(int f=0; f<nFaces; f++)
		{

                    fluxB += gbc.applyNeumannBC(f,bForce[f]/faceAreaMag[f]);

		}
	    }
            else if (bc.bcType == "SpecifiedShear")
	    {
                FloatValEvaluator<VectorT3>
                  bShear(bc.getVal("specifiedXShear"),
			 bc.getVal("specifiedYShear"),
			 bc.getVal("specifiedZShear"),
			 faces);
                for(int f=0; f<nFaces; f++)
		{
		        
                    fluxB += gbc.applyCantileverBC(f,bShear[f]);

		}
	    }
            else if (bc.bcType != "Interface")
              throw CException(bc.bcType + " not implemented for StructureModel");
            //cout << "force sum for " << fg.id  << " = " << fluxB << endl;  
	    }
        }
	/*
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
	    const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            GenericBCS<VectorT3,DiagTensorT3,DiagTensorT3> gbc(faces,mesh,
                                                    _geomFields,
                                                    _structureFields.deformation,
                                                    _structureFields.deformationFlux,
                                                    ls.getMatrix(), ls.getX(), ls.getB());

            gbc.applyInterfaceBC();
	}
	*/
    }
        
	
#ifdef FVM_PARALLEL
     int count = 1;
     int allNeumannInt = int( allNeumann);
     MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &allNeumannInt, count, MPI::INT, MPI::PROD);
     allNeumann = bool(allNeumannInt);
#endif
	
	
	        
    if(allNeumann && !_options.transient)
    {
        const Mesh& mesh = *_meshes[0];
	const StorageSite& cells = mesh.getCells();
	MultiField& b = ls.getB();
	MultiFieldMatrix& matrix = ls.getMatrix();
	MultiField::ArrayIndex wIndex(&_plateFields.deformation,&cells);
	VectorT3Array& rCell = dynamic_cast<VectorT3Array&>(b[wIndex]);
	VectorT3Array& w = dynamic_cast<VectorT3Array&>
	  (_plateFields.deformation[cells]);
	VVMatrix& vvMatrix =
	  dynamic_cast<VVMatrix&>(matrix.getMatrix(wIndex,wIndex));
	rCell[0] = T(0);
        w[0] = T(0);
        vvMatrix.setDirichlet(0);
    }

#if 0
    shared_ptr<Discretization>
      ud(new Underrelaxer<VectorT3,DiagTensorT3,DiagTensorT3>
         (_meshes,_structureFields.deformation,
          _options["deformationURF"]));
    
    DiscrList discretizations2;
    discretizations2.push_back(ud);

    linearizer.linearize(discretizations2,_meshes,ls.getMatrix(),
                         ls.getX(), ls.getB());
#endif

  }


  MFRPtr solveDeformation()
  {
    LinearSystem ls;
    initDeformationLinearization(ls);
    ls.initAssembly();
    linearizeDeformation(ls);
    ls.initSolve();    
    
    //AMG solver(ls);
    MFRPtr rNorm = _options.getDeformationLinearSolver().solve(ls);
    if (!_initialDeformationNorm) _initialDeformationNorm = rNorm;
        
    _options.getDeformationLinearSolver().cleanup();
    ls.postSolve();
    ls.updateSolution();
    //postPlateSolve(ls);
    if (_options.transient){
      calculatePlateVelocity();
      calculatePlateAcceleration();
    }
    return rNorm;
  }

  void calculatePlateVelocity()
  {
    const int numMeshes = _meshes.size();
    for(int n=0;n<numMeshes;n++)
    {
        const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCountLevel1();
	const VectorT3Array& w =
	  dynamic_cast<const VectorT3Array&>(_plateFields.deformation[cells]);
	const VectorT3Array& wN1 =
	  dynamic_cast<const VectorT3Array&>(_plateFields.deformationN1[cells]);
	VectorT3Array& velocity =
	  dynamic_cast<VectorT3Array&>(_plateFields.velocity[cells]);
	const T timeStep = _options["timeStep"];
	for (int c=0; c<nCells; c++){
	    velocity[c] = (w[c]-wN1[c])/timeStep;
	}
	
    }
  }
	

  void calculatePlateAcceleration()
  {

    T_Scalar two(2.0);
    T_Scalar three(3.0);
    T_Scalar five(5.0);
    T_Scalar four(4.0);
    T_Scalar twelve(12.0);
    T_Scalar dT = _options["timeStep"];
    
    const T_Scalar dT2 = dT*dT;

    const int numMeshes = _meshes.size();
    for(int n=0;n<numMeshes;n++)
    {
        const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCountLevel1();

	TArray& acceleration =
	  dynamic_cast<TArray&>(_plateFields.acceleration[cells]); 
	const VectorT3Array& w =
	  dynamic_cast<const VectorT3Array&>(_plateFields.deformation[cells]);
	const VectorT3Array& wN1 =
	  dynamic_cast<const VectorT3Array&>(_plateFields.deformationN1[cells]);
	const VectorT3Array& wN2 =
	  dynamic_cast<const VectorT3Array&>(_plateFields.deformationN2[cells]);
	const TArray& density = 
	  dynamic_cast<const TArray&> (_plateFields.density[cells]);
	//constant time step
	if(!_options.variableTimeStep)
	  {
	    if (_options.timeDiscretizationOrder > 1)        //second order
	      {
		const VectorT3Array& wN3 =
		  dynamic_cast<const VectorT3Array&>(_plateFields.deformationN3[cells]);
		for(int c=0; c<nCells; c++)
		  {		    
		    const T_Scalar rhobydT2 = density[c]/dT2;
		    acceleration[c] = rhobydT2*(two*w[c][2] - five*wN1[c][2] + four*wN2[c][2]
						- wN3[c][2]);
		  }
	      }
	    else                                              // first order
	      {
		for(int c=0; c<nCells; c++)
		  {
		    const T_Scalar rhobydT2 = density[c]/dT2;
		    acceleration[c] = rhobydT2*(w[c][2]- two*wN1[c][2] + wN2[c][2]);
		  }
	      }
	  }
	//variable time step
	else
	  {
	    T_Scalar dTN1 = _options.timeStepN1;
	    T_Scalar dTN2 = _options.timeStepN2;
	    T_Scalar a = (dT + dTN1)/dT;
	    T_Scalar b = (dT + dTN1 + dTN2)/dT;
	    T_Scalar one(1.0);
	    if (_options.timeDiscretizationOrder > 1)        //second order
	      {	
		const VectorT3Array& wN3 =
		  dynamic_cast<const VectorT3Array&>(_plateFields.deformationN3[cells]);
		T_Scalar c1 = (two*a*b*(pow(a,two)-pow(b,two))+two*b*(pow(b,two)-one)-two*a*(pow(a,two)-one))/
		  (a*b*(a-one)*(b-one)*(a-b));
		T_Scalar c2 = -two*(a+b)/((a-1)*(b-1));
		T_Scalar c3 = -two*(b+one)/(a*(a-b)*(a-one));
		T_Scalar c4 = two*(a+one)/(b*(a-b)*(b-one));		
		for(int c=0; c<nCells; c++)
		  {
		    const T_Scalar rhobydT2 = density[c]/dT2;
		    acceleration[c] = rhobydT2*(c1*w[c][2] + c2*wN1[c][2] + c3*wN2[c][2]
					      + c4*wN3[c][2]);
		  }
	      }
	    else
	      {
		T_Scalar c1 = two/a;
		T_Scalar c2 = -two/(a-one);
		T_Scalar c3 = two/(a*(a-one));
		for(int c=0; c<nCells; c++)
		  {
		    const T_Scalar rhobydT2 = density[c]/dT2;
		    acceleration[c] = rhobydT2*(c1*w[c][2] + c2*wN1[c][2]
						+ c3*wN2[c][2]);
		  }
	      }
	  }
    }
  }

  void postPlateSolve(LinearSystem& ls)
  {
    MultiField& sField = ls.getDelta();

    const int numMeshes = _meshes.size();
    for(int n=0;n<numMeshes;n++)
    {
        const Mesh& mesh = *_meshes[n];

	const StorageSite& cells = mesh.getCells();

	MultiField::ArrayIndex sIndex(&_plateFields.deformation,&cells);
	VectorT3Array& w = dynamic_cast<VectorT3Array&>
	  (_plateFields.deformation[cells]);
        const VectorT3Array& ww = dynamic_cast<const VectorT3Array&>
          (sField[sIndex]);
	//const T deformationURF(_options["deformationURF"]);

	const int nCells = cells.getCountLevel1();
	for(int c=0;c<nCells;c++)
	{
	    w[c] += ww[c];
	}
    }
  }

  bool advance(const int niter)
  {

    for(int n=0; n<niter; n++)
    { 
        MFRPtr dNorm = solveDeformation();

        if (_niters < 5)
        {
            _initialDeformationNorm->setMax(*dNorm);
        }
        
        MFRPtr dNormRatio(dNorm->normalize(*_initialDeformationNorm));
        
        if (_options.printNormalizedResiduals)
          cout << _niters << ": " << *dNormRatio <<  endl;
        else
          cout << _niters << ": " << *dNorm <<  endl;

        _niters++;
        if (*dNormRatio < _options.deformationTolerance)
          return true;
    }
    return false;
  }

  void printBCs()
  {
    foreach(typename PlateBCMap::value_type& pos, _bcMap)
    {
        cout << "Face Group " << pos.first << ":" << endl;
        cout << "    bc type " << pos.second->bcType << endl;
        foreach(typename PlateBC<T>::value_type& vp, *pos.second)
        {
            cout << "   " << vp.first << " "  << vp.second.constant <<  endl;
        }
    }
  }

  void getMoment(const Mesh& mesh)
  {
    const StorageSite& cells = mesh.getCells();

    const int nCells = cells.getCountLevel1();

    shared_ptr<VectorT3Array> momentPtr(new VectorT3Array(nCells));
    momentPtr->zero();
    _plateFields.moment.addArray(cells,momentPtr);
    VectorT3Array& moment = *momentPtr;

    _deformationGradientModel.compute();

    const VGradArray& wGrad =
      dynamic_cast<const VGradArray&>(_plateFields.deformationGradient[cells]);
    
    const TArray& ym = dynamic_cast<const TArray&>(_plateFields.ym[cells]);
    const TArray& nu = dynamic_cast<const TArray&>(_plateFields.nu[cells]);
    
    const TArray& thickness = dynamic_cast<const TArray&>(_plateFields.thickness[cells]);
    
    TArray& VMStress = dynamic_cast<TArray&>(_plateFields.VMStress[cells]);
    TArray& VMStressOut = dynamic_cast<TArray&>(_plateFields.VMStressOut[cells]);
    VectorT3Array& cellStress =
      dynamic_cast<VectorT3Array&>(_plateFields.stress[cells]);      
    VectorT4Array& devStress =
      dynamic_cast<VectorT4Array&>(_plateFields.devStress[cells]);
    VectorT4Array& pg =
      dynamic_cast<VectorT4Array&>(_plateFields.plasticStrain[cells]);    

    VectorT3Array& plasticMoment =
      dynamic_cast<VectorT3Array&>(_plateFields.plasticMoment[cells]);

    VectorT3Array& cellStrain =
      dynamic_cast<VectorT3Array&>(_plateFields.strain[cells]);
    
    const T onethird(1.0/3.0);
    const T one(1.0);
    const T two(2.0);
    const T three(3.0);
    const T six(6.0);
    const T twelve(12.0);
    const T zero(0.0);

    for(int n=0; n<nCells; n++)
    {
        T cellD = ym[n]*pow(thickness[n],three)/(twelve*(one - pow(nu[n],two)));
        const VGradType& wg = wGrad[n];
        VectorT3 stress;
        
        moment[n][0] = cellD*(wg[0][0]+nu[n]*wg[1][1]);
        moment[n][1] = cellD*(wg[1][1]+nu[n]*wg[0][0]);
        moment[n][2] = cellD*((one-nu[n])/two)*(wg[1][0]+wg[0][1]);
        if(_options.creep)
          moment[n] = moment[n] - plasticMoment[n];
	
        cellStress[n][0] = six*moment[n][0]/pow(thickness[n],two);
        cellStress[n][1] = six*moment[n][1]/pow(thickness[n],two);
        cellStress[n][2] = six*moment[n][2]/pow(thickness[n],two);
       
        cellStrain[n][0] = (thickness[n]/two)*wg[0][0];
        cellStrain[n][1] = (thickness[n]/two)*wg[1][1];
        cellStrain[n][2] = (thickness[n]/two)*(wg[1][0]+wg[0][1]);

        T cellE = ym[n]/(one - pow(nu[n],two));

        int nn = n*(_options.nz+1);
        for(int k=0; k<=_options.nz; k++)
	{
            T zz(thickness[n]*(T(k)-T(_options.nz)/T(2))/T(_options.nz));
            stress[0] = (twelve*zz/pow(thickness[n],three))*cellD*(wg[0][0]+nu[n]*wg[1][1])-
              cellE*(pg[nn+k][0]+nu[n]*pg[nn+k][1]);
            stress[1] = (twelve*zz/pow(thickness[n],three))*cellD*(wg[1][1]+nu[n]*wg[0][0])-
              cellE*(pg[nn+k][1]+nu[n]*pg[nn+k][0]);
            stress[2] = (twelve*zz/pow(thickness[n],three))*cellD*((one-nu[n])/two)*(wg[1][0]+wg[0][1])-
              cellE*(one-nu[n])*(pg[nn+k][3]);
            
	    T trace = stress[0] + stress[1];
            devStress[nn+k][0] = stress[0];
            devStress[nn+k][1] = stress[1];
            devStress[nn+k][2] = zero;
            devStress[nn+k][3] = stress[2];
            devStress[nn+k][0] = devStress[nn+k][0] - onethird*trace;
            devStress[nn+k][1] = devStress[nn+k][1] - onethird*trace;
            devStress[nn+k][2] = devStress[nn+k][2] - onethird*trace;
            VMStress[nn+k] = sqrt(pow(stress[0],2.0) + pow(stress[1],2.0) -
				  stress[0]*stress[1] + three*pow(stress[2],2.0));
            //VMStress[nn+k] = sqrt(pow(stress[0],2.0) + three*pow(stress[2],2.0));
	}
	VMStressOut[n] = VMStress[nn+_options.nz];
      }
  }

  map<string,shared_ptr<ArrayBase> >&
  getPersistenceData()
  {
    _persistenceData.clear();
    
    Array<int>* niterArray = new Array<int>(1);
    (*niterArray)[0] = _niters;
    _persistenceData["niters"]=shared_ptr<ArrayBase>(niterArray);

    if (_initialDeformationNorm)
    {
      _persistenceData["initialDeformationNorm"] =
	_initialDeformationNorm->getArrayPtr(_plateFields.deformation);
    }
    else
    {
      Array<Vector<T,3> >* xArray = new Array<Vector<T,3> >(1);
      xArray->zero();
      _persistenceData["initialDeformationNorm"]=shared_ptr<ArrayBase>(xArray);
    }
    return _persistenceData;
  }

  void restart()
  {
    if (_persistenceData.find("niters") != _persistenceData.end())
    {
        shared_ptr<ArrayBase> rp = _persistenceData["niters"];
        ArrayBase& r = *rp;
        Array<int>& niterArray = dynamic_cast<Array<int>& >(r);
        _niters = niterArray[0];
    }
    if (_persistenceData.find("initialDeformationNorm") != _persistenceData.end())
    {
        shared_ptr<ArrayBase>  r = _persistenceData["initialDeformationNorm"];
        _initialDeformationNorm = MFRPtr(new MultiFieldReduction());
        _initialDeformationNorm->addArray(_plateFields.deformation,r);
    }
  }

  #if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC))
  
  void dumpMatrix(const string fileBase)
  {
    LinearSystem ls;
    initDeformationLinearization(ls);
    ls.initAssembly();
    linearizeDeformation(ls);
    ls.initSolve();
    
    MultiFieldMatrix& mfmatrix = ls.getMatrix();
    MultiField& b = ls.getB();

    const Mesh& mesh = *_meshes[0];
    const StorageSite& cells = mesh.getCells();
    
    MultiField::ArrayIndex vIndex(&_plateFields.deformation,&cells);

    VVMatrix& matrix =
      dynamic_cast<VVMatrix&>(mfmatrix.getMatrix(vIndex,vIndex));

    VVDiagArray& diag = matrix.getDiag();
    VVDiagArray& coeff = matrix.getOffDiag();

    VectorT3Array& rCell = dynamic_cast<VectorT3Array&>(b[vIndex]);

    const CRConnectivity& cr = matrix.getConnectivity();

    const Array<int>& row = cr.getRow();
    const Array<int>& col = cr.getCol();
    
    const int nCells = cells.getSelfCount();

    const int dimension = 3;

    const int dim2 = dimension*dimension;
    
    int nFlatRows = dimension*nCells;

    int nFlatCoeffs = nCells*dim2;

    for(int i=0; i<nCells; i++)
      for(int jp=row[i]; jp<row[i+1]; jp++)
      {
          const int j = col[jp];
          if (j<nCells) nFlatCoeffs += dim2;
      }
    
    string matFileName = fileBase + ".mat";
    FILE *matFile = fopen(matFileName.c_str(),"wb");
    
    fprintf(matFile,"%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(matFile,"%d %d %d\n", nFlatRows,nFlatRows,nFlatCoeffs);

    for(int i=0; i<nCells; i++)
    {
        for(int ndr=0; ndr<dimension; ndr++)
          for(int ndc=0; ndc<dimension; ndc++)
          {
              const int nfr = i*dimension + ndr;
              const int nfc = i*dimension + ndc;
              T ap = diag[i](ndr,ndc);
              //if (fabs(ap) < 1.0) ap =0.;
              fprintf(matFile,"%d %d %22.15le\n", nfr+1, nfc+1, ap);
          }
        
        for(int jp=row[i]; jp<row[i+1]; jp++)
        {
            const int j = col[jp];
            if (j<nCells)
            {
                for(int ndr=0; ndr<dimension; ndr++)
                  for(int ndc=0; ndc<dimension; ndc++)
                  {
                      const int nfr = i*dimension + ndr;
                      const int nfc = j*dimension + ndc;
                      T anb = coeff[jp](ndr,ndc);
                      //if (fabs(anb) < 1.0) anb =0.;
                      fprintf(matFile,"%d %d %22.15le\n", nfr+1, nfc+1, anb);
                  }
            }
        }
    }
    
    fclose(matFile);

    string rhsFileName = fileBase + ".rhs";
    FILE *rhsFile = fopen(rhsFileName.c_str(),"wb");
    
    for(int i=0; i<nCells; i++)
      for(int nd=0; nd < dimension; nd++)
            fprintf(rhsFile,"%22.15le\n",-rCell[i][nd]);

    fclose(rhsFile);
  }
#endif


private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  PlateFields& _plateFields;

  PlateBCMap _bcMap;
  PlateVCMap _vcMap;
  
  PlateModelOptions<T> _options;
  GradientModel<VectorT3> _deformationGradientModel;
  
  MFRPtr _initialDeformationNorm;
  int _niters;

  map<string,shared_ptr<ArrayBase> > _persistenceData;
};

template<class T>
PlateModel<T>::PlateModel(const GeomFields& geomFields,
                        PlateFields& plateFields,
                        const MeshList& meshes) :
  Model(meshes),
  _impl(new Impl(geomFields,plateFields,meshes))
{
  logCtor();
}


template<class T>
PlateModel<T>::~PlateModel()
{
  logDtor();
}

template<class T>
void
PlateModel<T>::init()
{
  _impl->init();
}
  
template<class T>
typename PlateModel<T>::PlateBCMap&
PlateModel<T>::getBCMap() {return _impl->getBCMap();}

template<class T>
typename PlateModel<T>::PlateVCMap&
PlateModel<T>::getVCMap() {return _impl->getVCMap();}

template<class T>
PlateModelOptions<T>&
PlateModel<T>::getOptions() {return _impl->getOptions();}


template<class T>
void
PlateModel<T>::printBCs()
{
  _impl->printBCs();
}

template<class T>
bool
PlateModel<T>::advance(const int niter)
{
  return _impl->advance(niter);
}

template<class T>
void
PlateModel<T>::updateTime()
{
  _impl->updateTime();
}

template<class T>
void
PlateModel<T>::getMoment(const Mesh& mesh)
{
  return  _impl->getMoment(mesh);
}

template<class T>
void
PlateModel<T>::dumpMatrix(const string fileBase)
{
  _impl->dumpMatrix(fileBase);
}


template<class T>
map<string,shared_ptr<ArrayBase> >&
PlateModel<T>::getPersistenceData()
{
  return _impl->getPersistenceData();
}

template<class T>
void
PlateModel<T>::restart()
{
  _impl->restart();
}
