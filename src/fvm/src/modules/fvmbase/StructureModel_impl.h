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
#include "TimeDerivativeStructureDiscretization.h"
#include "StructureSourceDiscretization.h"
#include "StructurePlasticDiscretization.h"
#include "SourceDiscretization.h"
#include "Underrelaxer.h"
#include "AMG.h"
#include "Linearizer.h"
#include "GradientModel.h"
#include "GenericIBDiscretization.h"
#include "StressTensor.h"
#include "SquareTensor.h"



template<class X, class Diag, class OffDiag>
class StructureBCS
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

  
  StructureBCS(const StorageSite& faces,
               const Mesh& mesh,
               const GeomFields& geomFields,
               Field& varField,
               MultiFieldMatrix& matrix,
               MultiField& xField,
               MultiField& rField,
               const bool explicitMode
               ) :
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
    _explicitMode(explicitMode)
  {}
  
  X applyDirichletBC(int f, const X& bValue) const
  {
    const int c1 = _faceCells(f,1);

    const X fluxB = -_r[c1];
    
    const X dXC1 = bValue - _x[c1];

    _dRdX.eliminateDirichlet(c1,_r,dXC1, _explicitMode);
    _x[c1] = bValue;
    _r[c1] = NumTypeTraits<X>::getZero();
    if (!_explicitMode)
      _dRdX.setDirichlet(c1);
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
    for(int i=0; i<_faces.getCount(); i++)
      applyInterfaceBC(i);
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
  const bool _explicitMode;
};

template<class T>
class StructureModel<T>::Impl
{
public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef VectorTranspose<T,3> VectorT3T;

  typedef Array<VectorT3> VectorT3Array;
  //typedef DiagonalTensor<T,3> DiagTensorT3;
  //typedef SquareTensor<T,3>  SquareTensorT3;
  typedef SquareTensor<T,3>  DiagTensorT3;

  typedef CRMatrix<DiagTensorT3,DiagTensorT3,VectorT3> VVMatrix;
  typedef typename VVMatrix::DiagArray VVDiagArray;


  typedef Gradient<VectorT3> VGradType;
  typedef Array<Gradient<VectorT3> > VGradArray;
  
  Impl(const GeomFields& geomFields,
       StructureFields& structureFields,
       const MeshList& meshes) :
    _meshes(meshes),
    _geomFields(geomFields),
    _structureFields(structureFields),
    _deformationGradientModel(_meshes,_structureFields.deformation,
                           _structureFields.deformationGradient,_geomFields),
    _initialDeformationNorm(),
    _niters(0)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        StructureVC<T> *vc(new StructureVC<T>());
        vc->vcType = "structure";
        _vcMap[mesh.getID()] = vc;
        foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            if ((_bcMap.find(fg.id) == _bcMap.end())&&(fg.groupType != "interior"))
            {
                StructureBC<T> *bc(new StructureBC<T>());
                
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
                  throw CException("StructuralModel: unknown face group type "
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

        const StructureVC<T>& vc = *_vcMap[mesh.getID()];
            
        const StorageSite& cells = mesh.getCells();
	//        const StorageSite& faces = mesh.getFaces();

        shared_ptr<VectorT3Array> sCell(new VectorT3Array(cells.getCountLevel1()));

        VectorT3 initialDeformation;
        initialDeformation[0] = _options["initialXDeformation"];
        initialDeformation[1] = _options["initialYDeformation"];
        initialDeformation[2] = _options["initialZDeformation"];
        *sCell = initialDeformation;

        _structureFields.deformation.addArray(cells,sCell);

        if (_options.transient)
        {
	    _structureFields.volume0.addArray(cells,
					 dynamic_pointer_cast<ArrayBase>
					 (_geomFields.volume[cells].newCopy()));
	    _structureFields.deformationN1.addArray(cells,
                                            dynamic_pointer_cast<ArrayBase>(sCell->newCopy()));
            _structureFields.deformationN2.addArray(cells,
					    dynamic_pointer_cast<ArrayBase>(sCell->newCopy()));            
            if (_options.timeDiscretizationOrder > 1)
	      _structureFields.deformationN3.addArray(cells,
						      dynamic_pointer_cast<ArrayBase>(sCell->newCopy()));
            if(_options.variableTimeStep)
            {
                _options.timeStepN1 = _options["timeStep"];
                _options.timeStepN2 = _options["timeStep"];
            }
        }
        
        if (_options.creep)
	{
            _structureFields.elasticDeformation.addArray(cells,
                                                         dynamic_pointer_cast<ArrayBase>(sCell->newCopy()));
            shared_ptr<VGradArray> devStressField(new VGradArray(cells.getCountLevel1()));
            devStressField->zero();
            _structureFields.devStress.addArray(cells,devStressField);
            shared_ptr<TArray> VMStressField(new TArray(cells.getCountLevel1()));
            VMStressField->zero();
            _structureFields.VMStress.addArray(cells,VMStressField);
            shared_ptr<VGradArray> plasticStrainField(new VGradArray(cells.getCountLevel1()));
            plasticStrainField->zero();
            _structureFields.plasticStrain.addArray(cells,plasticStrainField);
	    shared_ptr<TArray> acCell(new TArray(cells.getCountLevel1()));
	    *acCell = _options.A;
	    _structureFields.creepConstant.addArray(cells, acCell);
	}    

        shared_ptr<TArray> rhoCell(new TArray(cells.getCountLevel1()));
        *rhoCell = vc["density"];
        _structureFields.density.addArray(cells,rhoCell);

        shared_ptr<TArray> etaCell(new TArray(cells.getCountLevel1()));
        *etaCell = vc["eta"];
        _structureFields.eta.addArray(cells,etaCell);

        shared_ptr<TArray> eta1Cell(new TArray(cells.getCountLevel1()));
        *eta1Cell = vc["eta1"];
        _structureFields.eta1.addArray(cells,eta1Cell);

        shared_ptr<TArray> alphaCell(new TArray(cells.getCountLevel1()));
        *alphaCell = vc["alpha"];
        _structureFields.alpha.addArray(cells,alphaCell);

        shared_ptr<TArray> tCell(new TArray(cells.getCountLevel1()));
        *tCell = _options["operatingTemperature"];
        _structureFields.temperature.addArray(cells,tCell);

        // compute values of deformation flux

        // store deformation flux at interfaces
	/*
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            shared_ptr<VectorT3Array> deformationFlux(new VectorT3Array(faces.getCount()));
            deformationFlux->zero();
            _structureFields.deformationFlux.addArray(faces,deformationFlux);
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
	        _structureFields.deformationFlux.addArray(faces,deformationFlux);
	    }
	}

    }
     _structureFields.eta.syncLocal();
     _structureFields.eta1.syncLocal();
     _structureFields.alpha.syncLocal();
     _structureFields.density.syncLocal();
     
     

    _niters  =0;
    _initialDeformationNorm = MFRPtr();
  }
  
  StructureBCMap& getBCMap() {return _bcMap;}
  StructureVCMap& getVCMap() {return _vcMap;}
  StructureModelOptions<T>& getOptions() {return _options;}

  void updateTime()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        VectorT3Array& w =
          dynamic_cast<VectorT3Array&>(_structureFields.deformation[cells]);
        VectorT3Array& wN1 =
          dynamic_cast<VectorT3Array&>(_structureFields.deformationN1[cells]);
        VectorT3Array& wN2 =
          dynamic_cast<VectorT3Array&>(_structureFields.deformationN2[cells]);
        if (_options.timeDiscretizationOrder > 1)
	{
            VectorT3Array& wN3 =
              dynamic_cast<VectorT3Array&>(_structureFields.deformationN3[cells]);
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

    }
  }

  void creepInit()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        VectorT3Array& w =
          dynamic_cast<VectorT3Array&>(_structureFields.deformation[cells]);
        VectorT3Array& ew =
          dynamic_cast<VectorT3Array&>(_structureFields.elasticDeformation[cells]);
        ew = w;

        _deformationGradientModel.compute();
        const VGradArray& wGradCell =
          dynamic_cast<const VGradArray&>(_structureFields.deformationGradient[cells]);

        const TArray& muCell =
          dynamic_cast<const TArray&>(_structureFields.eta[cells]);

        const TArray& lambdaCell =
          dynamic_cast<const TArray&>(_structureFields.eta1[cells]);


        VGradArray& devStress =
          dynamic_cast<VGradArray&>(_structureFields.devStress[cells]);
        TArray& VMStress =
          dynamic_cast<TArray&>(_structureFields.VMStress[cells]);

        const int nCells = cells.getCount();
        const T onethird(1.0/3.0);
        const T half(0.5);
        const T twothirds(2.0/3.0);
        const T six(6.0);
        const T zero(0.0);

        const VectorT3Array& xc =
          dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);
        T xm(1.0);
        T xv(1.0);
        T xI(1.0);
        const T one(1.0);
        const T two(2.0);
        const T three(3.0);
        const T four(4.0);
        const T eight(8.0);
        const T twelve(12.0);
        const T xl(400.e-6);
        const T mid(200.e-6);
        const T th(4.e-6);
        const T midh(2.e-6);
        xI = pow(th,three)/twelve;

        for(int k=0; k<nCells; k++)
          {
            const VGradType& wg = wGradCell[k];
            VGradType wgPlusTranspose = wGradCell[k];
            VGradType stress = wGradCell[k];
            
            for(int i=0;i<3;i++)
              for(int j=0;j<3;j++)
                wgPlusTranspose[i][j] += wg[j][i];

            stress[0][0] = wgPlusTranspose[0][0]*muCell[k]+
              (wg[0][0]+wg[1][1]+wg[2][2])*lambdaCell[k];
            stress[0][1] = wgPlusTranspose[0][1]*muCell[k];
            stress[0][2] = wgPlusTranspose[0][2]*muCell[k];
            
            stress[1][0] = wgPlusTranspose[1][0]*muCell[k];
            stress[1][1] = wgPlusTranspose[1][1]*muCell[k]+
              (wg[0][0]+wg[1][1]+wg[2][2])*lambdaCell[k];
            stress[1][2] = wgPlusTranspose[1][2]*muCell[k];

            stress[2][0] = wgPlusTranspose[2][0]*muCell[k];
            stress[2][1] = wgPlusTranspose[2][1]*muCell[k];
            stress[2][2] = wgPlusTranspose[2][2]*muCell[k]+
              (wg[0][0]+wg[1][1]+wg[2][2])*lambdaCell[k];

            if (mesh.getDimension() == 2)
	      {         
                stress[2][0] = zero;
                stress[2][1] = zero;
                stress[2][2] = zero;
	      }
            
            /*            
            if(xc[k][0]<=mid)
            {
                xm = -xl*two/eight+two*xc[k][0]/two;
                xv = one;
            }
            else
            {
                xm = (two/eight)*(three*xl-four*xc[k][0]);
                xv = -one;
            }

            stress[0][0] = -xm*(xc[k][1]-midh)/xI;
            stress[0][1] = (three/two)*(xv/th)*(pow(midh,two)-
                                                pow((xc[k][1]-midh),two))/pow(midh,two);
            stress[0][2] = zero;
            stress[1][0] = (three/two)*(xv/th)*(pow(midh,two)-
                                                pow((xc[k][1]-midh),two))/pow(midh,two);
            stress[1][1] = zero;
            stress[1][2] = zero;
            */

            devStress[k] = stress;
            const T trace = stress[0][0]+stress[1][1]+stress[2][2];
            (devStress[k])[0][0] =  (devStress[k])[0][0] - onethird*trace;
            (devStress[k])[1][1] =  (devStress[k])[1][1] - onethird*trace;
            (devStress[k])[2][2] =  (devStress[k])[2][2] - onethird*trace;
            
            VMStress[k] = sqrt(half*(pow((stress[0][0]-stress[1][1]),2.0) +
                                     pow((stress[1][1]-stress[2][2]),2.0) +
                                     pow((stress[2][2]-stress[0][0]),2.0) +
                                     six*(pow((stress[0][1]),2.0) +
                                          pow((stress[1][2]),2.0) +
                                          pow((stress[2][0]),2.0))));
          }
    
      }
  }

  void computeVMStress()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        VectorT3Array& w =
          dynamic_cast<VectorT3Array&>(_structureFields.deformation[cells]);
        _deformationGradientModel.compute();
        const VGradArray& wGradCell =
          dynamic_cast<const VGradArray&>(_structureFields.deformationGradient[cells]);

        const TArray& muCell =
          dynamic_cast<const TArray&>(_structureFields.eta[cells]);

        const TArray& lambdaCell =
          dynamic_cast<const TArray&>(_structureFields.eta1[cells]);

        VGradArray& devStress =
          dynamic_cast<VGradArray&>(_structureFields.devStress[cells]);
        const VGradArray& plasticStrainCell =
          dynamic_cast<VGradArray&>(_structureFields.plasticStrain[cells]);
        TArray& VMStress =
          dynamic_cast<TArray&>(_structureFields.VMStress[cells]);

        const int nCells = cells.getCount();
        const T onethird(1.0/3.0);
        const T half(0.5);
        const T twothirds(2.0/3.0);
        const T two(2.0);
        const T three(3.0);
        const T six(6.0);
        const T zero(0.0);

        if (mesh.getDimension() == 3)
          {
            for(int k=0; k<nCells; k++)
              {
                const VGradType& wg = wGradCell[k];
                VGradType wgPlusTranspose = wGradCell[k];
                VGradType stress = wGradCell[k];
                const VGradType pS = plasticStrainCell[k];

                for(int i=0;i<3;i++)
                  for(int j=0;j<3;j++)
                    wgPlusTranspose[i][j] += wg[j][i];

                stress[0][0] = wgPlusTranspose[0][0]*muCell[k]-
                  two*pS[0][0]*muCell[k]+
                  (wg[0][0]+wg[1][1]+wg[2][2])*lambdaCell[k];
                stress[0][1] = wgPlusTranspose[0][1]*muCell[k]-
                  two*pS[0][1]*muCell[k];
                stress[0][2] = wgPlusTranspose[0][2]*muCell[k]-
                  two*pS[0][2]*muCell[k];

                stress[1][0] = wgPlusTranspose[1][0]*muCell[k]-
                  two*pS[1][0]*muCell[k];
                stress[1][1] = wgPlusTranspose[1][1]*muCell[k]-
                  two*pS[1][1]*muCell[k]+
                  (wg[0][0]+wg[1][1]+wg[2][2])*lambdaCell[k];
                stress[1][2] = wgPlusTranspose[1][2]*muCell[k]-
                  two*pS[1][2]*muCell[k];

                stress[2][0] = wgPlusTranspose[2][0]*muCell[k]-
                  two*pS[2][0]*muCell[k];
                stress[2][1] = wgPlusTranspose[2][1]*muCell[k]-
                  two*pS[2][1]*muCell[k];
                stress[2][2] = wgPlusTranspose[2][2]*muCell[k]-
                  two*pS[2][2]*muCell[k]+
                  (wg[0][0]+wg[1][1]+wg[2][2])*lambdaCell[k];

                if(_options.creepModel!=3)
                  {
                    stress[0][0]-=(pS[0][0]+pS[1][1]+pS[2][2])*lambdaCell[k];
                    stress[1][1]-=(pS[0][0]+pS[1][1]+pS[2][2])*lambdaCell[k];
                    stress[2][2]-=(pS[0][0]+pS[1][1]+pS[2][2])*lambdaCell[k];
                  }

                devStress[k] = stress;
                const T trace = stress[0][0]+stress[1][1]+stress[2][2];
                (devStress[k])[0][0] =  (devStress[k])[0][0] - onethird*trace;
                (devStress[k])[1][1] =  (devStress[k])[1][1] - onethird*trace;
                (devStress[k])[2][2] =  (devStress[k])[2][2] - onethird*trace;

                VMStress[k] = sqrt(half*(pow((stress[0][0]-stress[1][1]),2.0) +
                                         pow((stress[1][1]-stress[2][2]),2.0) +
                                         pow((stress[2][2]-stress[0][0]),2.0) +
                                         six*(pow((stress[0][1]),2.0) +
                                              pow((stress[1][2]),2.0) +
                                              pow((stress[2][0]),2.0))));
              }
          }
        else
          {
            for(int k=0; k<nCells; k++)
              {
                const VGradType& wg = wGradCell[k];
                VGradType wgPlusTranspose = wGradCell[k];
                VGradType stress = wGradCell[k];
                const VGradType pS = plasticStrainCell[k];

                for(int i=0;i<3;i++)
                  for(int j=0;j<3;j++)
                    wgPlusTranspose[i][j] += wg[j][i];

                stress[0][0] = wgPlusTranspose[0][0]*muCell[k]-
                  two*pS[0][0]*muCell[k]+
                  (wg[0][0]+wg[1][1]+wg[2][2])*lambdaCell[k];
                stress[0][1] = wgPlusTranspose[0][1]*muCell[k]-
                  two*pS[0][1]*muCell[k];
                stress[0][2] = wgPlusTranspose[0][2]*muCell[k]-
                  two*pS[0][2]*muCell[k];

                stress[1][0] = wgPlusTranspose[1][0]*muCell[k]-
                  two*pS[1][0]*muCell[k];
                stress[1][1] = wgPlusTranspose[1][1]*muCell[k]-
                  two*pS[1][1]*muCell[k]+
                  (wg[0][0]+wg[1][1]+wg[2][2])*lambdaCell[k];
                stress[1][2] = wgPlusTranspose[1][2]*muCell[k]-
                  two*pS[1][2]*muCell[k];


                if(_options.creepModel!=3)
                  {
                    stress[0][0]+=(pS[2][2])*lambdaCell[k];
                    stress[1][1]+=(pS[2][2])*lambdaCell[k];
                  }

                stress[2][0] = zero;
                stress[2][1] = zero;
                stress[2][2] = zero;

                devStress[k] = stress;
                const T trace = stress[0][0]+stress[1][1]+stress[2][2];
                (devStress[k])[0][0] =  (devStress[k])[0][0] - onethird*trace;
                (devStress[k])[1][1] =  (devStress[k])[1][1] - onethird*trace;
                (devStress[k])[2][2] =  (devStress[k])[2][2] - onethird*trace;
                //devStress[k][2][2] = zero;
                /*
                VMStress[k] = sqrt((three/two)*(pow(devStress[k][0][0],2.0) +
                                                pow(devStress[k][1][1],2.0) +
                                                pow(devStress[k][2][2],2.0) +
                                                two*(pow((devStress[k][0][1]),2.0) +
                                                     pow((devStress[k][1][2]),2.0) +
                                                     pow((devStress[k][2][0]),2.0))));
		*/
                VMStress[k] = sqrt(half*(pow((stress[0][0]-stress[1][1]),2.0) +
                                         pow((stress[1][1]-stress[2][2]),2.0) +
                                         pow((stress[2][2]-stress[0][0]),2.0) +
                                         six*(pow((stress[0][1]),2.0) +
                                              pow((stress[1][2]),2.0) +
                                              pow((stress[2][0]),2.0))));

              }       
          }
      }
  }
  
  void initDeformationLinearization(LinearSystem& ls)
  {
  
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex wIndex(&_structureFields.deformation,&cells);

        ls.getX().addArray(wIndex,_structureFields.deformation.getArrayPtr(cells));
	
        const CRConnectivity& cellCells2 = mesh.getCellCells2();
        
        shared_ptr<Matrix> m(new CRMatrix<DiagTensorT3,DiagTensorT3,VectorT3>(cellCells2));

        ls.getMatrix().addMatrix(wIndex,wIndex,m);

    }
  }

  void applyBC(LinearSystem& ls, bool explicitMode)
  {
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
            const StructureBC<T>& bc = *_bcMap[fg.id];
            

            StructureBCS<VectorT3,DiagTensorT3,DiagTensorT3> gbc(faces,mesh,
                                                    _geomFields,
                                                    _structureFields.deformation,
                                                                 ls.getMatrix(), ls.getX(), ls.getB(),
                                                                 explicitMode);
	    
            VectorT3 fluxB(NumTypeTraits<VectorT3>::getZero());

	    if (bc.bcType == "SpecifiedDeformation")
            {
	        FloatValEvaluator<VectorT3>
		  bDeformation(bc.getVal("specifiedXDeformation"),
			       bc.getVal("specifiedYDeformation"),
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
            else if (bc.bcType == "Interface")
            {
                gbc.applyInterfaceBC();
	    }
            else if (bc.bcType == "ZeroDerivative")
            {
                gbc.applyZeroDerivativeBC();
	    }
            else if (bc.bcType == "SpecifiedTraction")
            {
	        TractionValEvaluator<T>
		  bTraction(bc.getVal("specifiedXXTraction"),
			    bc.getVal("specifiedXYTraction"),
			    bc.getVal("specifiedXZTraction"),
			    bc.getVal("specifiedYXTraction"),
			    bc.getVal("specifiedYYTraction"),
			    bc.getVal("specifiedYZTraction"),
			    bc.getVal("specifiedZXTraction"),
			    bc.getVal("specifiedZYTraction"),
			    bc.getVal("specifiedZZTraction"),
			    _geomFields,
			    faces);
	        for(int f=0; f<nFaces; f++)
                {
                    
		    fluxB += gbc.applyNeumannBC(f,bTraction[f]);
        
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
            else if (bc.bcType == "SpecifiedDistForce")
	    {
                FloatValEvaluator<VectorT3>
                  bDistForce(bc.getVal("specifiedXDistForce"),
			     bc.getVal("specifiedYDistForce"),
			     bc.getVal("specifiedZDistForce"),
			     faces);
                for(int f=0; f<nFaces; f++)
		{
                    fluxB += gbc.applyNeumannBC(f,bDistForce[f]);

		}
	    }
            else
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


    if(allNeumann && !explicitMode)
    {
        const Mesh& mesh = *_meshes[0];
	const StorageSite& cells = mesh.getCells();
	MultiField& b = ls.getB();
	MultiFieldMatrix& matrix = ls.getMatrix();
	MultiField::ArrayIndex wIndex(&_structureFields.deformation,&cells);
	VectorT3Array& rCell = dynamic_cast<VectorT3Array&>(b[wIndex]);
	VectorT3Array& w = dynamic_cast<VectorT3Array&>
	  (_structureFields.deformation[cells]);
	VVMatrix& vvMatrix =
	  dynamic_cast<VVMatrix&>(matrix.getMatrix(wIndex,wIndex));
	rCell[0] = T(0);
        w[0] = T(0);
        vvMatrix.setDirichlet(0);
    }

  }

  
  void linearizeDeformation(LinearSystem& ls, bool explicitMode)
  {
    _deformationGradientModel.compute();
    DiscrList discretizations;
    
    if(!_options.creep)
    {        
	shared_ptr<Discretization>
	  sd(new StructureSourceDiscretization<T,DiagTensorT3,DiagTensorT3>
	     (_meshes,_geomFields,
	      _structureFields.deformation,
	      _structureFields.eta,
	      _structureFields.eta1,
	      _structureFields.alpha,
	      _structureFields.deformationGradient,
	      _structureFields.temperature,
	      _options["operatingTemperature"],
	      _options["residualXXStress"],
	      _options["residualYYStress"],
	      _options["residualZZStress"],
	      _options.thermo,
	      _options.residualStress));
      
	//    shared_ptr<Discretization>
	//  bfd(new SourceDiscretization<VectorT3>
	//     (_meshes,_geomFields,
	//    _structureFields.deformation,
	//    _structureFields.bodyForce));

	discretizations.push_back(sd);
	//discretizations.push_back(bfd);
    }
   
    if (_options.creep)
    {
        shared_ptr<Discretization>
          scd(new StructurePlasticDiscretization<T,DiagTensorT3,DiagTensorT3>
              (_meshes,_geomFields,
               _structureFields.deformation,
               _structureFields.eta,
               _structureFields.eta1,
	       _structureFields.alpha,
               _structureFields.deformationGradient,
	       _structureFields.temperature,
	       _options["operatingTemperature"],
	       _options["residualXXStress"],
	       _options["residualYYStress"],
	       _options["residualZZStress"],
	       _options.thermo,
	       _options.residualStress,
               _structureFields.devStress,
               _structureFields.VMStress,
               _structureFields.plasticStrain,
	       _structureFields.creepConstant,
               _options.A,
               _options.B,
               _options.m,
               _options.n,
               _options.Sy0,
               _options["timeStep"],
               _options.creepModel));

        discretizations.push_back(scd);
    }
 
    if (_options.transient && !explicitMode)
    {
        shared_ptr<Discretization>
          td(new TimeDerivativeStructureDiscretization
	     <VectorT3,DiagTensorT3,DiagTensorT3>
             (_meshes,_geomFields,
              _structureFields.deformation,
              _structureFields.deformationN1,
              _structureFields.deformationN2,
              _structureFields.deformationN3,
              _structureFields.density,
	      _structureFields.volume0,
	      _options.variableTimeStep,
              _options["timeStep"],
              _options.timeStepN1,
              _options.timeStepN2));
        
        discretizations.push_back(td);
    }
 
    /*
    shared_ptr<Discretization>
      ibm(new GenericIBDiscretization<VectorT3,DiagTensorT3,T>
	  (_meshes,_geomFields,_structureFields.deformation));
      
    discretizations.push_back(ibm);
    */
    Linearizer linearizer;

    linearizer.linearize(discretizations,_meshes,ls.getMatrix(),
                         ls.getX(), ls.getB());

    if (!explicitMode)
      applyBC(ls,explicitMode);

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
    linearizeDeformation(ls,false);
    ls.initSolve();    
    //AMG solver(ls);
    MFRPtr rNorm = _options.getDeformationLinearSolver().solve(ls);
    if (!_initialDeformationNorm) _initialDeformationNorm = rNorm;
    _options.getDeformationLinearSolver().cleanup();
    ls.postSolve();
    ls.updateSolution();
    //postStructureSolve(ls);
    return rNorm;
  }

  void initExplicitAdvance()
  {
    this->_lsK = shared_ptr<LinearSystem>(new LinearSystem());
    LinearSystem& ls = *(this->_lsK);
    initDeformationLinearization(ls);
    ls.initAssembly();
    linearizeDeformation(ls,true);
    // add delta explicitly since we won't be calling initSolve
    ls.replaceDelta(dynamic_pointer_cast<MultiField>(ls.getX().newClone()));
    ls.replaceResidual(dynamic_pointer_cast<MultiField>(ls.getX().newClone()));
    
  }
  
  void finishExplicitAdvance()
  {
    this->_lsK = shared_ptr<LinearSystem>();
  }

  void advanceExplicit(const int nSteps, const double deltaT)
  {
    const int nMeshes = this->_meshes.size();
    LinearSystem& ls = *(this->_lsK);
    
    TimeDerivativeStructureDiscretization
      <VectorT3,DiagTensorT3,DiagTensorT3>
      td(_meshes,_geomFields,
         _structureFields.deformation,
         _structureFields.deformationN1,
         _structureFields.deformationN2,
         _structureFields.deformationN3,
         _structureFields.density,
         _structureFields.volume0,
	 _options.variableTimeStep,
         deltaT,
	 _options.timeStepN1,
	 _options.timeStepN2);
    
        
    for(int n=0; n<nSteps; n++)
    {

        ls.getMatrix().multiply(ls.getB(), ls.getX());

        
        applyBC(ls, true);

        ls.getDelta().zero();

        for(int n=0; n<nMeshes; n++)
        {
            const Mesh& mesh = *(this->_meshes[n]);
            td.explicitAdvance(mesh,ls.getX(),ls.getB(),ls.getDelta());
        }

        // this will update the boundary delta's
        ls.postSolve();

        ls.updateSolution();

        updateTime();
    }
  }
  
#if 0

  void postStructureSolve(LinearSystem& ls)
  {
    MultiField& sField = ls.getDelta();

    const int numMeshes = _meshes.size();
    for(int n=0;n<numMeshes;n++)
    {
        const Mesh& mesh = *_meshes[n];

	const StorageSite& cells = mesh.getCells();

	MultiField::ArrayIndex sIndex(&_structureFields.deformation,&cells);
	VectorT3Array& w = dynamic_cast<VectorT3Array&>
	  (_structureFields.deformation[cells]);
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

#endif

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
        

#ifdef FVM_PARALLEL	
        if ( MPI::COMM_WORLD.Get_rank() == 0 ){  //only root process
          if (_options.printNormalizedResiduals)
             cout << _niters << ": " << *dNormRatio <<  endl;
          else
             cout << _niters << ": " << *dNorm <<  endl;
        }	     
#endif

#ifndef FVM_PARALLEL	
        if (_options.printNormalizedResiduals)
          cout << _niters << ": " << *dNormRatio <<  endl;
        else
          cout << _niters << ": " << *dNorm <<  endl;
#endif


        _niters++;
        if (*dNormRatio < _options.deformationTolerance)
          return true;
    }
    return false;
  }

#if 0
  void printBCs()
  {
    foreach(typename StructureBCMap::value_type& pos, _bcMap)
    {
        cout << "Face Group " << pos.first << ":" << endl;
        cout << "    bc type " << pos.second->bcType << endl;
        foreach(typename StructureBC<T>::value_type& vp, *pos.second)
        {
            cout << "   " << vp.first << " "  << vp.second.constant <<  endl;
        }
    }
  }


  VectorT3 getDeformationFluxIntegral(const Mesh& mesh, const int faceGroupId)
  {
    VectorT3 r(VectorT3::getZero());
    bool found = false;
    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
        const FaceGroup& fg = *fgPtr;
        if (fg.id == faceGroupId)
        {
            const StorageSite& faces = fg.site;
            const int nFaces = faces.getCount();
            const VectorT3Array& deformationFlux =
              dynamic_cast<const VectorT3Array&>(_structureFields.deformationFlux[faces]);
            for(int f=0; f<nFaces; f++)
              r += deformationFlux[f];
            found=true;
        }
    }
    if (!found)
      throw CException("getDeformationFluxIntegral: invalid faceGroupID");
    return r;
  }

  VectorT3 getDeformationDerivativeIntegral(const Mesh& mesh)
  {
    VectorT3 r(VectorT3::getZero());
    const StorageSite& cells = mesh.getCells();
    const int nCells = cells.getSelfCount();

    const TArray& density =
          dynamic_cast<const TArray&>(_structureFields.density[cells]);
    const VectorT3Array& w =
          dynamic_cast<const VectorT3Array&>(_structureFields.deformation[cells]);
    const VectorT3Array& wN1 =
      dynamic_cast<const VectorT3Array&>(_structureFields.deformationN1[cells]);
    const VectorT3Array& wN2 =
      dynamic_cast<const VectorT3Array&>(_structureFields.deformationN2[cells]);

    const TArray& volume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    
    const T deltaT = _options["timeStep"];
    
    T onePointFive(1.5);
    T two(2.0);
    T pointFive(0.5);

    for(int c=0; c<nCells; c++)
    {
	const T rhoVbydT = density[c]*volume[c]/deltaT;
	r += rhoVbydT*(onePointFive*w[c]- two*wN1[c]
		       + pointFive*wN2[c]);
    }
    return r;
  }

  boost::shared_ptr<ArrayBase> getStressTensor(const Mesh& mesh, const ArrayBase& gcellIds)
  {
    typedef Array<StressTensor<T> > StressTensorArray;
    
    const StorageSite& cells = mesh.getCells();
    
    const Array<int>& cellIds = dynamic_cast<const Array<int> &>(gcellIds);
    const int nCells = cellIds.getLength();

    _deformationGradientModel.compute();

    const VGradArray& wGrad =
      dynamic_cast<const VGradArray&>(_structureFields.deformationGradient[cells]);

    /*    const TArray& pCell =
      dynamic_cast<const TArray&>(_flowFields.pressure[cells]);
    */
    const TArray& eta = dynamic_cast<const TArray&>(_structureFields.eta[cells]);
    const TArray& eta1 = dynamic_cast<const TArray&>(_structureFields.eta1[cells]);

    boost::shared_ptr<StressTensorArray> stressTensorPtr( new StressTensorArray(nCells));
    StressTensorArray& stressTensor = *stressTensorPtr;

    for(int n=0; n<nCells; n++)
    {
        const int c = cellIds[n];
        const VGradType& wg = wGrad[c];
        VGradType wgPlusTranspose = wGrad[c];

        for(int i=0;i<3;i++)
          for(int j=0;j<3;j++)
            wgPlusTranspose[i][j] += wg[j][i];
        
        stressTensor[n][0] = wgPlusTranspose[0][0]*eta[c]+
	  (wg[0][0]+wg[1][1]+wg[2][2])*eta1[c];
        stressTensor[n][1] = wgPlusTranspose[1][1]*eta[c]+
	  (wg[0][0]+wg[1][1]+wg[2][2])*eta1[c]; 
        stressTensor[n][2] = wgPlusTranspose[2][2]*eta[c]+
	  (wg[0][0]+wg[1][1]+wg[2][2])*eta1[c];
        stressTensor[n][3] = wgPlusTranspose[0][1]*eta[c];
        stressTensor[n][4] = wgPlusTranspose[1][2]*eta[c];
        stressTensor[n][5] = wgPlusTranspose[2][0]*eta[c];
    }

    return stressTensorPtr;
  }

#endif
    
  void getTraction(const Mesh& mesh)
  {
      const StorageSite& cells = mesh.getCells();

      const int nCells = cells.getCount();

      shared_ptr<VectorT3Array> tractionXPtr(new VectorT3Array(nCells));
      tractionXPtr->zero();
      _structureFields.tractionX.addArray(cells,tractionXPtr);
      VectorT3Array& tractionX = *tractionXPtr;

      shared_ptr<VectorT3Array> tractionYPtr(new VectorT3Array(nCells));
      tractionYPtr->zero();
      _structureFields.tractionY.addArray(cells,tractionYPtr);
      VectorT3Array& tractionY = *tractionYPtr;

      shared_ptr<VectorT3Array> tractionZPtr(new VectorT3Array(nCells));
      tractionZPtr->zero();
      _structureFields.tractionZ.addArray(cells,tractionZPtr);
      VectorT3Array& tractionZ = *tractionZPtr;

      _deformationGradientModel.compute();

      const VGradArray& wGrad =
        dynamic_cast<const VGradArray&>(_structureFields.deformationGradient[cells]);

      const VGradArray& plasticStrainCell =
        dynamic_cast<VGradArray&>(_structureFields.plasticStrain[cells]);

      const TArray& eta = dynamic_cast<const TArray&>(_structureFields.eta[cells]);
      const TArray& eta1 = dynamic_cast<const TArray&>(_structureFields.eta1[cells]);
      const TArray& alpha = dynamic_cast<const TArray&>(_structureFields.alpha[cells]);

      const TArray& temperature = dynamic_cast<const TArray&>(_structureFields.temperature[cells]);
      
      const T two(2.0);
      const T three(3.0);
      const T zero(0.0);
      
      for(int n=0; n<nCells; n++)
      {
	  const VGradType& wg = wGrad[n];
	  VGradType wgPlusTranspose = wGrad[n];
	  const VGradType& pS =  plasticStrainCell[n];
	  for(int i=0;i<3;i++)
	    for(int j=0;j<3;j++)
	      wgPlusTranspose[i][j] += wg[j][i];
	  
	  tractionX[n][0] = wgPlusTranspose[0][0]*eta[n]+
	 	    (wg[0][0]+wg[1][1]+wg[2][2])*eta1[n];
	  tractionX[n][1] = wgPlusTranspose[0][1]*eta[n];
	  tractionX[n][2] = wgPlusTranspose[0][2]*eta[n];

          tractionY[n][0] = wgPlusTranspose[1][0]*eta[n];
          tractionY[n][1] = wgPlusTranspose[1][1]*eta[n]+
	    (wg[0][0]+wg[1][1]+wg[2][2])*eta1[n];
          tractionY[n][2] = wgPlusTranspose[1][2]*eta[n];

          tractionZ[n][0] = wgPlusTranspose[2][0]*eta[n];
          tractionZ[n][1] = wgPlusTranspose[2][1]*eta[n];
          tractionZ[n][2] = wgPlusTranspose[2][2]*eta[n]+
	    (wg[0][0]+wg[1][1]+wg[2][2])*eta1[n];	  

	  if(_options.residualStress)
	  {
	      tractionX[n][0] += _options["residualXXStress"];
	      tractionY[n][1] += _options["residualYYStress"];
	      tractionZ[n][2] += _options["residualZZStress"];
	  }

	  if(_options.thermo)
	  {
	      if(mesh.getDimension()==2)
	      {
		  tractionX[n][0] -= (two*eta1[n]+two*eta[n])*alpha[n]*(temperature[n]-_options["operatingTemperature"]);
		  tractionY[n][1] -= (two*eta1[n]+two*eta[n])*alpha[n]*(temperature[n]-_options["operatingTemperature"]);
	      }
	      else
	      {
		  tractionX[n][0] -= (three*eta1[n]+two*eta[n])*alpha[n]*(temperature[n]-_options["operatingTemperature"]);
		  tractionY[n][1] -= (three*eta1[n]+two*eta[n])*alpha[n]*(temperature[n]-_options["operatingTemperature"]);
		  tractionZ[n][2] -= (three*eta1[n]+two*eta[n])*alpha[n]*(temperature[n]-_options["operatingTemperature"]);
	      }
	  }

          if(_options.creep)
          {
              tractionX[n][0]-=(two*eta[n]*pS[0][0] + 
                                (pS[0][0] + pS[1][1])*eta1[n]);
              tractionX[n][1]-=(two*eta[n]*pS[0][1]);
              tractionX[n][2]-=(two*eta[n]*pS[0][2]);

              tractionY[n][0]-=(two*eta[n]*pS[1][0]);
              tractionY[n][1]-=(two*eta[n]*pS[1][1] +
                                (pS[0][0] + pS[1][1])*eta1[n]);
              tractionY[n][2]-=(two*eta[n]*pS[1][2]);
	  }

          if (mesh.getDimension() == 2)
          {
              tractionZ[n] = zero;
	  }
      }
  }

  void getStrain(const Mesh& mesh)
  {
    const StorageSite& cells = mesh.getCells();

    const int nCells = cells.getCount();

    shared_ptr<VectorT3Array> strainXPtr(new VectorT3Array(nCells));
    strainXPtr->zero();
    _structureFields.strainX.addArray(cells,strainXPtr);
    VectorT3Array& strainX = *strainXPtr;

    shared_ptr<VectorT3Array> strainYPtr(new VectorT3Array(nCells));
    strainYPtr->zero();
    _structureFields.strainY.addArray(cells,strainYPtr);
    VectorT3Array& strainY = *strainYPtr;

    shared_ptr<VectorT3Array> strainZPtr(new VectorT3Array(nCells));
    strainZPtr->zero();
    _structureFields.strainZ.addArray(cells,strainZPtr);
    VectorT3Array& strainZ = *strainZPtr;

    _deformationGradientModel.compute();

    const VGradArray& wGrad =
      dynamic_cast<const VGradArray&>(_structureFields.deformationGradient[cells]);
    
    const T half(1./2.);
    
    for(int n=0; n<nCells; n++)
      {
        const VGradType& wg = wGrad[n];
        VGradType wgPlusTranspose = wGrad[n];

        for(int i=0;i<3;i++)
          for(int j=0;j<3;j++)
            wgPlusTranspose[i][j] += wg[j][i];

        strainX[n][0] = half*wgPlusTranspose[0][0];
        strainX[n][1] = half*wgPlusTranspose[0][1];
        strainX[n][2] = half*wgPlusTranspose[0][2];

        strainY[n][0] = half*wgPlusTranspose[1][0];
        strainY[n][1] = half*wgPlusTranspose[1][1];
        strainY[n][2] = half*wgPlusTranspose[1][2];

        strainZ[n][0] = half*wgPlusTranspose[2][0];
        strainZ[n][1] = half*wgPlusTranspose[2][1];
        strainZ[n][2] = half*wgPlusTranspose[2][2];
      }
  }

  void getPlasticDiagStrain(const Mesh& mesh)
  {
    const StorageSite& cells = mesh.getCells();

    const int nCells = cells.getSelfCount();

    shared_ptr<VectorT3Array> plasticDiagStrainPtr(new VectorT3Array(nCells));
    plasticDiagStrainPtr->zero();
    _structureFields.plasticDiagStrain.addArray(cells,plasticDiagStrainPtr);
    VectorT3Array& plasticDiagStrain = *plasticDiagStrainPtr;

    const VGradArray& plasticStrain =
      dynamic_cast<const VGradArray&>(_structureFields.plasticStrain[cells]);

    for(int n=0; n<nCells; n++)
      {
        plasticDiagStrain[n][0]=plasticStrain[n][0][0];
        plasticDiagStrain[n][1]=plasticStrain[n][1][1];
        plasticDiagStrain[n][2]=plasticStrain[n][0][1];
      }
  }

void updateForceOnBoundary(const StorageSite& faceSite, const ArrayBase& bforceA, const map<int,int>& commonFacesMap, 
                           ArrayBase& fxA, ArrayBase& fyA, ArrayBase& fzA)
{
    //bforce came from fluid+elec side
    const VectorT3Array& bforce = dynamic_cast<const VectorT3Array&>(bforceA);
    //following will be updated
    TArray& fx = dynamic_cast<TArray&> (fxA);
    TArray& fy = dynamic_cast<TArray&> (fyA);
    TArray& fz = dynamic_cast<TArray&> (fzA);
    const int offset = faceSite.getOffset();
    for (int i = 0; i < faceSite.getCount(); i++ ){
        const int faceID = i + offset; //localface ID
        //commonFacesMap will get right index in bforce
        const int indx = commonFacesMap.find(faceID)->second;
        fx[i] = bforce[indx][0];
        fy[i] = bforce[indx][1];
        fz[i] = bforce[indx][2];
    }
        
}

  /*
  VectorT3 getMomentumFluxIntegralonIBFaces(const Mesh& mesh)
  {
    VectorT3 r(VectorT3::getZero());
    const StorageSite& ibFaces = mesh.getIBFaces();
    const StorageSite& faces = mesh.getFaces();
    const Array<int>& ibFaceIndices = mesh.getIBFaceList();
    const int nibf = ibFaces.getCount();
    const VectorT3Array& deformationFlux =
              dynamic_cast<const VectorT3Array&>(_structureFields.deformationFlux[faces]);
    for ( int f = 0; f < nibf; f ++){
      const int ibFaceIndex = ibFaceIndices[f];
      r += deformationFlux[ibFaceIndex];
    }
    if (nibf == 0)
      throw CException("getDeformationFluxIntegralonIBFaces:  no ibFaces found!");
    return r;
  }
  

  void printPressureIntegrals()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            
            VectorT3 r(VectorT3::getZero());
            
            const StorageSite& faces = fg.site;
            const VectorT3Array& faceArea =
              dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
            const int nFaces = faces.getCount();
            const TArray& facePressure = dynamic_cast<const TArray&>(_flowFields.pressure[faces]);
            for(int f=0; f<nFaces; f++)
              r += faceArea[f]*facePressure[f];

            cout << "Mesh " << mesh.getID() << " faceGroup " << fg.id << " : " << r <<  endl;
        }
    }
  }
  
  */

#if 0

  void printDeformationFluxIntegrals()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            
            VectorT3 r(VectorT3::getZero());
            
            const StorageSite& faces = fg.site;
            const int nFaces = faces.getCount();
            const VectorT3Array& deformationFlux =
              dynamic_cast<const VectorT3Array&>(_structureFields.deformationFlux[faces]);
            for(int f=0; f<nFaces; f++)
              r += deformationFlux[f];

            cout << "Mesh " << mesh.getID() << " faceGroup " << fg.id << " : " << r <<  endl;
        }
    }
  }

#endif

  /*  
  void printMassFluxIntegrals()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            
            T r(0.);
            
            const StorageSite& faces = fg.site;
            const int nFaces = faces.getCount();
            const TArray& massFlux = dynamic_cast<const TArray&>(_flowFields.massFlux[faces]);
            for(int f=0; f<nFaces; f++)
              r += massFlux[f];

            cout << "Mesh " << mesh.getID() << " faceGroup " << fg.id << " : " << r <<  endl;
        }
    }
  }

#include "FlowModelPressureBC.h"
#include "FlowModelVelocityBC.h"
#include "FlowModelInterior.h"
#include "FlowModelSlipJump.h"  

  */


private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  StructureFields& _structureFields;

  StructureBCMap _bcMap;
  StructureVCMap _vcMap;
  
  StructureModelOptions<T> _options;
  GradientModel<VectorT3> _deformationGradientModel;
  //  GradientModel<T> _pressureGradientModel;
  
  MFRPtr _initialDeformationNorm;
//  MFRPtr _initialCoupledNorm;
  int _niters;

  shared_ptr<LinearSystem> _lsK;
  //  shared_ptr<Field> _previousVelocity;
  //  shared_ptr<Field> _momApField;

  //  bool _useReferencePressure;
  //  T _referencePP;
  //AMG _momSolver;
  //AMG _continuitySolver;
};

template<class T>
StructureModel<T>::StructureModel(const GeomFields& geomFields,
                        StructureFields& structureFields,
                        const MeshList& meshes) :
  Model(meshes),
  _impl(new Impl(geomFields,structureFields,meshes))
{
  logCtor();
}


template<class T>
StructureModel<T>::~StructureModel()
{
  logDtor();
}

template<class T>
void
StructureModel<T>::init()
{
  _impl->init();
}
  
template<class T>
typename StructureModel<T>::StructureBCMap&
StructureModel<T>::getBCMap() {return _impl->getBCMap();}

template<class T>
typename StructureModel<T>::StructureVCMap&
StructureModel<T>::getVCMap() {return _impl->getVCMap();}

template<class T>
StructureModelOptions<T>&
StructureModel<T>::getOptions() {return _impl->getOptions();}

#if 0

template<class T>
void
StructureModel<T>::printBCs()
{
  _impl->printBCs();
}

#endif

template<class T>
bool
StructureModel<T>::advance(const int niter)
{
  return _impl->advance(niter);
}

template<class T>
void
StructureModel<T>::advanceExplicit(const int nsteps, const double deltaT)
{
  _impl->advanceExplicit(nsteps,deltaT);
}

template<class T>
void
StructureModel<T>::initExplicitAdvance()
{
  return _impl->initExplicitAdvance();
}


template<class T>
void
StructureModel<T>::updateForceOnBoundary(const StorageSite& faceSite, const ArrayBase& bforceA, const map<int,int>& commonFacesMap, 
                           ArrayBase& fxA, ArrayBase& fyA, ArrayBase& fzA)

{
  _impl->updateForceOnBoundary(faceSite, bforceA, commonFacesMap, fxA, fyA, fzA);
;
}

template<class T>
void
StructureModel<T>::finishExplicitAdvance()
{
  return _impl->finishExplicitAdvance();
}

template<class T>
void
StructureModel<T>::creepInit()
{
  _impl->creepInit();
}

template<class T>
void
StructureModel<T>::computeVMStress()
{
  _impl->computeVMStress();
}

template<class T>
void
StructureModel<T>::getStrain(const Mesh& mesh)
{
  return  _impl->getStrain(mesh);
}

template<class T>
void
StructureModel<T>::getPlasticDiagStrain(const Mesh& mesh)
{
  return  _impl->getPlasticDiagStrain(mesh);
}

template<class T>
void
StructureModel<T>::updateTime()
{
  _impl->updateTime();
}

#if 0

template<class T>
void
StructureModel<T>::printDeformationFluxIntegrals()
{
  _impl->printDeformationFluxIntegrals();
}


template<class T>
Vector<T,3>
StructureModel<T>::getDeformationFluxIntegral(const Mesh& mesh, const int faceGroupId)
{
 return  _impl->getDeformationFluxIntegral(mesh,faceGroupId);
}

template<class T>
Vector<T,3>
StructureModel<T>::getDeformationDerivativeIntegral(const Mesh& mesh)
{
 return  _impl->getDeformationDerivativeIntegral(mesh);
}

#endif

template<class T>
void
StructureModel<T>::getTraction(const Mesh& mesh)
{
  return  _impl->getTraction(mesh);
}

#if 0

template<class T>
boost::shared_ptr<ArrayBase>
StructureModel<T>::getStressTensor(const Mesh& mesh, const ArrayBase& cellIds)
{
  return _impl->getStressTensor(mesh, cellIds);
}

#endif
