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
    _faceArea(dynamic_cast<const VectorT3Array&>(_areaField[_faces]))
  {}
  
  X applyDirichletBC(int f, const X& bValue) const
  {
    const int c1 = _faceCells(f,1);

    const X fluxB = -_r[c1];
    
    const X dXC1 = bValue - _x[c1];

    _dRdX.eliminateDirichlet(c1,_r,dXC1);
    _x[c1] = bValue;
    _r[c1] = NumTypeTraits<X>::getZero();
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
    const int c0 = _faceCells(f,0);
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
};

template<class T>
class PlateModel<T>::Impl
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

        shared_ptr<VectorT3Array> sCell(new VectorT3Array(cells.getCount()));

        VectorT3 initialDeformation;
        initialDeformation[0] = _options["initialXRotation"];
        initialDeformation[1] = _options["initialYRotation"];
        initialDeformation[2] = _options["initialZDeformation"];
        *sCell = initialDeformation;

        _plateFields.deformation.addArray(cells,sCell);

        if (_options.transient)
        {
	    _plateFields.volume0.addArray(cells,
					 dynamic_pointer_cast<ArrayBase>
					 (_geomFields.volume[cells].newCopy()));
	    _plateFields.deformationN1.addArray(cells,
                                            dynamic_pointer_cast<ArrayBase>(sCell->newCopy()));
            _plateFields.deformationN2.addArray(cells,
					    dynamic_pointer_cast<ArrayBase>(sCell->newCopy()));            
            if (_options.timeDiscretizationOrder > 1)
	      _plateFields.deformationN3.addArray(cells,
						      dynamic_pointer_cast<ArrayBase>(sCell->newCopy()));
        }
        

        shared_ptr<TArray> rhoCell(new TArray(cells.getCount()));
        *rhoCell = vc["density"];
        _plateFields.density.addArray(cells,rhoCell);

        shared_ptr<TArray> ymCell(new TArray(cells.getCount()));
        *ymCell = vc["ym"];
        _plateFields.ym.addArray(cells,ymCell);

        shared_ptr<TArray> nuCell(new TArray(cells.getCount()));
        *nuCell = vc["nu"];
        _plateFields.nu.addArray(cells,nuCell);

        shared_ptr<TArray> forceCell(new TArray(cells.getCount()));
        forceCell->zero();
        _plateFields.force.addArray(cells,forceCell);

        shared_ptr<TArray> thicknessCell(new TArray(cells.getCount()));
        thicknessCell->zero();
        _plateFields.thickness.addArray(cells,thicknessCell);

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
            
    shared_ptr<Discretization>
      sd(new PlateSourceDiscretization<T,DiagTensorT3,DiagTensorT3>
         (_meshes,_geomFields,
          _plateFields.deformation,
          _plateFields.ym,
	  _plateFields.nu,
          _plateFields.deformationGradient,
	  _plateFields.thickness,
	  _plateFields.force,
	  _options.scf));
      
    //    shared_ptr<Discretization>
    //  bfd(new SourceDiscretization<VectorT3>
    //     (_meshes,_geomFields,
    //    _structureFields.deformation,
    //    _structureFields.bodyForce));

    discretizations.push_back(sd);
    //discretizations.push_back(bfd);
    /*    
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
	      _plateFields.volume0,
              _options["timeStep"]));
        
        discretizations.push_back(td);
    }
    */ 
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
                
    if(allNeumann)
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
    return rNorm;
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
	const T deformationURF(_options["deformationURF"]);

	const int nCells = cells.getCount();
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
  PlateModel<T>::dumpMatrix(const string fileBase)
  {
    _impl->dumpMatrix(fileBase);
  }
