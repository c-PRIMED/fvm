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
             MultiField& xField, MultiField& rField) :
    _faces(faces),
    _cells(mesh.getCells()),
    _ibType(mesh.getIBType()),
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
  
  void applyDirichletBC(int f, const X& bValue) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);

    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    
    const X dXC1 = bValue - _x[c1];
    _dRdX.eliminateDirichlet(c1,_r,dXC1);
    _x[c1] = bValue;
    _r[c1] = NumTypeTraits<X>::getZero();
    _dRdX.setDirichlet(c1);
  }

  void applyDirichletBC(const X& bValue) const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyDirichletBC(i,bValue);
  }
  
  void applyDirichletBC(const FloatValEvaluator<X>& bValue) const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyDirichletBC(i,bValue[i]);
  }
  
  void applyNeumannBC(const int f,
                      const X& specifiedFlux) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);

    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    // the current value of flux and its Jacobians
    const X fluxB = -_r[c1];
    const X dFlux = specifiedFlux*_faceAreaMag[f] - fluxB;

    // setup the equation for the boundary value; the coefficients
    // are already computed so just need to set the rhs
    _r[c1] = dFlux;

    // _dRdX.eliminate(c1,_r);
    // mark this row as a "boundary" row so that we will update it
    // after the overall system is solved
    _dRdX.setBoundary(c1);
  }


  void applyNeumannBC(const X& bFlux) const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyNeumannBC(i,bFlux);
  }
  
  void applyNeumannBC(const FloatValEvaluator<X>& bFlux) const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyNeumannBC(i,bFlux[i]);
  }

  void applyInterfaceBC(const int f) const
  {
    throw;
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
  const Array<int>& _ibType;
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
                    bc->bcType = "SpecifiedTraction";
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

        shared_ptr<VectorT3Array> sCell(new VectorT3Array(cells.getCount()));

        VectorT3 initialDeformation;
        initialDeformation[0] = _options["initialXDeformation"];
        initialDeformation[1] = _options["initialYDeformation"];
        initialDeformation[2] = _options["initialZDeformation"];
        *sCell = initialDeformation;

        _structureFields.deformation.addArray(cells,sCell);

        if (_options.transient)
        {
	    _structureFields.deformationN1.addArray(cells,
                                            dynamic_pointer_cast<ArrayBase>(sCell->newCopy()));
            _structureFields.deformationN2.addArray(cells,
					    dynamic_pointer_cast<ArrayBase>(sCell->newCopy()));            
        }
        

        shared_ptr<TArray> rhoCell(new TArray(cells.getCount()));
        *rhoCell = vc["density"];
        _structureFields.density.addArray(cells,rhoCell);

        shared_ptr<TArray> etaCell(new TArray(cells.getCount()));
        *etaCell = vc["eta"];
        _structureFields.eta.addArray(cells,etaCell);

        shared_ptr<TArray> eta1Cell(new TArray(cells.getCount()));
        *eta1Cell = vc["eta1"];
        _structureFields.eta1.addArray(cells,eta1Cell);

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
        MultiField::ArrayIndex wIndex(&_structureFields.deformation,&cells);

        ls.getX().addArray(wIndex,_structureFields.deformation.getArrayPtr(cells));

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
      sd(new StructureSourceDiscretization<T,DiagTensorT3,DiagTensorT3>
         (_meshes,_geomFields,
          _structureFields.deformation,
          _structureFields.eta,
	  _structureFields.eta1,
          _structureFields.deformationGradient));
      
    //    discretizations.push_back(dd);
    discretizations.push_back(sd);
    
    if (_options.transient)
    {
        shared_ptr<Discretization>
          td(new TimeDerivativeStructureDiscretization
	     <VectorT3,DiagTensorT3,DiagTensorT3>
             (_meshes,_geomFields,
              _structureFields.deformation,
              _structureFields.deformationN1,
              _structureFields.deformationN2,
              _structureFields.density,
              _options["timeStep"]));
        
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
                                                    ls.getMatrix(), ls.getX(), ls.getB());
	    
	    if (bc.bcType == "SpecifiedDeformation")
            {
	        FloatValEvaluator<VectorT3>
		  bDeformation(bc.getVal("specifiedXDeformation"),
			       bc.getVal("specifiedYDeformation"),
			       bc.getVal("specifiedZDeformation"),
			       faces);
                gbc.applyDirichletBC(bDeformation);
		allNeumann = false;
            }
            else if (bc.bcType == "Symmetry")
            {
                gbc.applySymmetryBC();
                allNeumann = false;
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
                    
		    gbc.applyNeumannBC(f,bTraction[f]);
        
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

                    gbc.applyNeumannBC(f,bForce[f]/faceAreaMag[f]);

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
		    
                    gbc.applyNeumannBC(f,bDistForce[f]);

		}
	    }
            else
              throw CException(bc.bcType + " not implemented for StructureModel");
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
    //postStructureSolve(ls);
    return rNorm;
  }

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
        
        MFRPtr dNormRatio((*dNorm)/(*_initialDeformationNorm));
        
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

  /*
  VectorT3 getPressureIntegral(const Mesh& mesh, const int faceGroupId)
  {
    VectorT3 r(VectorT3::getZero());
    bool found = false;
    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
        const FaceGroup& fg = *fgPtr;
        if (fg.id == faceGroupId)
        {
            const StorageSite& faces = fg.site;
            const VectorT3Array& faceArea =
              dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
            const int nFaces = faces.getCount();
            const TArray& facePressure = dynamic_cast<const TArray&>(_flowFields.pressure[faces]);
            for(int f=0; f<nFaces; f++)
              r += faceArea[f]*facePressure[f];

            found=true;
        }
    }
    if (!found)
      throw CException("getPressureIntegral: invalid faceGroupID");
    return r;
  }

  VectorT3 getPressureIntegralonIBFaces(const Mesh& mesh)
  {
    VectorT3 r(VectorT3::getZero());
    const StorageSite& ibFaces = mesh.getIBFaces();
    const StorageSite& faces = mesh.getFaces();
    const Array<int>& ibFaceIndices = mesh.getIBFaceList();
    const int nibf = ibFaces.getCount();
    const CRConnectivity& faceCells = mesh.getFaceCells(faces);
    const Array<int>& ibType = mesh.getIBType();
    const VectorT3Array& faceArea =
              dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
    const TArray& facePressure = dynamic_cast<const TArray&>(_flowFields.pressure[faces]);

    for ( int f = 0; f < nibf; f ++)
    {
        const int ibFaceIndex = ibFaceIndices[f];
        const int c0 = faceCells(ibFaceIndex,0);

        // need to check whether the face points in or out of the
        // fluid cell and adjust sign of area accordingly
        
        if (ibType[c0] == Mesh::IBTYPE_FLUID)
          r += faceArea[ibFaceIndex]*facePressure[ibFaceIndex];
        else
          r -= faceArea[ibFaceIndex]*facePressure[ibFaceIndex];
    }
    
    if (nibf == 0)
      throw CException("getPressureIntegralonIBFaces: no IBFaces found!");
    return r;
  }

  VectorT3 getPVIntegral(const Field& velCoeffField,
                         const Mesh& mesh, const int faceGroupId)
  {
    VectorT3 r(VectorT3::getZero());
    bool found = false;
    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
        const FaceGroup& fg = *fgPtr;
        if (fg.id == faceGroupId)
        {
            const StorageSite& faces = fg.site;
            const VectorT3Array& faceArea =
              dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
            const VectorT3Array& velCoeff =
              dynamic_cast<const VectorT3Array&>(velCoeffField[faces]);
            const int nFaces = faces.getCount();
            const TArray& facePressure = dynamic_cast<const TArray&>(_flowFields.pressure[faces]);
            for(int f=0; f<nFaces; f++)
              r += velCoeff[f]*(faceArea[f]*facePressure[f]);

            found=true;
        }
    }
  if (!found)
    throw CException("getPVIntegral: invalid faceGroupID");
    return r;
  }  
  */

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
    
  void getTractionX(const Mesh& mesh)
  {
      const StorageSite& cells = mesh.getCells();

      const int nCells = cells.getSelfCount();

      shared_ptr<VectorT3Array> tractionXPtr(new VectorT3Array(nCells));
      tractionXPtr->zero();
      _structureFields.tractionX.addArray(cells,tractionXPtr);
      VectorT3Array& tractionX = *tractionXPtr;

      _deformationGradientModel.compute();

      const VGradArray& wGrad =
        dynamic_cast<const VGradArray&>(_structureFields.deformationGradient[cells]);

      const TArray& eta = dynamic_cast<const TArray&>(_structureFields.eta[cells]);
      const TArray& eta1 = dynamic_cast<const TArray&>(_structureFields.eta1[cells]);
      
      for(int n=0; n<nCells; n++)
      {
	  const VGradType& wg = wGrad[n];
	  VGradType wgPlusTranspose = wGrad[n];
	  	  	    
	  for(int i=0;i<3;i++)
	    for(int j=0;j<3;j++)
	      wgPlusTranspose[i][j] += wg[j][i];
	  
	  tractionX[n][0] = wgPlusTranspose[0][0]*eta[n]+
	 	    (wg[0][0]+wg[1][1]+wg[2][2])*eta1[n];
	  tractionX[n][1] = wgPlusTranspose[1][0]*eta[n];
	  tractionX[n][2] = wgPlusTranspose[1][1]*eta[n]+
	  	    (wg[0][0]+wg[1][1]+wg[2][2])*eta1[n];
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


template<class T>
void
StructureModel<T>::printBCs()
{
  _impl->printBCs();
}

template<class T>
bool
StructureModel<T>::advance(const int niter)
{
  return _impl->advance(niter);
}

/*
template<class T>
bool
StructureModel<T>::advanceCoupled(const int niter)
{
  return _impl->advanceCoupled(niter);
}
*/


template<class T>
void
StructureModel<T>::updateTime()
{
  _impl->updateTime();
}

/*template<class T>
void
FlowModel<T>::printPressureIntegrals()
{
  _impl->printPressureIntegrals();
}
*/


template<class T>
void
StructureModel<T>::printDeformationFluxIntegrals()
{
  _impl->printDeformationFluxIntegrals();
}

/*template<class T>
void
FlowModel<T>::printMassFluxIntegrals()
{
  _impl->printMassFluxIntegrals();
}

template<class T>
Vector<T,3>
FlowModel<T>::getPressureIntegral(const Mesh& mesh, const int faceGroupId)
{
 return  _impl->getPressureIntegral(mesh,faceGroupId);
}

template<class T>
Vector<T,3>
FlowModel<T>::getPVIntegral(const Field& velCoeff,const Mesh& mesh, const int faceGroupId)
{
  return  _impl->getPVIntegral(velCoeff,mesh,faceGroupId);
}
*/

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

template<class T>
void
StructureModel<T>::getTractionX(const Mesh& mesh)
{
  return  _impl->getTractionX(mesh);
}

/*
template<class T>
Vector<T,3>
FlowModel<T>::getPressureIntegralonIBFaces(const Mesh& mesh)
{
 return  _impl->getPressureIntegralonIBFaces(mesh);
}


template<class T>
Vector<T,3>
StructureModel<T>::getMomentumFluxIntegralonIBFaces(const Mesh& mesh)
{
 return  _impl->getDeformationFluxIntegralonIBFaces(mesh);
}


template<class T>
void
FlowModel<T>::computeIBFaceVelocity(const StorageSite& particles)
{
  return _impl->computeIBFaceVelocity(particles);
}


template<class T>
void
FlowModel<T>::computeIBandSolidVelocity(const StorageSite& particles)
{
  return _impl->computeIBandSolidVelocity(particles);
}
*/


template<class T>
boost::shared_ptr<ArrayBase>
StructureModel<T>::getStressTensor(const Mesh& mesh, const ArrayBase& cellIds)
{
  return _impl->getStressTensor(mesh, cellIds);
}

