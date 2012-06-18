#ifndef _COMETINTERFACE_H_
#define _COMETINTERFACE_H_


template<class T>
class COMETInterface
{
 public:
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef Array<int> IntArray;
  typedef Array<T> TArray;
  typedef KSConnectivity<T> TKConn;
  typedef vector<TKConn*> TKClist;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef Kspace<T> Tkspace;
  typedef Tkspace* TkspPtr;
  typedef vector<Tkspace*> TkspList;
  typedef kvol<T> Tkvol;
  typedef pmode<T> Tmode;
  typedef map<int,int> MeshKspaceMap;
  typedef typename StorageSite::ScatterMap ScatterMap;
  typedef typename StorageSite::CommonMap CommonMap;
  typedef pair<T,T> FreqBin;
  typedef shared_ptr<FreqBin> BinPtr;
  typedef vector<BinPtr> BinList;
  typedef DensityOfStates<T> DOST;
  
 COMETInterface(const MeshList& meshes, TkspList& klist, MeshKspaceMap& MKMap, PhononMacro& macro,
		const GeomFields& geomFields):
  _meshes(meshes),
    _KList(klist),
    _MeshKspaceMap(MKMap),
    _macro(macro),
    _geomFields(geomFields)
    {}

  void makeDMMcoeffs(COMETIC<T>& ic)
  {
    const int Mid0=ic.MeshID0;
    const int Mid1=ic.MeshID1;
    const int Fid0=ic.FgID0;
    const int Fid1=ic.FgID1;
    const Mesh& mesh0=*_meshes[Mid0];
    const Mesh& mesh1=*_meshes[Mid1];
    Tkspace& kspace0=*_KList[_MeshKspaceMap[Mid0]];
    Tkspace& kspace1=*_KList[_MeshKspaceMap[Mid1]];
    const int k0len=kspace0.gettotmodes();
    const int k1len=kspace1.gettotmodes();
    const StorageSite& faces0=mesh0.getFaceGroup(Fid0).site;
    const StorageSite& faces1=mesh1.getFaceGroup(Fid1).site;
    const StorageSite& cells0=mesh0.getCells();
    const StorageSite& cells1=mesh1.getCells();
    const CRConnectivity& faceCells0=mesh0.getFaceCells(faces0);
    const CRConnectivity& faceCells1=mesh1.getFaceCells(faces1);
    const int faceCount=faces0.getCount();
    TKClist& kclist01=ic.getKConnectivity(Fid0);
    TKClist& kclist10=ic.getKConnectivity(Fid1);
    const Field& faceAreaField=_geomFields.area;
    const Field& coordsField=_geomFields.coordinate;
    const VectorT3Array& faceAreaArray=dynamic_cast<const VectorT3Array&>(faceAreaField[faces0]);
    const VectorT3Array& coord0Array=dynamic_cast<const VectorT3Array&>(coordsField[cells0]);
    TArray& TL0Array=dynamic_cast<TArray&>(_macro.temperature[cells0]);
    TArray& TL1Array=dynamic_cast<TArray&>(_macro.temperature[cells1]);

    DOST& dos0=*kspace0.getDOSptr();
    DOST& dos1=*kspace1.getDOSptr();
    TArray& freqMids=dos0.getFreqMidsT();
    const T binNos=freqMids.getLength();
    const TArray& transArray01=kspace0.getTransArray(kspace1);

    ic.clearConnections();

    for(int f=0;f<faceCount;f++)
      {
	const VectorT3& An=faceAreaArray[f];
	VectorT3 n=An/sqrt(pow(An[0],2)+pow(An[1],2)+pow(An[2],2));
	const int cell0=faceCells0(f,0);
	const int cell0ghost=faceCells0(f,1);
	const int cell1=faceCells1(f,0);
	const int cell1ghost=faceCells1(f,1);
	const VectorT3 c0Pos=coord0Array[cell0];
	const VectorT3 cGhstPos=coord0Array[cell0ghost];
	const VectorT3 from0toGhst=cGhstPos-c0Pos;
	T TL1(TL1Array[cell1]), TL0(TL0Array[cell0]);

	if((n[0]*from0toGhst[0]+n[1]*from0toGhst[1]+n[2]*from0toGhst[2])<0.)
	  n*=-1.;

	findWallTemp(Mid0,cell0,cell0ghost,An,TL0);
	findWallTemp(Mid1,cell1,cell1ghost,-An,TL1);
	
	TKConn* conn01=new TKConn(k0len+1);  //must include lattice connection (reason for +1)
	TKConn* conn10=new TKConn(k1len+1);
	conn01->setColumnLength(k1len+1);
	conn10->setColumnLength(k0len+1);

	for(int binIndx=0;binIndx<binNos;binIndx++)
	  {
	    IntArray& k0kpts=dos0.getKIndices(binIndx);
	    IntArray& k1kpts=dos1.getKIndices(binIndx);
	    IntArray& k0mpts=dos0.getMIndices(binIndx);
	    IntArray& k1mpts=dos1.getMIndices(binIndx);
	    IntArray pos0self(k0kpts.getLength());
	    IntArray pos1self(k1kpts.getLength());
	    IntArray pos0other(k0kpts.getLength());
	    IntArray pos1other(k1kpts.getLength());
	    const T t01=transArray01[binIndx];
	    const T r01=1.-t01;
	    const T t10=r01;
	    const T r10=1-t10;

	    int out0(0);
	    int out1(0);
	    T in0sum(0.);
	    T in1sum(0.);
	    pos0self.zero();
	    pos1self.zero();
	    pos0other.zero();
	    pos1other.zero();

	    for(int kI=0;kI<k0kpts.getLength();kI++)
	      {
		Tkvol& kvol0=kspace0.getkvol(k0kpts[kI]);
		T dk30=kvol0.getdk3();
		Tmode& mode0=kvol0.getmode(k0mpts[kI]);
		VectorT3 vg0=mode0.getv();
		T v0dotN=vg0[0]*n[0]+vg0[1]*n[1]+vg0[2]*n[2];
		
		if(v0dotN>=0)
		  {
		    out0++;
		    int count0=mode0.getIndex()-1;
		    conn01->makeSelf(count0,1);
		    conn01->setSelf(count0,0,count0,1.);
		  }
		else
		  {
		    in0sum-=v0dotN*dk30*mode0.calce0(TL0);
		    //incp0sum+=mode0.calcde0dT(TL0)*dk30;
		  }
	      }

	    for(int kI=0;kI<k1kpts.getLength();kI++)
	      {
		Tkvol& kvol1=kspace1.getkvol(k1kpts[kI]);
		T dk31=kvol1.getdk3();
		Tmode& mode1=kvol1.getmode(k1mpts[kI]);
		VectorT3 vg1=mode1.getv();
		T v1dotN=vg1[0]*n[0]+vg1[1]*n[1]+vg1[2]*n[2];

		if(v1dotN<=0)
		  {
		    out1++;
		    int count1=mode1.getIndex()-1;
		    conn10->makeSelf(count1,1);
		    conn10->setSelf(count1,0,count1,1.);
		  }
		else
		  {
		    in1sum+=v1dotN*dk31*mode1.calce0(TL1);
		    //incp1sum+=mode1.calcde0dT(TL1)*dk31;
		  }
	      }

	    for(int kI=0;kI<k0kpts.getLength();kI++)
	      {
		Tkvol& kvol0=kspace0.getkvol(k0kpts[kI]);
		Tmode& mode0=kvol0.getmode(k0mpts[kI]);
		int count0=mode0.getIndex()-1;
		VectorT3 vg0=mode0.getv();
		T v0dotN=vg0[0]*n[0]+vg0[1]*n[1]+vg0[2]*n[2];
		
		if(v0dotN<0)
		  {
		    conn01->makeOther(count0,out1);
		    conn01->makeSelf(count0,out0);
		  }
	      }

	    for(int kI=0;kI<k1kpts.getLength();kI++)
	      {
		Tkvol& kvol1=kspace1.getkvol(k1kpts[kI]);
		Tmode& mode1=kvol1.getmode(k1mpts[kI]);
		int count1=mode1.getIndex()-1;
		VectorT3 vg1=mode1.getv();
		T v1dotN=vg1[0]*n[0]+vg1[1]*n[1]+vg1[2]*n[2];
		
		if(v1dotN>0)
		  {
		    conn10->makeOther(count1,out0);
		    conn10->makeSelf(count1,out1);
		  }

	      }

	    for(int kI=0;kI<k0kpts.getLength();kI++)
	      {
		Tkvol& kvol0=kspace0.getkvol(k0kpts[kI]);
		Tmode& mode0=kvol0.getmode(k0mpts[kI]);
		//T dk30=kvol0.getdk3();
		T de0dT0=mode0.calce0(TL0);
		int count0=mode0.getIndex()-1;
		VectorT3 vg0=mode0.getv();
		T v0dotN=vg0[0]*n[0]+vg0[1]*n[1]+vg0[2]*n[2];

		if(v0dotN<0)
		  {
		    for(int kkI=0;kkI<k0kpts.getLength();kkI++)
		      {
			Tkvol& kkvol0=kspace0.getkvol(k0kpts[kkI]);
			T ddk30=kkvol0.getdk3();
			Tmode& mmode0=kkvol0.getmode(k0mpts[kkI]);
			int ccount0=mmode0.getIndex()-1;
			VectorT3 vvg0=mmode0.getv();
			T vv0dotN=vvg0[0]*n[0]+vvg0[1]*n[1]+vvg0[2]*n[2];
			
			if(vv0dotN>=0)
			  {
			    T coeff=vv0dotN*ddk30/in0sum*r01*de0dT0;
			    conn01->setSelf(count0,pos0self[kI],ccount0,coeff);
			    pos0self[kI]++;
			  }
		      }

		    for(int kkI=0;kkI<k1kpts.getLength();kkI++)
		      {
			Tkvol& kkvol1=kspace1.getkvol(k1kpts[kkI]);
			T ddk31=kkvol1.getdk3();
			Tmode& mmode1=kkvol1.getmode(k1mpts[kkI]);
			int ccount1=mmode1.getIndex()-1;
			VectorT3 vvg1=mmode1.getv();
			T vv1dotN=vvg1[0]*n[0]+vvg1[1]*n[1]+vvg1[2]*n[2];
			
			if(vv1dotN<=0)
			  {
			    T coeff=-vv1dotN*ddk31/in0sum*t10*de0dT0;
			    conn01->setOther(count0,pos0other[kI],ccount1,coeff);
			    pos0other[kI]++;
			  }
		      }
		  }
	      }

	    for(int kI=0;kI<k1kpts.getLength();kI++)
	      {
		Tkvol& kvol1=kspace1.getkvol(k1kpts[kI]);
		Tmode& mode1=kvol1.getmode(k1mpts[kI]);
		//T dk31=kvol1.getdk3();
		T de0dT1=mode1.calce0(TL1);
		int count1=mode1.getIndex()-1;
		VectorT3 vg1=mode1.getv();
		T v1dotN=vg1[0]*n[0]+vg1[1]*n[1]+vg1[2]*n[2];

		if(v1dotN>0)
		  {
		    for(int kkI=0;kkI<k1kpts.getLength();kkI++)
		      {
			Tkvol& kkvol1=kspace1.getkvol(k1kpts[kkI]);
			T ddk31=kkvol1.getdk3();
			Tmode& mmode1=kkvol1.getmode(k1mpts[kkI]);
			int ccount1=mmode1.getIndex()-1;
			VectorT3 vvg1=mmode1.getv();
			T vv1dotN=vvg1[0]*n[0]+vvg1[1]*n[1]+vvg1[2]*n[2];
			
			if(vv1dotN<=0)
			  {
			    T coeff=-vv1dotN*ddk31/in1sum*r10*de0dT1;
			    conn10->setSelf(count1,pos1self[kI],ccount1,coeff);
			    pos1self[kI]++;
			  }
		      }

		    for(int kkI=0;kkI<k0kpts.getLength();kkI++)
		      {
			Tkvol& kkvol0=kspace0.getkvol(k0kpts[kkI]);
			T ddk30=kkvol0.getdk3();
			Tmode& mmode0=kkvol0.getmode(k0mpts[kkI]);
			int ccount0=mmode0.getIndex()-1;
			VectorT3 vvg0=mmode0.getv();
			T vv0dotN=vvg0[0]*n[0]+vvg0[1]*n[1]+vvg0[2]*n[2];
			
			if(vv0dotN>=0)
			  {
			    T coeff=vv0dotN*ddk30/in1sum*t01*de0dT1;
			    conn10->setOther(count1,pos1other[kI],ccount0,coeff);
			    pos1other[kI]++;
			  }
		      }
		  }
	      }
	  }

	kclist01.push_back(conn01);
	kclist10.push_back(conn10);
      }

  }

  void remakeDMMcoeffs(COMETIC<T>& ic)
  {
    const int Mid0=ic.MeshID0;
    const int Mid1=ic.MeshID1;
    const int Fid0=ic.FgID0;
    const int Fid1=ic.FgID1;
    const Mesh& mesh0=*_meshes[Mid0];
    const Mesh& mesh1=*_meshes[Mid1];
    Tkspace& kspace0=*_KList[_MeshKspaceMap[Mid0]];
    Tkspace& kspace1=*_KList[_MeshKspaceMap[Mid1]];
    const StorageSite& faces0=mesh0.getFaceGroup(Fid0).site;
    const StorageSite& faces1=mesh1.getFaceGroup(Fid1).site;
    const StorageSite& cells0=mesh0.getCells();
    const StorageSite& cells1=mesh1.getCells();
    const CRConnectivity& faceCells0=mesh0.getFaceCells(faces0);
    const CRConnectivity& faceCells1=mesh1.getFaceCells(faces1);
    const int faceCount=faces0.getCount();
    TKClist& kclist01=ic.getKConnectivity(Fid0);
    TKClist& kclist10=ic.getKConnectivity(Fid1);
    const Field& faceAreaField=_geomFields.area;
    const Field& coordsField=_geomFields.coordinate;
    const VectorT3Array& faceAreaArray=dynamic_cast<const VectorT3Array&>(faceAreaField[faces0]);
    const VectorT3Array& coord0Array=dynamic_cast<const VectorT3Array&>(coordsField[cells0]);
    TArray& TL0Array=dynamic_cast<TArray&>(_macro.temperature[cells0]);
    TArray& TL1Array=dynamic_cast<TArray&>(_macro.temperature[cells1]);

    DOST& dos0=*kspace0.getDOSptr();
    DOST& dos1=*kspace1.getDOSptr();
    TArray& freqMids=dos0.getFreqMidsT();
    const T binNos=freqMids.getLength();
    const TArray& transArray01=kspace0.getTransArray(kspace1);

    for(int f=0;f<faceCount;f++)
      {
	const VectorT3& An=faceAreaArray[f];
	VectorT3 n=An/sqrt(pow(An[0],2)+pow(An[1],2)+pow(An[2],2));
	const int cell0=faceCells0(f,0);
	const int cell0ghost=faceCells0(f,1);
	const int cell1=faceCells1(f,0);
	const int cell1ghost=faceCells1(f,1);
	const VectorT3 c0Pos=coord0Array[cell0];
	const VectorT3 cGhstPos=coord0Array[cell0ghost];
	const VectorT3 from0toGhst=cGhstPos-c0Pos;
	T TL1(TL1Array[cell1]), TL0(TL0Array[cell0]);

	if((n[0]*from0toGhst[0]+n[1]*from0toGhst[1]+n[2]*from0toGhst[2])<0.)
	  n*=-1.;

	findWallTemp(Mid0,cell0,cell0ghost,An,TL0);
	findWallTemp(Mid1,cell1,cell1ghost,-An,TL1);
	
	TKConn* conn01=kclist01[f];
	TKConn* conn10=kclist10[f];

	for(int binIndx=0;binIndx<binNos;binIndx++)
	  {
	    IntArray& k0kpts=dos0.getKIndices(binIndx);
	    IntArray& k1kpts=dos1.getKIndices(binIndx);
	    IntArray& k0mpts=dos0.getMIndices(binIndx);
	    IntArray& k1mpts=dos1.getMIndices(binIndx);
	    IntArray pos0self(k0kpts.getLength());
	    IntArray pos1self(k1kpts.getLength());
	    IntArray pos0other(k0kpts.getLength());
	    IntArray pos1other(k1kpts.getLength());

	    /*
	    const T t01=transArray01[binIndx];
	    const T r01=1.-t01;
	    const T t10=r01;
	    const T r10=1-t10;
	    */

	    VectorT3 normal;
	    normal[0]=1.;
	    normal[0]=0.;
	    normal[0]=0.;
	    
	    T eout0T0=dos0.sumOutgoing(normal,binIndx,(TL0+TL1)/2.);
	    T eout1T0=dos1.sumOutgoing(-normal,binIndx,(TL0+TL1)/2.);
	    const T t01=1./(1.+eout0T0/eout1T0);
	    const T r01=1-t01;

	    //T eout0T1=dos0.sumOutgoing(normal,binIndx,TL1);
	    //T eout1T1=dos1.sumOutgoing(-normal,binIndx,TL1);
	    const T t10=r01;//1./(1.+eout1T1/eout0T1);
	    const T r10=1-t10;

	    T in0sum(0.);
	    T in1sum(0.);
	    pos0self.zero();
	    pos1self.zero();
	    pos0other.zero();
	    pos1other.zero();

	    for(int kI=0;kI<k0kpts.getLength();kI++)
	      {
		Tkvol& kvol0=kspace0.getkvol(k0kpts[kI]);
		T dk30=kvol0.getdk3();
		Tmode& mode0=kvol0.getmode(k0mpts[kI]);
		VectorT3 vg0=mode0.getv();
		T v0dotN=vg0[0]*n[0]+vg0[1]*n[1]+vg0[2]*n[2];
		
		if(v0dotN<0)
		  in0sum-=v0dotN*dk30*mode0.calce0(TL0);

	      }

	    for(int kI=0;kI<k1kpts.getLength();kI++)
	      {
		Tkvol& kvol1=kspace1.getkvol(k1kpts[kI]);
		T dk31=kvol1.getdk3();
		Tmode& mode1=kvol1.getmode(k1mpts[kI]);
		VectorT3 vg1=mode1.getv();
		T v1dotN=vg1[0]*n[0]+vg1[1]*n[1]+vg1[2]*n[2];

		if(v1dotN>0)
		  in1sum+=v1dotN*dk31*mode1.calce0(TL1);

	      }


	    for(int kI=0;kI<k0kpts.getLength();kI++)
	      {
		Tkvol& kvol0=kspace0.getkvol(k0kpts[kI]);
		Tmode& mode0=kvol0.getmode(k0mpts[kI]);
		T de0dT0=mode0.calce0(TL0);
		int count0=mode0.getIndex()-1;
		VectorT3 vg0=mode0.getv();
		T v0dotN=vg0[0]*n[0]+vg0[1]*n[1]+vg0[2]*n[2];

		if(v0dotN<0)
		  {
		    for(int kkI=0;kkI<k0kpts.getLength();kkI++)
		      {
			Tkvol& kkvol0=kspace0.getkvol(k0kpts[kkI]);
			T ddk30=kkvol0.getdk3();
			Tmode& mmode0=kkvol0.getmode(k0mpts[kkI]);
			int ccount0=mmode0.getIndex()-1;
			VectorT3 vvg0=mmode0.getv();
			T vv0dotN=vvg0[0]*n[0]+vvg0[1]*n[1]+vvg0[2]*n[2];
			
			if(vv0dotN>=0)
			  {
			    T coeff=vv0dotN*ddk30/in0sum*r01*de0dT0;
			    conn01->resetSelf(count0,pos0self[kI],ccount0,coeff);
			    pos0self[kI]++;
			  }
		      }

		    for(int kkI=0;kkI<k1kpts.getLength();kkI++)
		      {
			Tkvol& kkvol1=kspace1.getkvol(k1kpts[kkI]);
			T ddk31=kkvol1.getdk3();
			Tmode& mmode1=kkvol1.getmode(k1mpts[kkI]);
			int ccount1=mmode1.getIndex()-1;
			VectorT3 vvg1=mmode1.getv();
			T vv1dotN=vvg1[0]*n[0]+vvg1[1]*n[1]+vvg1[2]*n[2];
			
			if(vv1dotN<=0)
			  {
			    T coeff=-vv1dotN*ddk31/in0sum*t10*de0dT0;
			    conn01->resetOther(count0,pos0other[kI],ccount1,coeff);
			    pos0other[kI]++;
			  }
		      }
		  }
	      }

	    for(int kI=0;kI<k1kpts.getLength();kI++)
	      {
		Tkvol& kvol1=kspace1.getkvol(k1kpts[kI]);
		Tmode& mode1=kvol1.getmode(k1mpts[kI]);
		T de0dT1=mode1.calce0(TL1);
		int count1=mode1.getIndex()-1;
		VectorT3 vg1=mode1.getv();
		T v1dotN=vg1[0]*n[0]+vg1[1]*n[1]+vg1[2]*n[2];

		if(v1dotN>0)
		  {
		    for(int kkI=0;kkI<k1kpts.getLength();kkI++)
		      {
			Tkvol& kkvol1=kspace1.getkvol(k1kpts[kkI]);
			T ddk31=kkvol1.getdk3();
			Tmode& mmode1=kkvol1.getmode(k1mpts[kkI]);
			int ccount1=mmode1.getIndex()-1;
			VectorT3 vvg1=mmode1.getv();
			T vv1dotN=vvg1[0]*n[0]+vvg1[1]*n[1]+vvg1[2]*n[2];
			
			if(vv1dotN<=0)
			  {
			    T coeff=-vv1dotN*ddk31/in1sum*r10*de0dT1;
			    conn10->resetSelf(count1,pos1self[kI],ccount1,coeff);
			    pos1self[kI]++;
			  }
		      }

		    for(int kkI=0;kkI<k0kpts.getLength();kkI++)
		      {
			Tkvol& kkvol0=kspace0.getkvol(k0kpts[kkI]);
			T ddk30=kkvol0.getdk3();
			Tmode& mmode0=kkvol0.getmode(k0mpts[kkI]);
			int ccount0=mmode0.getIndex()-1;
			VectorT3 vvg0=mmode0.getv();
			T vv0dotN=vvg0[0]*n[0]+vvg0[1]*n[1]+vvg0[2]*n[2];
			
			if(vv0dotN>=0)
			  {
			    T coeff=vv0dotN*ddk30/in1sum*t01*de0dT1;
			    conn10->resetOther(count1,pos1other[kI],ccount0,coeff);
			    pos1other[kI]++;
			  }
		      }
		  }
	      }
	  }

      }

  }

  void updateOtherGhost(const COMETIC<T>& ic, const int Mid0, const bool plusFAS)
  {//always going to be scattering the values of mesh with id0

    int Mid1,Fid0,Fid1;

    if(Mid0==ic.MeshID1)
      {
	Mid1=ic.MeshID0;
	Fid1=ic.FgID0;
	Fid0=ic.FgID1;
      }
    else
      {
	Mid1=ic.MeshID1;
	Fid1=ic.FgID1;
	Fid0=ic.FgID0;
      }

    const Mesh& mesh0=*_meshes[Mid0];
    const Mesh& mesh1=*_meshes[Mid1];
    Tkspace& kspace0=*_KList[_MeshKspaceMap[Mid0]];
    Tkspace& kspace1=*_KList[_MeshKspaceMap[Mid1]];    
    const StorageSite& faces0=mesh0.getFaceGroup(Fid0).site;
    const StorageSite& faces1=mesh1.getFaceGroup(Fid1).site;
    const ScatterMap& ScatMap0=faces0.getScatterMap();
    const IntArray& scatter01=*ScatMap0.find(&faces1)->second;
    const CRConnectivity& faceCells0=mesh0.getFaceCells(faces0);
    const CRConnectivity& faceCells1=mesh1.getFaceCells(faces1);
    const TKClist& KconnList1=ic.getKConnectivity(Fid1);
    const int numFaces=faces0.getCount();
    const int klen0=kspace0.gettotmodes()+1;
    const int klen1=kspace1.gettotmodes()+1;

    TArray transmitted(klen1);
    TArray reflected(klen1);
    TArray vals0(klen0);
    TArray vals1(klen1);

    for(int f0=0;f0<numFaces;f0++)
      {
	int cell0=faceCells0(f0,0);
	int cell1ghost=scatter01[f0];
	int cell1=faceCells1(f0,0);
	const TKConn& Kconn1=*KconnList1[f0];
	
	vals1.zero();
	vals0.zero();
	
	makeValueArray(Mid0,cell0,vals0);
	makeValueArray(Mid1,cell1,vals1);
	
	Kconn1.multiplyOther(vals0,transmitted);
	Kconn1.multiplySelf(vals1,reflected);
	transmitted+=reflected;
	if(plusFAS)
	  kspace1.addFASint(cell1ghost,transmitted);
	Distribute(Mid1,cell1ghost,transmitted);	
      }
  }

  void updateResid(const COMETIC<T>& ic, const bool plusFAS)
  {
    int Mid0,Mid1,Fid0,Fid1;

    Mid0=ic.MeshID0;
    Mid1=ic.MeshID1;
    Fid1=ic.FgID1;
    Fid0=ic.FgID0;
    
    const Mesh& mesh0=*_meshes[Mid0];
    const Mesh& mesh1=*_meshes[Mid1];
    Tkspace& kspace0=*_KList[_MeshKspaceMap[Mid0]];
    Tkspace& kspace1=*_KList[_MeshKspaceMap[Mid1]];    
    const StorageSite& faces0=mesh0.getFaceGroup(Fid0).site;
    const StorageSite& faces1=mesh1.getFaceGroup(Fid1).site;
    const ScatterMap& ScatMap0=faces0.getScatterMap();
    const IntArray& scatter01=*ScatMap0.find(&faces1)->second;
    const CRConnectivity& faceCells0=mesh0.getFaceCells(faces0);
    const CRConnectivity& faceCells1=mesh1.getFaceCells(faces1);
    const TKClist& KconnList1=ic.getKConnectivity(Fid1);
    const TKClist& KconnList0=ic.getKConnectivity(Fid0);
    const int numFaces=faces0.getCount();
    const int klen0=kspace0.gettotmodes()+1;
    const int klen1=kspace1.gettotmodes()+1;

    TArray transmitted01(klen1);
    TArray reflected10(klen1);
    TArray transmitted10(klen0);
    TArray reflected01(klen0);
    TArray vals0(klen0);
    TArray vals1(klen1);
    TArray currentSol0(klen0);
    TArray currentSol1(klen1);

    for(int f0=0;f0<numFaces;f0++)
      {
	int cell0=faceCells0(f0,0);
	int cell0ghost=faceCells0(f0,1);
	int cell1ghost=scatter01[f0];
	int cell1=faceCells1(f0,0);
	const TKConn& Kconn1=*KconnList1[f0];
	const TKConn& Kconn0=*KconnList0[f0];
	
	vals1.zero();
	vals0.zero();
	
	makeValueArray(Mid0,cell0,vals0);
	makeValueArray(Mid1,cell1,vals1);
	makeValueArray(Mid0,cell0ghost,currentSol0);
	makeValueArray(Mid1,cell1ghost,currentSol1);
	
	Kconn1.multiplyOther(vals0,transmitted01);
	Kconn1.multiplySelf(vals1,reflected10);
	transmitted01+=reflected10;
	transmitted01-=currentSol1; //now transmitted01 is the residual
	if(plusFAS)
	  kspace1.addFASint(cell1ghost,transmitted01);
	DistributeResid(Mid1,cell1ghost,transmitted01);

	Kconn0.multiplyOther(vals1,transmitted10);
	Kconn0.multiplySelf(vals0,reflected01);
	transmitted10+=reflected01;
	transmitted10-=currentSol0; //now transmitted10 is the residual
	if(plusFAS)
	  kspace0.addFASint(cell0ghost,transmitted10);
	DistributeResid(Mid0,cell0ghost,transmitted10);
      }
  }

  void makeValueArray(const int msh, const int c, TArray& o)
  {
    const Mesh& mesh=*_meshes[msh];
    const StorageSite& cells=mesh.getCells();
    Tkspace& kspace=*_KList[_MeshKspaceMap[msh]];
    int klen=kspace.gettotmodes();
    int kpts=kspace.getlength();

    if(klen+1==o.getLength())
      {
	kspace.geteCellVals(c,o);
	o[klen]=0.;
      }
    else
      throw CException("makeValueArray: Array not the same size as the k-space!");
  }

  void makeEquilibriumArray(const int msh, const T Temp, TArray& o)
  {
    const Mesh& mesh=*_meshes[msh];
    Tkspace& kspace=*_KList[_MeshKspaceMap[msh]];
    int klen=kspace.gettotmodes();
    int kpts=kspace.getlength();
    
    if(klen+1==o.getLength())
      {
	for(int k=0;k<kpts;k++)
	  {
	    Tkvol& kvol=kspace.getkvol(k);
	    const int numModes=kvol.getmodenum();
	    for(int m=0;m<numModes;m++)
	      {
		Tmode& mode=kvol.getmode(m);
		const int count=mode.getIndex()-1;
		o[count]=mode.calce0(Temp);
	      }
	  }
	o[klen]=0.;
      }
    else
      throw CException("makeEquilibriumArray: Array not the same size as the k-space!");
  }
  
  void Distribute(const int msh, const int cell, TArray& BVec)
  {
    const Mesh& mesh=*_meshes[msh];
    const StorageSite& cells=mesh.getCells();
    Tkspace& kspace=*_KList[_MeshKspaceMap[msh]];
    const int klen=kspace.gettotmodes();
    const int kpts=kspace.getlength();
    
    if(klen+1==BVec.getLength())
      {
	kspace.seteCellVals(cell, BVec);
	TArray& TlArray=dynamic_cast<TArray&>(_macro.temperature[cells]);
	TlArray[cell]=BVec[klen];
      }
    else
      throw CException("Distribute: Array not the same size as the k-space!");
  }

  void DistributeResid(const int msh, const int cell, TArray& BVec)
  {
    const Mesh& mesh=*_meshes[msh];
    const StorageSite& cells=mesh.getCells();
    Tkspace& kspace=*_KList[_MeshKspaceMap[msh]];
    const int klen=kspace.gettotmodes();
    const int kpts=kspace.getlength();
    
    if(klen+1==BVec.getLength())
      {
	kspace.setResidCell(cell, BVec);
      }
    else
      throw CException("DistributeResid: Array not the same size as the k-space!");
  }

  void ZeroGhost(const int msh, const int cell)
  {
    const Mesh& mesh=*_meshes[msh];
    const StorageSite& cells=mesh.getCells();
    Tkspace& kspace=*_KList[_MeshKspaceMap[msh]];
    const int klen=kspace.getlength();
    
    for(int k=0;k<klen;k++)
      {
	Tkvol& kvol=kspace.getkvol(k);
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    Field& efield=mode.getfield();
	    TArray& eArray=dynamic_cast<TArray&>(efield[cells]);
	    eArray[cell]=0.;
	  }
      }
    
    TArray& TlArray=dynamic_cast<TArray&>(_macro.temperature[cells]);
    TlArray[cell]=0;
  }

  void findWallTemp(const int meshID, const int cell0, const int cell0ghost, const VectorT3 Af, T& Tguess)
  {
    const Mesh& mesh=*_meshes[meshID];
    const StorageSite& cells=mesh.getCells();
    Tkspace& kspace=*_KList[_MeshKspaceMap[meshID]];
    int kpts=kspace.getlength();
    const VectorT3Array& coordArray=dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);

    const VectorT3 c0Pos=coordArray[cell0];
    const VectorT3 cGhstPos=coordArray[cell0ghost];
    const VectorT3 from0toGhst=cGhstPos-c0Pos;

    VectorT3 n=Af/sqrt(pow(Af[0],2)+pow(Af[1],2)+pow(Af[2],2));
    
    if((n[0]*from0toGhst[0]+n[1]*from0toGhst[1]+n[2]*from0toGhst[2])<0.)
      n*=-1.;
    
    T sume(0.);
    TArray eArray(kspace.gettotmodes());
    kspace.geteCellVals(cell0, eArray);

    for(int k=0;k<kpts;k++)
      {
	Tkvol& kvol=kspace.getkvol(k);
	const int numModes=kvol.getmodenum();
	T dk3=kvol.getdk3();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    VectorT3 vg=mode.getv();
	    T VdotN=vg[0]*n[0]+vg[1]*n[1]+vg[2]*n[2];
	    if(VdotN>0)
	      {
		int index=mode.getIndex()-1;
		sume+=eArray[index]*dk3;
	      }
	  }
      }

    kspace.calcTemp(Tguess,sume,n);
  }

  void addFAS(const int msh, const int c, TArray& bVec)
  {
    const Mesh& mesh=*_meshes[msh];
    const StorageSite& cells=mesh.getCells();
    Tkspace& kspace=*_KList[_MeshKspaceMap[msh]];
    const int klen=kspace.gettotmodes();
    const int kpts=kspace.getlength();

    if(klen+1==bVec.getLength())
      {
	for(int k=0;k<kpts;k++)
	  {
	    Tkvol& kvol=kspace.getkvol(k);
	    const int numModes=kvol.getmodenum();
	    for(int m=0;m<numModes;m++)
	      {
		Tmode& mode=kvol.getmode(m);
		const int count=mode.getIndex()-1;
		Field& fasField=mode.getFASfield();
		TArray& fasArray=dynamic_cast<TArray&>(fasField[cells]);
		bVec[count]+=fasArray[c];
	      }
	  }
	TArray& fasArray=dynamic_cast<TArray&>(_macro.TlFASCorrection[cells]);
	bVec[klen]+=fasArray[c];
      }
  }
  
 private:
  const MeshList& _meshes;
  TkspList& _KList;
  MeshKspaceMap _MeshKspaceMap;
  PhononMacro& _macro;
  const GeomFields& _geomFields;

};


#endif
