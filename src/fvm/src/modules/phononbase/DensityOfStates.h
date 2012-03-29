#ifndef _DENSITYOFSTATES_H_
#define _DENSITYOFSTATES_H_

#include <iostream>
#include <fstream>

template<class T>
class DensityOfStates
{

 public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef shared_ptr<TArray> TArrPtr;
  typedef vector<TArrPtr> TArrList;
  typedef Array<int> IntArray;
  typedef shared_ptr<IntArray> IntArrayPtr;
  typedef vector<IntArrayPtr> IntArrList;
  typedef Kspace<T> Tkspace;
  typedef kvol<T> Tkvol;
  typedef pmode<T> Tmode;

 DensityOfStates(const Tkspace& kspace):
  _FreqMids(0),
    _FreqBounds(0),
    _Density(0),
    _BinKpts(),
    _BinModes(),
    _ModeFractions(),
    _kspace(kspace)
    {}

  void binMode(const int mode, const int noBins, const T minw, const T maxw)
  {
    const T dwStar=1./T(noBins);
    const T Deltaw=maxw-minw;
    TArray modeBounds(noBins+1);
    TArray modeMids(noBins);
    IntArray BinKcount(noBins);
    BinKcount.zero();
    IntArrList modeKpts;

    for(int i=0;i<noBins+1;i++)
      modeBounds[i]=i*dwStar*Deltaw+minw;
    
    for(int i=0;i<noBins;i++)
      modeMids[i]=(modeBounds[i]+modeBounds[i+1])/2.;
    
    const int klen=_kspace.getlength();

    for(int k=0;k<klen;k++)
      {
	T omega=_kspace.getkvol(k).getmode(mode).getomega();
	int bin=findBin(modeBounds, omega);
	BinKcount[bin]++;
      }

    for(int i=0;i<noBins;i++)
      {
	IntArrayPtr newArrPtr=IntArrayPtr(new IntArray(BinKcount[i]));
	newArrPtr->zero();
	modeKpts.push_back(newArrPtr);
      }
    
    BinKcount.zero();
    
    for(int k=0;k<klen;k++)
      {
	T omega=_kspace.getkvol(k).getmode(mode).getomega();
	int bin=findBin(modeBounds, omega);
	IntArray& binArray=*modeKpts[bin];
	binArray[BinKcount[bin]]=k;
	BinKcount[bin]++;
      }
    
    addMode(modeBounds, modeMids, modeKpts, mode);

  }
  
  void copyDOS(DensityOfStates& otherDOS)
  {
    _FreqMids.resize(otherDOS.getFreqMidsT().getLength());
    _FreqMids.copyFrom(otherDOS.getFreqMidsT());
    _FreqBounds.resize(otherDOS.getFreqBinsT().getLength());
    _FreqMids.copyFrom(otherDOS.getFreqBinsT());
    _BinKpts=otherDOS.getKptsList();
    _BinModes=otherDOS.getMList();
    setDensity();
  }

  void setDensity()
  {
    _Density.resize(_FreqMids.getLength());
    _Density.zero();
    for(int i=0;i<_Density.getLength();i++)
      {
	IntArray& Kpts=*_BinKpts[i];
	for(int j=0;j<Kpts.getLength();j++)
	  {
	    const T dk3=_kspace.getkvol(Kpts[j]).getdk3();
	    _Density[i]+=dk3;
	  }
      }
  }

  const Tkspace& getKspace() {return _kspace;}
  
  void saveNormDOS(const char* filename)
  {
    const T wMax=_FreqBounds[_FreqBounds.getLength()-1];
    const T wMin=_FreqBounds[0];
    const T freqDif=wMax-wMin;
    T maxDen=0.0;
    
    for(int i=0;i<_FreqMids.getLength();i++)
      if(_Density[i]>maxDen)
	maxDen=_Density[i];
    
    ofstream DOSfile;
    DOSfile.open(filename);

    for(int i=0;i<_FreqMids.getLength();i++)
      DOSfile<<(_FreqMids[i]-wMin)/freqDif<<" "<<_Density[i]/maxDen<<endl;

    DOSfile.close();
  }
  
  ArrayBase* getFreqMids() {return &_FreqMids;}
  ArrayBase* getFreqBins() {return &_FreqBounds;}
  TArray& getFreqMidsT() {return _FreqMids;}
  TArray& getFreqBinsT() {return _FreqBounds;}

  int findBin(const TArray& bounds, const T freq)
  {
    for(int i=0;i<bounds.getLength()-1;i++)
      {
	if(freq>bounds[i] && freq<=bounds[i+1])
	  return i;
      }
    throw CException("Frequency not in given discretization!");
  }

  const int findBin(const T freq) const
  {
    for(int i=0;i<_FreqBounds.getLength()-1;i++)
      {
	if(freq>_FreqBounds[i] && freq<=_FreqBounds[i+1])
	  return i;
      }
    return -1;
  }
  
  IntArray& getKIndices(const int fBin) {return *_BinKpts[fBin];}
  IntArray& getMIndices(const int fBin) {return *_BinModes[fBin];}
  IntArrList& getKptsList() {return _BinKpts;}
  IntArrList& getMList() {return _BinModes;}
  TArray& getModeFractions(const int fBin) {return *_ModeFractions[fBin];}
  
  T sumOutgoing(const VectorT3& n, const int fBin, const T Temp)
  {
    T outSum(0.);
    
    if(fBin>=0)
      {
	const IntArray& Kpts=getKIndices(fBin);
	const IntArray& Mpts=getMIndices(fBin);
	for(int i=0;i<Kpts.getLength();i++)
	  {
	    const int k=Kpts[i];
	    const int m=Mpts[i];
	    const VectorT3 vg=_kspace.getkvol(k).getmode(m).getv();
	    const T vdot=vg[0]*n[0]+vg[1]*n[1]+vg[2]*n[2];
	    if(vdot>0.)
	      {
		const T dk3=_kspace.getkvol(k).getdk3();
		outSum+=dk3*vdot*_kspace.getkvol(k).getmode(m).calce0(Temp);
	      }
	  }
      }
    
    return outSum;
  }
  
  ArrayBase* makeDMMtransmission(DensityOfStates& otherDOS, const T Temp)
  {
    mergeBins(otherDOS, true);
    
    const int freqCount=_FreqMids.getLength();
    TArray* trans=new TArray(freqCount);
    VectorT3 n;

    n[0]=1.;
    n[1]=0.;
    n[2]=0.;

    for(int f0=0;f0<freqCount;f0++)
      {
	T w0=_FreqMids[f0];
	const int f1=otherDOS.findBin(w0);
	const T mat0sum=sumOutgoing(n,f0,Temp);
	const T mat1sum=otherDOS.sumOutgoing(-n,f1,Temp);
	(*trans)[f0]=mat1sum/(mat1sum+mat0sum);
      }
    return trans;
  }

  T calcBinFlux(const T Temp, const int fBin, const T tau)
  {
    T outSum(0.);
    VectorT3 n;
    n[0]=1.;
    n[1]=0.;
    n[2]=0.;
    
    if(fBin>=0)
      {
	const IntArray& Kpts=getKIndices(fBin);
	const IntArray& Mpts=getMIndices(fBin);
	
	for(int i=0;i<Kpts.getLength();i++)
	  {
	    const int k=Kpts[i];
	    const int m=Mpts[i];
	    const VectorT3 vg=_kspace.getkvol(k).getmode(m).getv();
	    const T vdot=vg[0]*n[0]+vg[1]*n[1]+vg[2]*n[2];
	    const T e0=_kspace.getkvol(k).getmode(m).calce0(Temp);
	    if(vdot>0.)
	      {
		const T dk3=_kspace.getkvol(k).getdk3();
		outSum+=dk3*vdot*e0*tau;
	      }
	  }
      }
    
    return outSum;
  }

  void mergeBins(DensityOfStates &otherDOS, const bool original)
  {
    DensityOfStates oldSelf(getKspace());
    oldSelf.copyDOS(*this);
    TArray old(_FreqBounds.getLength());
    TArray& mBounds=otherDOS.getFreqBinsT();
    old.copyFrom(_FreqBounds);
    int newBinCount(0);
    
    int i(0), j(0);
    T ceil=min(old[0], mBounds[0]);
    T o,m;
    
    while(i<old.getLength() || j<mBounds.getLength())
      {
	if(i<old.getLength())
	  o=old[i];
	else
	  o=mBounds[j]+1.;
	
	if(j<mBounds.getLength())
	  m=mBounds[j];
	else
	  m=old[i]+1.;
	
	if(o==m)
	  {
	    if(o>ceil)
	      {
		newBinCount++;
		ceil=o;
	      }
	    i++;
	    j++;
	  }
	else if( o>=ceil && o<m)
	  {
	    if(o>ceil)
	      {
		ceil=o;
		newBinCount++;
	      }
	    i++;
	  }
	else if(m>=ceil && m<o)
	  {
	    if(m>ceil)
	      {
		ceil=m;
		newBinCount++;
	      }
	    j++;
	  }
      }

    TArray newBounds(newBinCount+1);
    ceil=min(old[0], mBounds[0]);
    i=0;
    j=0;
    newBinCount=1;
    newBounds[0]=ceil;
    
    while(i<old.getLength() || j<mBounds.getLength())
      {
	if(i<old.getLength())
	  o=old[i];
	else
	  o=mBounds[j]+1.;
	
	if(j<mBounds.getLength())
	  m=mBounds[j];
	else
	  m=old[i]+1.;
	
	if(o==m)
	  {
	    if(o>ceil)
	      {
		newBounds[newBinCount]=o;
		newBinCount++;
	      }
	    i++;
	    j++;
	  }
	else if( o>=ceil && o<m)
	  {
	    if(o>ceil)
	      {
		ceil=o;
		newBounds[newBinCount]=o;
		newBinCount++;
	      }
	    i++;
	  }
	else if(m>=ceil && m<o)
	  {
	    if(m>ceil)
	      {
		ceil=m;
		newBounds[newBinCount]=m;
		newBinCount++;
	      }
	    j++;
	  }
      }

    _FreqBounds.resize(newBounds.getLength());
    _FreqBounds.copyFrom(newBounds);
    refineBins();
    setMids();

    IntArray BinPop(_FreqMids.getLength());
    BinPop.zero();
    int searchStart(0);

    for(int k=0; k<old.getLength()-1;k++)
      {
	const T oTop=old[k+1];
	const T oBot=old[k];
	IntArray& klist=*_BinKpts[k];
	const int plusKs=klist.getLength();
	
	for(int l=searchStart;l<_FreqBounds.getLength();l++)
	  {
	    if(_FreqBounds[l]>=oBot && _FreqBounds[l]<oTop)
	      BinPop[l]+=plusKs;
	    if(_FreqBounds[l]>=oTop)
	      {
		searchStart=l;
		break;
	      }
	  }
      }

    IntArrList oldBinKpts=_BinKpts;
    IntArrList oldBinModes=_BinModes;
    TArrList oldModeFractions=_ModeFractions;
    _BinKpts.resize(_FreqMids.getLength());
    _BinModes.resize(_FreqMids.getLength());
    _ModeFractions.resize(_FreqMids.getLength());

    for(int fInd=0;fInd<_FreqMids.getLength();fInd++)
      {
	const int len=BinPop[fInd];
	IntArrayPtr newKArray=IntArrayPtr(new IntArray(len));
	IntArrayPtr newMArray=IntArrayPtr(new IntArray(len));
	TArrPtr newFArray=TArrPtr(new TArray(len));
	_BinKpts[fInd]=newKArray;
	_BinModes[fInd]=newMArray;
	_ModeFractions[fInd]=newFArray;
      }

    BinPop.zero();
    searchStart=0;
    for(int k=0; k<old.getLength()-1;k++)
      {
	const T oTop=old[k+1];
	const T oBot=old[k];
	const T dif=oTop-oBot;
	IntArray& klist=*oldBinKpts[k];
	IntArray& mlist=*oldBinModes[k];
	TArray& fraclist=*oldModeFractions[k];
	
	for(int l=searchStart;l<_FreqBounds.getLength();l++)
	  {
	    if(_FreqBounds[l]>=oBot && _FreqBounds[l]<oTop)
	      {
		const T masterDif=_FreqBounds[l+1]-_FreqBounds[l];
		const T coeff=masterDif/dif;
		IntArray& newKArray=*_BinKpts[l];
		IntArray& newMArray=*_BinModes[l];
		TArray& newFArray=*_ModeFractions[l];
		
		for(int addInd=0;addInd<klist.getLength();addInd++)
		  {
		    newKArray[BinPop[l]]=klist[addInd];
		    newMArray[BinPop[l]]=mlist[addInd];
		    newFArray[BinPop[l]]=coeff*fraclist[addInd];
		    BinPop[l]++;
		  }
	      }
	    if(_FreqBounds[l]==oTop)
	      {
		searchStart=l;
		break;
	      }
	  }
      }

    setDensity();
    if(original)
      otherDOS.mergeBins(oldSelf,false);
    
  }

 private:

  void addMode(const TArray& mBounds, const TArray& mMids, const IntArrList& mKlist, const int mode)
  {
    if(_BinKpts.empty())
      {
	_FreqMids.resize(mMids.getLength());
	_FreqMids.copyFrom(mMids);
	_FreqBounds.resize(mBounds.getLength());
	_FreqBounds.copyFrom(mBounds);
	_BinKpts=mKlist;
	_BinModes.resize(mMids.getLength());
	_ModeFractions.resize(mMids.getLength());
	for(int i=0;i<mMids.getLength();i++)
	  {
	    IntArray& Karray=*_BinKpts[i];
	    IntArrayPtr MarrayPtr=IntArrayPtr(new IntArray(Karray.getLength()));
	    TArrPtr FracArrPtr=TArrPtr(new TArray(Karray.getLength()));
	    *MarrayPtr=mode;
	    *FracArrPtr=1.;
	    _BinModes[i]=MarrayPtr;
	    _ModeFractions[i]=FracArrPtr;
	  }
      }
    else
      {
	TArray old(_FreqBounds.getLength());
	old.copyFrom(_FreqBounds);
	int newBinCount(0);
	
	int i(0), j(0);
	T ceil=min(old[0], mBounds[0]);
	T o,m;

	while(i<old.getLength() || j<mBounds.getLength())
	  {
	    if(i<old.getLength())
	      o=old[i];
	    else
	      o=mBounds[j]+1.;

	    if(j<mBounds.getLength())
	      m=mBounds[j];
	    else
	      m=old[i]+1.;

	    if(o==m)
	      {
		if(o>ceil)
		  {
		    newBinCount++;
		    ceil=o;
		  }
		i++;
		j++;
	      }
	    else if( o>=ceil && o<m)
	      {
		if(o>ceil)
		  {
		    ceil=o;
		    newBinCount++;
		  }
		i++;
	      }
	    else if(m>=ceil && m<o)
	      {
		if(m>ceil)
		  {
		    ceil=m;
		    newBinCount++;
		  }
		j++;
	      }
	  }

	TArray newBounds(newBinCount+1);
	ceil=min(old[0], mBounds[0]);
	i=0;
	j=0;
	newBinCount=1;
	newBounds[0]=ceil;
	
	while(i<old.getLength() || j<mBounds.getLength())
	  {
	    if(i<old.getLength())
	      o=old[i];
	    else
	      o=mBounds[j]+1.;
	    
	    if(j<mBounds.getLength())
	      m=mBounds[j];
	    else
	      m=old[i]+1.;
	    
	    if(o==m)
	      {
		if(o>ceil)
		  {
		    newBounds[newBinCount]=o;
		    newBinCount++;
		  }
		i++;
		j++;
	      }
	    else if( o>=ceil && o<m)
	      {
		if(o>ceil)
		  {
		    ceil=o;
		    newBounds[newBinCount]=o;
		    newBinCount++;
		  }
		i++;
	      }
	    else if(m>=ceil && m<o)
	      {
		if(m>ceil)
		  {
		    ceil=m;
		    newBounds[newBinCount]=m;
		    newBinCount++;
		  }
		j++;
	      }
	  }

	_FreqBounds.resize(newBounds.getLength());
	_FreqBounds.copyFrom(newBounds);
	setMids();

	IntArray BinPop(_FreqMids.getLength());
	BinPop.zero();
	int searchStart(0);
	for(int k=0; k<mBounds.getLength()-1;k++)
	  {
	    const T mTop=mBounds[k+1];
	    const T mBot=mBounds[k];
	    IntArray& klist=*mKlist[k];
	    const int plusKs=klist.getLength();
	    
	    for(int l=searchStart;l<_FreqBounds.getLength();l++)
	      {
		if(_FreqBounds[l]>=mBot && _FreqBounds[l]<mTop)
		  BinPop[l]+=plusKs;
		if(_FreqBounds[l]==mTop)
		  {
		    searchStart=l;
		    break;
		  }
	      }
	  }

	searchStart=0;
	for(int k=0; k<old.getLength()-1;k++)
	  {
	    const T oTop=old[k+1];
	    const T oBot=old[k];
	    IntArray& klist=*_BinKpts[k];
	    const int plusKs=klist.getLength();
	    
	    for(int l=searchStart;l<_FreqBounds.getLength();l++)
	      {
		if(_FreqBounds[l]>=oBot && _FreqBounds[l]<oTop)
		  BinPop[l]+=plusKs;
		if(_FreqBounds[l]>=oTop)
		  {
		    searchStart=l;
		    break;
		  }
	      }
	  }
	
	IntArrList oldBinKpts=_BinKpts;
	IntArrList oldBinModes=_BinModes;
	TArrList oldModeFractions=_ModeFractions;
	_BinKpts.resize(_FreqMids.getLength());
	_BinModes.resize(_FreqMids.getLength());
	_ModeFractions.resize(_FreqMids.getLength());

	for(int fInd=0;fInd<_FreqMids.getLength();fInd++)
	  {
	    const int len=BinPop[fInd];
	    IntArrayPtr newKArray=IntArrayPtr(new IntArray(len));
	    IntArrayPtr newMArray=IntArrayPtr(new IntArray(len));
	    TArrPtr newFArray=TArrPtr(new TArray(len));
	    _BinKpts[fInd]=newKArray;
	    _BinModes[fInd]=newMArray;
	    _ModeFractions[fInd]=newFArray;
	  }

	searchStart=0;
	BinPop.zero();
	for(int k=0; k<mBounds.getLength()-1;k++)
	  {
	    const T mTop=mBounds[k+1];
	    const T mBot=mBounds[k];
	    const T dif=mTop-mBot;
	    IntArray& klist=*mKlist[k];
	    
	    for(int l=searchStart;l<_FreqBounds.getLength();l++)
	      {
		if(_FreqBounds[l]>=mBot && _FreqBounds[l]<mTop)
		  {
		    const T masterDif=_FreqBounds[l+1]-_FreqBounds[l];
		    const T coeff=masterDif/dif;
		    IntArray& newKArray=*_BinKpts[l];
		    IntArray& newMArray=*_BinModes[l];
		    TArray& newFArray=*_ModeFractions[l];
		    
		    for(int addInd=0;addInd<klist.getLength();addInd++)
		      {
			newKArray[BinPop[l]]=klist[addInd];
			newMArray[BinPop[l]]=mode;
			newFArray[BinPop[l]]=coeff;
			BinPop[l]++;
		      }

		  }
		if(_FreqBounds[l]==mTop)
		  {
		    searchStart=l;
		    break;
		  }
	      }
	  }

	searchStart=0;
	for(int k=0; k<old.getLength()-1;k++)
	  {
	    const T oTop=old[k+1];
	    const T oBot=old[k];
	    const T dif=oTop-oBot;
	    IntArray& klist=*oldBinKpts[k];
	    IntArray& mlist=*oldBinModes[k];
	    TArray& fraclist=*oldModeFractions[k];
	    
	    for(int l=searchStart;l<_FreqBounds.getLength();l++)
	      {
		if(_FreqBounds[l]>=oBot && _FreqBounds[l]<oTop)
		  {
		    const T masterDif=_FreqBounds[l+1]-_FreqBounds[l];
		    const T coeff=masterDif/dif;
		    IntArray& newKArray=*_BinKpts[l];
		    IntArray& newMArray=*_BinModes[l];
		    TArray& newFArray=*_ModeFractions[l];
		    
		    for(int addInd=0;addInd<klist.getLength();addInd++)
		      {
			newKArray[BinPop[l]]=klist[addInd];
			newMArray[BinPop[l]]=mlist[addInd];
			newFArray[BinPop[l]]=coeff*fraclist[addInd];
			BinPop[l]++;
		      }
		  }
		if(_FreqBounds[l]==oTop)
		  {
		    searchStart=l;
		    break;
		  }
	      }
	  }

      }
  }
  
  void setMids()
  {
    _FreqMids.resize(_FreqBounds.getLength()-1);
    for(int i=0;i<_FreqBounds.getLength()-1;i++)
      _FreqMids[i]=(_FreqBounds[i]+_FreqBounds[i+1])/2.;
  }

  void refineBins()
  {
    int newBinCount(0);
    TArray refined(_FreqBounds.getLength());
    int i(0);

    while(i<_FreqBounds.getLength()-1)
      {
	if(_FreqBounds[i+1]-_FreqBounds[i]<1.)
	  {
	    refined[newBinCount]=_FreqBounds[i];
	    newBinCount++;
	    i+=2;
	  }
	else
	  {
	    refined[newBinCount]=_FreqBounds[i];
	    newBinCount++;
	    i++;
	  }
      }
    refined[newBinCount]=_FreqBounds[_FreqBounds.getLength()-1];

    _FreqBounds.resize(newBinCount+1);
    for(int j=0;j<newBinCount+1;j++)
      _FreqBounds[j]=refined[j];
  }
  
  TArray _FreqMids;
  TArray _FreqBounds;
  TArray _Density;
  IntArrList _BinKpts;
  IntArrList _BinModes;
  TArrList _ModeFractions;
  const Tkspace& _kspace;
  
};



#endif
