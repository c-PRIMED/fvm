%{
#include "DensityOfStates.h"
  %}

%include "std_vector.i"

template<class T>
class DensityOfStates
{
 public:
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef Array<T_Scalar> TArray;
  typedef shared_ptr<TArray> TArrPtr;
  typedef vector<TArrPtr> TArrList;
  typedef Array<int> IntArray;
  typedef shared_ptr<IntArray> IntArrayPtr;
  typedef vector<IntArrayPtr> IntArrList;
  typedef Kspace<T> Tkspace;

  DensityOfStates(const Tkspace& kspace);
  void binMode(const int mode, const int noBins, const T minw, const T maxw);
  void binEntireKspace(const int noBins, const T minw, const T maxw);
  void setDensity();
  void saveNormDOS(const char* filename);
  ArrayBase* makeDMMtransmission(DensityOfStates& otherDOS, const T Temp, const bool merge);
  ArrayBase* getFreqMids();
  ArrayBase* getFreqBins();
  T calcBinFlux(const T Temp, const int fBin, const T tau);
  ArrayBase* makeDMMreflection(ArrayBase* trans);
  
 private:
  TArray _FreqMids;
  TArray _FreqBounds;
  TArray _Density;
  IntArrList _BinKpts;
  IntArrList _BinModes;
  TArrList _ModeFractions;
  Tkspace& _kspace;
  
};

%template(DOS) DensityOfStates<ATYPE_STR>;
