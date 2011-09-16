class AMG : public LinearSolver
{
public:

  AMG();

  void setMergeLevelSize( int level_size);
  void redirectPrintToFile( const string& fname );
  void redirectPrintToScreen( );
  int getTotalIterations() const;
  enum CycleType
    {
      V_CYCLE,
      W_CYCLE,
      F_CYCLE
    };
  
  enum SmootherType
    {
      GAUSS_SEIDEL,
      JACOBI
    };

  int maxCoarseLevels;
  int nPreSweeps;
  int nPostSweeps;
  int coarseGroupSize;
  double weightRatioThreshold;
  CycleType cycleType;
  SmootherType smootherType;
  bool scaleCorrections;
private:
  AMG(const AMG&);
};

