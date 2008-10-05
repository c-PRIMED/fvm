
class AMG
{
public:

  enum CycleType
    {
      V_CYCLE,
      W_CYCLE,
      F_CYCLE
    };
  
  // these parameters can be tuned.
  int nMaxCycles;
  int maxCoarseLevels;
  int nPreSweeps;
  int nPostSweeps;
  int coarseGroupSize;
  double weightRatioThreshold;
  CycleType cycleType;
  int verbosity;
  double relativeTolerance;
  double absoluteTolerance;
  
private:

  AMG();
  AMG(const AMG&);
};

