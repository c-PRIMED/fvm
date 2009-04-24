class AMG : public LinearSolver
{
public:

  AMG();
  
  enum CycleType
    {
      V_CYCLE,
      W_CYCLE,
      F_CYCLE
    };
  
  int maxCoarseLevels;
  int nPreSweeps;
  int nPostSweeps;
  int coarseGroupSize;
  double weightRatioThreshold;
  CycleType cycleType;
private:
  AMG(const AMG&);
};

