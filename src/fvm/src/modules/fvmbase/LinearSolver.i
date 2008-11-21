
class LinearSolver
{
public:

  int nMaxIterations;
  int verbosity;
  double relativeTolerance;
  double absoluteTolerance;
  
private:

  LinearSolver();
  LinearSolver(const LinearSolver&);
};

