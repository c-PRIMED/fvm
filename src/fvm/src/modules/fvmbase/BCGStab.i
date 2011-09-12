class BCGStab : public LinearSolver
{
public:

  BCGStab();
  int getTotalIterations() const { return _totalIterations;}
  LinearSolver *preconditioner;
};

