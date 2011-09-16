class BCGStab : public LinearSolver
{
public:

  BCGStab();
  int getTotalIterations() const;
  LinearSolver *preconditioner;
};

