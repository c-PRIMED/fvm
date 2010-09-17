class CG : public LinearSolver
{
public:

  CG();
  
  LinearSolver *preconditioner;
};

