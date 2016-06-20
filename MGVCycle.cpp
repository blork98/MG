#include <MGVCycle.h>
#include <LAOperation.h>

MGSolverVCycle::MGSolverVCycle(unsigned int numLevels, double tolerance,
	const std::vector<unsigned int>& sweeps,
	const std::vector<std::shared_ptr<SparseLinearSolver>>& solvers,
	const std::shared_ptr<InterpolationOperator>& I,
	const std::shared_ptr<RestrictionOperator>& R)
	:MultiGridSolver(numLevels, tolerance, sweeps, solvers,I,R)
{
}

bool MGSolverVCycle::solve(bool interiorPointsOnly, std::vector<double>& sol) const
{
	verify_inputs();

	if (sol.size() != A_[0]->cols())
		sol.resize(A_[0]->cols(), 0.0);

	sols_[0] = sol;

	//down sweep
	for (unsigned int level = 0; level < (numLevels_-1); ++level)
	{
		//initial relaxation downsweeps
		solvers_[level]->set_A(A_[level]);
		solvers_[level]->set_rhs(rhs_[level]);
		solvers_[level]->solve(sols_[level]);

		//calculate residual
		calc_residual(*rhs_[level], sols_[level], *A_[level], resids_[level]);

		//apply restriction operator
		I_->apply_operator(resids_[level], *rhs_[level+1]);
	}

	//Coarsest Level
	solvers_[numLevels_ - 1]->set_A(A_[numLevels_ - 1]);
	solvers_[numLevels_ - 1]->set_rhs(rhs_[numLevels_ - 1]);
	solvers_[numLevels_ - 1]->solve(sols_[numLevels_ - 1]);

	//up sweep
	unsigned int levelCtr = 0;
	for (unsigned int level = 0; level < (numLevels_-1); ++level) //error is here
	{
		levelCtr = numLevels_ - level - 1;

		//apply interpolation operator
		I_->apply_operator(resids_[levelCtr],sols_[levelCtr+1]);

		//perform correction
		vector_add(sols_[levelCtr],resids_[levelCtr],sols_[levelCtr]);

		//perform relaxation upsweep
		solvers_[levelCtr]->solve(sols_[levelCtr]);
	}

	sol = sols_[0];

	//calculate final residual
	calc_residual(*rhs_[0], sols_[0], *A_[0], resids_[0]);
	double error = vector_norm_2(resids_[0]);

	if (error < tolerance_)
		return true;

	return false;
}