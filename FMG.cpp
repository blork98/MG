#include <FMG.h>
#include <LAOperation.h>

FMGSolver::FMGSolver(unsigned int numLevels, double tolerance,
	const std::vector<unsigned int>& sweeps,
	const std::vector<std::shared_ptr<SparseLinearSolver>>& solvers,
	const std::shared_ptr<InterpolationOperator>& I,
	const std::shared_ptr<RestrictionOperator>& R)
	:MultiGridSolver(numLevels, tolerance, sweeps, solvers, I, R)
{
}

void FMGSolver::create_rhs() const
{
	for (size_t rhsIt = 1; rhsIt < rhs_.size(); ++rhsIt)
	{
		std::vector<double>& rhsFine = *(rhs_[rhsIt-1]);
		std::vector<double>& rhsCoarse = *(rhs_[rhsIt]);

		R_->apply_operator(rhsFine, rhsCoarse);
	}
}

bool FMGSolver::solve(bool interiorPointsOnly, std::vector<double>& sol) const
{
	verify_inputs();

	//initialize rhs
	create_rhs();

	if (sol.size() != A_[0]->cols())
		sol.resize(A_[0]->cols(), 0.0);

	sols_[0] = sol;

	//solve on coarsest grid
	solvers_[numLevels_ - 1]->set_A(A_[numLevels_ - 1]);
	solvers_[numLevels_ - 1]->set_rhs(rhs_[numLevels_ - 1]);
	solvers_[numLevels_ - 1]->solve(sols_[numLevels_ - 1]);

	//for all levels above coarsest level, do a v-cycle on each
	for (int level = numLevels_- 2; level >= 0; --level)
	{
		//create initial guess of next level via interpolation of this level 
		I_->apply_operator(sols_[level], sols_[level+1]);

		//do V-Cycle starting at level
		v_cycle(level);
	}

	sol = sols_[0];

	//calculate final residual
	calc_residual(*rhs_[0], sols_[0], *A_[0], resids_[0]);
	double error = vector_norm_2(resids_[0]);

	if (error < tolerance_)
		return true;

	return false;
}

void FMGSolver::v_cycle(int start) const
{
	int startLevel = start;
	int endLevel = numLevels_-1;

	//down sweep
	for (int level = startLevel; level < (numLevels_ - 1); ++level)
	{
		//initial relaxation downsweeps
		solvers_[level]->set_A(A_[level]);
		solvers_[level]->set_rhs(rhs_[level]);
		solvers_[level]->solve(sols_[level]);

		//calculate residual
		calc_residual(*rhs_[level], sols_[level], *A_[level], resids_[level]);

		//apply restriction operator
		R_->apply_operator(resids_[level], *rhs_[level + 1]);
	}

	//Coarsest Level
	solvers_[numLevels_ - 1]->set_A(A_[numLevels_ - 1]);
	solvers_[numLevels_ - 1]->set_rhs(rhs_[numLevels_ - 1]);
	solvers_[numLevels_ - 1]->solve(sols_[numLevels_ - 1]);

	//up sweep
	for (int level = numLevels_ - 2; level >= startLevel; --level)
	{
		//apply interpolation operator
		I_->apply_operator(resids_[level], sols_[level + 1]);

		//perform correction
		vector_add(sols_[level], resids_[level], sols_[level]);

		//perform relaxation upsweep
		solvers_[level]->solve(sols_[level]);
	}

}