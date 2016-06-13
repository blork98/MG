#include <MGTwoGridCorrection.h>
#include <LAOperation.h>

#include <assert.h>

TwoGridCorrection::TwoGridCorrection(unsigned int numRelaxSweeps,
	const std::shared_ptr<SparseLinearSolver>& coarseGridSolver,
	const std::shared_ptr<SparseIterativeLinearSolver> fineGridSolver)
	: numRelaxSweeps_(numRelaxSweeps), fineA_(nullptr), coarseA_(nullptr), 
	rhs_(nullptr), R_(nullptr), I_(nullptr),
	coarseGridSolver_(coarseGridSolver), fineGridSolver_(fineGridSolver)
{}

void TwoGridCorrection::set_coarse_grid_solver(
	const std::shared_ptr<SparseLinearSolver>& coarseGridSolver)
{
	coarseGridSolver_ = coarseGridSolver;
}

void TwoGridCorrection::set_fine_grid_solver(
	const std::shared_ptr<SparseIterativeLinearSolver>& fineGridSolver)
{
	fineGridSolver_ = fineGridSolver;
}

void TwoGridCorrection::set_A_fine(const std::shared_ptr<SparseMatrix>& A)
{
	fineA_ = A;
}

void TwoGridCorrection::set_A_coarse(const std::shared_ptr<SparseMatrix>& A)
{
	coarseA_ = A;
}

void TwoGridCorrection::set_rhs(const std::shared_ptr<std::vector<double>>& rhs)
{
	rhs_ = rhs;
}

void TwoGridCorrection::set_restriction(const std::shared_ptr<RestrictionOperator>& R)
{
	R_ = R;
}
void TwoGridCorrection::set_interpolation(const std::shared_ptr<InterpolationOperator>& I)
{
	I_ = I;
}

void TwoGridCorrection::verify_inputs() const
{
	//make sure none of the pointers are null
	assert(fineA_ != nullptr);
	assert(coarseA_ != nullptr);
	assert(rhs_ != nullptr);
	assert(coarseGridSolver_ != nullptr);
	assert(fineGridSolver_ != nullptr);
	assert(R_ != nullptr);
	assert(I_ != nullptr);

	//dimension checking
	assert(fineA_->cols() == fineA_->rows());
	assert(coarseA_->cols() == coarseA_->rows());
	assert(fineA_->cols() == rhs_->size());

}

const std::shared_ptr<SparseMatrix>& TwoGridCorrection::A_fine() const
{
	return fineA_;
}

const std::shared_ptr<SparseMatrix>& TwoGridCorrection::A_coarse() const
{
	return coarseA_;
}

const std::shared_ptr<std::vector<double>>& TwoGridCorrection::rhs() const
{
	return rhs_;
}

bool TwoGridCorrection::solve(bool interiorPointsOnly,
	std::vector<double>& sol) const
{
	verify_inputs();
	assert(sol.size() == fineA_->cols());

	size_t h = fineA_->cols() - 1;
	size_t h2 = coarseA_->cols() - 1;
	assert((h / 2) == h2);

	std::vector<double> coarseSol(h2+1,0.0);
	std::vector<double>& fineGridSol = sol;
	std::vector<double> coarseGridResidual(coarseSol.size(), 0.0),
		coarseGridSol(coarseSol.size(), 0.0);
	std::vector<double> fineGridResid(sol.size(), 0.0);

	//Initial Relaxation Sweep
	fineGridSolver_->set_A(fineA_);
	fineGridSolver_->set_rhs(rhs_);
	fineGridSolver_->set_max_iters(numRelaxSweeps_);
	fineGridSolver_->solve(fineGridSol);

	//Compute Fine Grid Residual
	fineA_->vec_multiply(fineGridSol, fineGridResid);
	vector_subtract(*rhs_, fineGridResid, fineGridResid);

	//Apply Restriction Operator to fine grid residual 
	I_->apply_operator(fineGridResid, coarseGridResidual);

	//Solve inner system
	std::shared_ptr<std::vector<double>> cGridRes; 
	*cGridRes = coarseGridResidual;

	coarseGridSolver_->set_A(coarseA_);
	coarseGridSolver_->set_rhs(cGridRes);
	coarseGridSolver_->solve(coarseGridSol);

	//Apply Interpolation Operator to coarse grid residual
	I_->apply_operator(fineGridResid, coarseGridResidual);

	//Add residual to initial approx
	vector_add(fineGridSol, fineGridResid, fineGridSol);

	//Final Relaxation Swwp
	fineGridSolver_->solve(fineGridSol);

	return true;
}