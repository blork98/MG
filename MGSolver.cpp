#include <MGSolver.h>

#include <assert.h>

using std::shared_ptr;

MultiGridSolver::MultiGridSolver(unsigned int numLevels, double tolerance,
	const std::vector<unsigned int>& sweeps,
	const std::vector<shared_ptr<SparseLinearSolver>>& solvers)
	: numLevels_(numLevels), tolerance_(tolerance), sweeps_(sweeps), solvers_(solvers),
	A_(0, nullptr), rhs_(numLevels_), resids_(numLevels_)
{
}

unsigned int MultiGridSolver::num_levels() const
{
	return numLevels_;
}

void MultiGridSolver::set_rhs(const std::vector<double>& rhs)
{
	rhs_[0] = rhs;
}

void MultiGridSolver::set_ith_A(size_t i, const std::shared_ptr<SparseMatrix>& A)
{
	if (i < A_.size())
		A_[i] = A;
}

void MultiGridSolver::set_A(const std::vector<std::shared_ptr<SparseMatrix>>& A)
{
	A_ = A;
	numLevels_ = A_.size();

	resize_rhs_resids();
}

void  MultiGridSolver::set_solvers(const std::vector<std::shared_ptr<SparseLinearSolver>>& solvers)
{
	solvers_ = solvers;
}

void MultiGridSolver::resize_rhs_resids()
{
	for (size_t i = 1; i < A_.size(); ++i)
	{
		rhs_[i] = std::move(std::vector<double>(A_[i]->cols(),0.0));
		resids_[i] = std::move(std::vector<double>(A_[i]->cols(),0.0));
	}
	resids_[0] = std::move(std::vector<double>(A_[0]->cols(), 0.0));
}

void  MultiGridSolver::verify_inputs() const
{
	//check level numbers
	assert(numLevels_ == A_.size());
	assert(A_.size() == rhs_.size());
	assert(rhs_.size() == resids_.size());
	assert(solvers_.size() == numLevels_);
	assert(sweeps_.size() == numLevels_);

	assert(tolerance_ >= 0.0);
	for (auto it = sweeps_.begin(); it != sweeps_.end(); ++it)
		assert(*it > 0);

	for (auto it = solvers_.begin(); it != solvers_.end(); ++it)
		assert(*it != nullptr);

	for (unsigned int i = 0; i < numLevels_; ++i)
	{
		assert(A_[i]->cols() == rhs_[i].size());
		assert(A_[i]->cols() == resids_[i].size());
	}
}
