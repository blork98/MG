#include <MGSolver.h>

#include <assert.h>

using std::shared_ptr;

MultiGridSolver::MultiGridSolver(unsigned int numLevels, double tolerance,
	const std::vector<unsigned int>& sweeps,
	const std::vector<shared_ptr<SparseLinearSolver>>& solvers,
	const std::shared_ptr<InterpolationOperator>& I,
	const std::shared_ptr<RestrictionOperator>& R)
	: numLevels_(numLevels), tolerance_(tolerance), solvers_(solvers),
	A_(0, nullptr), rhs_(0,nullptr), resids_(numLevels_),
	sols_(numLevels_), R_(R), I_(I)
{
	//resize_rhs_resids();
}

unsigned int MultiGridSolver::num_levels() const
{
	return numLevels_;
}

void MultiGridSolver::set_rhs(const std::shared_ptr<std::vector<double>>& rhs)
{
	rhs_[0] = rhs;

	if ( rhs->size() !=  rhs_[0]->size() )
		resize_rhs_resids();
}

void MultiGridSolver::set_ith_A(size_t i, const std::shared_ptr<SparseMatrix>& A)
{
	if (i < A_.size())
		A_[i] = A;
}

void MultiGridSolver::set_A(const std::vector<std::shared_ptr<SparseMatrix>>& A)
{	
	if (A.size() != A_.size() )
	{
		A_ = A;
		rhs_ = std::vector<std::shared_ptr<std::vector<double>>>(A_.size());
		resids_ = std::vector<std::vector<double>>(A.size());
		sols_ = std::vector<std::vector<double>>(A.size());
		resize_rhs_resids();
		numLevels_ = A_.size();
	} else {
		A_ = A;
	}

}

void  MultiGridSolver::set_solvers(const std::vector<std::shared_ptr<SparseLinearSolver>>& solvers)
{
	solvers_ = solvers;
}

void MultiGridSolver::resize_rhs_resids()
{
	rhs_.resize(numLevels_);
	resids_.resize(numLevels_);
	sols_.resize(numLevels_);

	for (size_t i = 1; i < A_.size(); ++i)
	{
		rhs_[i] = std::make_shared<std::vector<double>>(std::vector<double>(A_[i]->cols(),0.0));
		resids_[i] = std::move(std::vector<double>(A_[i]->cols(),0.0));
		sols_[i] = std::move(std::vector<double>(A_[i]->cols(), 0.0));
	}
	resids_[0] = std::move(std::vector<double>(A_[0]->cols(), 0.0));
	sols_[0] = std::move(std::vector<double>(A_[0]->cols(), 0.0));
}

void MultiGridSolver::set_interpolation(const std::shared_ptr<InterpolationOperator>& I)
{
	I_ = I;
}

void MultiGridSolver::set_restriction(const std::shared_ptr<RestrictionOperator>& R)
{
	R_ = R;
}

void  MultiGridSolver::verify_inputs() const
{
	//check level numbers
	assert(numLevels_ == A_.size());
	assert(A_.size() == rhs_.size());
	assert(rhs_.size() == resids_.size());
	assert(solvers_.size() == numLevels_);

	assert(I_ != nullptr);
	assert(R_ != nullptr);

	assert(tolerance_ >= 0.0);

	for (auto it = solvers_.begin(); it != solvers_.end(); ++it)
		assert(*it != nullptr);

	for (unsigned int i = 0; i < numLevels_; ++i)
	{
		assert(A_[i]->cols() == rhs_[i]->size());
		assert(A_[i]->cols() == resids_[i].size());
		assert(A_[i]->cols() == sols_[i].size());
	}

	for (size_t i = 1; i < numLevels_; ++i)
	{
		assert(rhs_[i]->size() == (rhs_[i-1]->size()/2 + 1) );
		assert(sols_[i].size() == (sols_[i - 1].size() / 2 + 1));
		assert(resids_[i].size() == (resids_[i - 1].size() / 2 + 1));
	}
}

std::shared_ptr<SparseLinearSolver>&  MultiGridSolver::solver(size_t i)
{
	assert(i < solvers_.size());
	return solvers_[i];
}

