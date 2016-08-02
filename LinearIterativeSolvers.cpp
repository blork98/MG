#include <LinearIterativeSolver.h>
#include <LAOperation.h>

#include <algorithm>
#include<assert.h>

SparseIterativeLinearSolver::SparseIterativeLinearSolver(int maxIters, double tolerance)
	:maxIters_(maxIters), tolerance_(tolerance),
	A_(nullptr), rhs_(nullptr), numIters_(0)
{
}

int SparseIterativeLinearSolver::num_iters() const
{
	return numIters_;
}

void SparseIterativeLinearSolver::set_max_iters(int maxIters)
{
	maxIters_ = maxIters;
}

void SparseIterativeLinearSolver::set_A(const std::shared_ptr<SparseMatrix>& A)
{
	A_ = A;
}

void SparseIterativeLinearSolver::set_rhs(const std::shared_ptr<std::vector<double>>& b)
{
	rhs_ = b;
}

bool SparseIterativeLinearSolver::convergence_achieved(const std::vector<double>& vec1,
	const std::vector<double>& vec2) const
{
	double result = 0.0;

	for (size_t elem = 0; elem < vec1.size(); ++elem)
	{
		result += std::pow((vec1[elem] - vec2[elem]), 2);
	}

	if (std::sqrt(result) < tolerance_)
		return true;
	else
		return false;
}

JacobiSolver::JacobiSolver( int maxIters, double tolerance )
	: SparseIterativeLinearSolver(maxIters,tolerance)
{
}

bool JacobiSolver::solve( std::vector<double>& sol )
{

	assert(A_ != nullptr);
	assert(rhs_ != nullptr);

	if( sol.size() != rhs_->size() )
		sol.resize(rhs_->size());

	if( A_->cols() != rhs_->size() )
		return false;

	double result = 0.0;
	std::vector<double> tempSol(sol.size(),0.0);
	std::vector<double> *solCurr, *solPrev;
	solPrev = &sol;
	solCurr = &tempSol;

	std::vector<double> diagElems(A_->rows(), 0.0);
	A_->get_diagonals(diagElems);

	for (int iter = 1; iter <= maxIters_; ++iter)
	{
		for (size_t i = 0; i < A_->rows(); ++i)
		{
			(*solCurr)[i] = ((*rhs_)[i] -
				(A_->vec_row_mult(i, *solPrev) - diagElems[i] * (*solPrev)[i]))
				/ diagElems[i];
		}

		//check for convergence
		if (convergence_achieved(*solPrev, *solCurr))
		{
			if ((iter % 2) != 0)
				sol = tempSol;
			
			numIters_ = iter;

			return true;
		}

		std::swap(solCurr, solPrev);
	}

	return false;
}

GaussSeidelSolver::GaussSeidelSolver(int maxIters, double tolerance)
	: SparseIterativeLinearSolver(maxIters, tolerance)
{
}

bool GaussSeidelSolver::solve(std::vector<double>& sol)
{
	assert(A_ != nullptr);
	assert(rhs_ != nullptr);

	if (sol.size() != rhs_->size())
		sol.resize(rhs_->size());

	if (A_->cols() != rhs_->size())
		return false;

	std::vector<double> diagElems(A_->rows(), 0.0);
	A_->get_diagonals(diagElems);

	std::vector<double> prevSol = sol;

	for (int iter = 1; iter <= maxIters_; ++iter)
	{

		for (size_t i = 0; i < A_->rows(); ++i)
		{
			sol[i] = ((*rhs_)[i] -
				(A_->vec_row_mult(i, sol) - diagElems[i] * sol[i]))
				/ diagElems[i];
		}

		//check for convergence
		if (convergence_achieved(sol, prevSol))
		{
			numIters_ = iter;
			return true;
		}

		prevSol = sol;

	}

	return false;
}

SORSolver::SORSolver(int maxIters, double tolerance, double w, bool isSymmetric)
	: SparseIterativeLinearSolver(maxIters, tolerance),
	w_(w), isSymmetric_(isSymmetric)
{
}

void SORSolver::set_w(double w)
{
	w_ = w;
}

void SORSolver::set_symmetric(bool isSymmetric)
{
	isSymmetric_ = isSymmetric;
}

bool SORSolver::solve(std::vector<double>& sol)
{
	assert(A_ != nullptr);
	assert(rhs_ != nullptr);

	if (sol.size() != rhs_->size())
		sol.resize(rhs_->size());

	if (A_->cols() != rhs_->size())
		return false;

	std::vector<double> diagElems(A_->rows(), 0.0);
	A_->get_diagonals(diagElems);

	std::vector<double> prevSol = sol;
	std::vector<double> solIntermed = sol;

	for (int iter = 1; iter <= maxIters_; ++iter)
	{
		//forward SOR
		for (size_t i = 0; i < A_->rows(); ++i)
		{
			sol[i] = ((*rhs_)[i] -
				(A_->vec_row_mult(i, sol) - diagElems[i] * sol[i]))
				/ diagElems[i];
			sol[i] = w_*sol[i] + (1 - w_)*prevSol[i];
		}

		//backward SOR
		if ( isSymmetric_ )
		{
			solIntermed = sol;
			for (size_t i = 0; i < A_->rows(); ++i)
			{
				size_t ctr = A_->rows() - 1 - i;
				//In Templates this is only A*x not b - A*x
				sol[ctr] = ((*rhs_)[ctr] -
					(A_->vec_row_mult(ctr, sol) - diagElems[ctr] * sol[ctr]))
					/ diagElems[i]; 
				sol[ctr] = w_*sol[ctr] + (1 - w_)*solIntermed[ctr];
			}
		}
		
		//check for convergence
		if (convergence_achieved(sol, prevSol))
		{
			numIters_ = iter;
			return true;
		}

		prevSol = sol;
	}

	return false;
}

WeightedJacobiSolver::WeightedJacobiSolver(int maxIters, double tolerance, double w)
	: SparseIterativeLinearSolver(maxIters, tolerance), w_(w)
{
}

bool WeightedJacobiSolver::solve(std::vector<double>& sol)
{
	assert(A_ != nullptr);
	assert(rhs_ != nullptr);

	if (sol.size() != rhs_->size())
		sol.resize(rhs_->size());

	if (A_->cols() != rhs_->size())
		return false;

	double result = 0.0;
	std::vector<double> tempSol(sol.size(), 0.0);
	std::vector<double> *solCurr, *solPrev;
	solPrev = &sol;
	solCurr = &tempSol;

	std::vector<double> diagElems(A_->rows(), 0.0);
	A_->get_diagonals(diagElems);

	for (int iter = 1; iter <= maxIters_; ++iter)
	{
		for (size_t i = 0; i < A_->rows(); ++i)
		{
			(*solCurr)[i] = ((*rhs_)[i] -
				(A_->vec_row_mult(i, *solPrev) - diagElems[i] * (*solPrev)[i]))
				/ diagElems[i];
			(*solCurr)[i] = w_*(*solCurr)[i] + (1 - w_)*((*solPrev)[i]);
		}

		//check for convergence
		if (convergence_achieved(*solPrev, *solCurr))
		{
			if ((iter % 2) != 0)
				sol = tempSol;

			numIters_ = iter;

			return true;
		}

		std::swap(solCurr, solPrev);
	}

	return false;
}

RedBlackGaussSeidelSolver::RedBlackGaussSeidelSolver(int maxIters, double tolerance)
	:SparseIterativeLinearSolver(maxIters, tolerance)
{
}

bool RedBlackGaussSeidelSolver::solve(std::vector<double>& sol)
{
	assert(A_ != nullptr);
	assert(rhs_ != nullptr);

	if (sol.size() != rhs_->size())
		sol.resize(rhs_->size());

	if (A_->cols() != rhs_->size())
		return false;

	std::vector<double> diagElems(A_->rows(), 0.0);
	A_->get_diagonals(diagElems);

	std::vector<double> prevSol = sol;

	for (int iter = 1; iter <= maxIters_; ++iter)
	{

		//update Red Nodes (Even)
		for (size_t i = 0; i < A_->rows(); i+=2)
		{
			sol[i] = ((*rhs_)[i] -
				(A_->vec_row_mult(i, sol) - diagElems[i] * sol[i]))
				/ diagElems[i];
		}

		//update Black Nodes (Odd)
		for (size_t i = 1; i < A_->rows(); i += 2)
		{
			sol[i] = ((*rhs_)[i] -
				(A_->vec_row_mult(i, sol) - diagElems[i] * sol[i]))
				/ diagElems[i];
		}

		//check for convergence
		if (convergence_achieved(sol, prevSol))
		{
			numIters_ = iter;
			return true;
		}

		prevSol = sol;

	}

	return false;
}