#ifndef LINEAR_ITERATIVE_SOLVER_H_
#define LINEAR_ITERATIVE_SOLVER_H_

#include <memory>
#include <vector>

#include <SparseMatrix.h>

class SparseLinearSolver
{
public:
	virtual bool solve( std::vector<double>& sol ) = 0;
	virtual void set_A(const std::shared_ptr<SparseMatrix>& A) = 0;
	virtual void set_rhs(const std::shared_ptr<std::vector<double>>& b) = 0; //change this to pass by reference
};

class SparseIterativeLinearSolver : public SparseLinearSolver
{
public:
	SparseIterativeLinearSolver(int maxIters, double tolerance);

	virtual bool solve(std::vector<double>& sol) = 0;
	virtual void set_A(const std::shared_ptr<SparseMatrix>& A);
	virtual void set_rhs(const std::shared_ptr<std::vector<double>>& b);

	void set_max_iters(int maxIters);
	int num_iters() const;

protected:
	int maxIters_, numIters_;
	double tolerance_;
	std::shared_ptr<SparseMatrix> A_;
	std::shared_ptr<std::vector<double>> rhs_;

	virtual bool convergence_achieved(const std::vector<double>& vec1,
		const std::vector<double>& vec2) const;
};

class JacobiSolver : public SparseIterativeLinearSolver
{
public:
	JacobiSolver( int maxIters, double tolerance );
	bool solve( std::vector<double>& sol );
};

class WeightedJacobiSolver : public SparseIterativeLinearSolver
{
public:
	WeightedJacobiSolver(int maxIters, double tolerance, double w);
	bool solve(std::vector<double>& sol);

private:
	double w_;
};

class GaussSeidelSolver : public SparseIterativeLinearSolver
{
public:
	GaussSeidelSolver(int maxIters, double tolerance);
	bool solve(std::vector<double>& sol);
};

class RedBlackGaussSeidelSolver : public SparseIterativeLinearSolver
{
public:
	RedBlackGaussSeidelSolver(int maxIters, double tolerance);
	bool solve(std::vector<double>& sol);
};

class SORSolver : public SparseIterativeLinearSolver
{
public:
	SORSolver(int maxIters, double tolerance, double w, bool isSymmetric);
	bool solve(std::vector<double>& sol);

	void set_w(double w);
	void set_symmetric(bool isSymmetric);

private:
	bool isSymmetric_;
	double w_;
};

#endif