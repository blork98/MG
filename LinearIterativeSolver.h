#ifndef LINEAR_ITERATIVE_SOLVER_H_
#define LINEAR_ITERATIVE_SOLVER_H_

#include <memory>
#include <vector>

#include <SparseMatrix.h>

class SparseLinearSolver
{
public:
	virtual bool solve( std::vector<double>& sol ) const = 0;
	virtual void set_A(const std::shared_ptr<SparseMatrix>& A) = 0;
	virtual void set_rhs(const std::shared_ptr<std::vector<double>>& b) = 0;
};

class SparseIterativeLinearSolver : public SparseLinearSolver
{
public:
	SparseIterativeLinearSolver(int maxIters, double tolerance);

	virtual bool solve(std::vector<double>& sol) const = 0;
	virtual void set_A(const std::shared_ptr<SparseMatrix>& A);
	virtual void set_rhs(const std::shared_ptr<std::vector<double>>& b);

protected:
	int maxIters_;
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
	bool solve( std::vector<double>& sol ) const;
};

class GaussSeidelSolver : public SparseIterativeLinearSolver
{
public:
	GaussSeidelSolver(int maxIters, double tolerance);
	bool solve(std::vector<double>& sol) const;
};

class SORSolver : public SparseIterativeLinearSolver
{
public:
	SORSolver(int maxIters, double tolerance, double w, bool isSymmetric);
	bool solve(std::vector<double>& sol) const;

	void set_w(double w);
	void set_symmetric(bool isSymmetric);

private:
	bool isSymmetric_;
	double w_;
};

#endif