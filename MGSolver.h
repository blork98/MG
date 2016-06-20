#ifndef MGSOLVER_H_
#define MGSOLVER_H_

#include<memory>
#include<vector>

#include<SparseMatrix.h>
#include<LinearIterativeSolver.h>
#include<GridOperators.h>

class MultiGridSolver
{
public:
	MultiGridSolver(unsigned int numLevels, double tolerance, 
		const std::vector<unsigned int>& sweeps, 
		const std::vector<std::shared_ptr<SparseLinearSolver>>& solvers,
		const std::shared_ptr<InterpolationOperator>& I,
		const std::shared_ptr<RestrictionOperator>& R);

	virtual bool solve( bool interiorPointsOnly, std::vector<double>& sol ) const = 0;
	unsigned int num_levels() const;
	void set_rhs(const std::shared_ptr<std::vector<double>>& rhs);
	void set_ith_A(size_t i, const std::shared_ptr<SparseMatrix>& A);
	void set_A(const std::vector<std::shared_ptr<SparseMatrix>>& A);
	void set_solvers(const std::vector<std::shared_ptr<SparseLinearSolver>>& solvers_);
	void set_interpolation(const std::shared_ptr<InterpolationOperator>& I);
	void set_restriction(const std::shared_ptr<RestrictionOperator>& R);
	std::shared_ptr<SparseLinearSolver>& solver(size_t i);

protected:
	unsigned int numLevels_;
	double tolerance_;
	std::vector<std::shared_ptr<SparseMatrix>> A_;
	mutable std::vector<std::vector<double>> resids_, sols_;
	std::vector<std::shared_ptr<std::vector<double>>> rhs_;
	std::vector<std::shared_ptr<SparseLinearSolver>> solvers_;
	std::shared_ptr<RestrictionOperator> R_;
	std::shared_ptr<InterpolationOperator> I_;

	void resize_rhs_resids();
	void verify_inputs() const;
};

#endif
