#ifndef MGSOLVER_H_
#define MGSOLVER_H_

#include<memory>
#include<vector>

#include<SparseMatrix.h>
#include<LinearIterativeSolver.h>

class MultiGridSolver
{
public:
	MultiGridSolver(unsigned int numLevels, double tolerance, 
		const std::vector<unsigned int>& sweeps, 
		const std::vector<std::shared_ptr<SparseLinearSolver>>& solvers);

	virtual bool solve( bool interiorPointsOnly, std::vector<double>& sol ) const = 0;

	unsigned int num_levels() const;
	void set_rhs(const std::vector<double>& rhs);
	void set_ith_A(size_t i, const std::shared_ptr<SparseMatrix>& A);
	void set_A(const std::vector<std::shared_ptr<SparseMatrix>>& A);
	void set_solvers(const std::vector<std::shared_ptr<SparseLinearSolver>>& solvers_);

protected:
	unsigned int numLevels_;
	double tolerance_;
	std::vector<std::shared_ptr<SparseMatrix>> A_;
	std::vector<std::vector<double>> rhs_, resids_;
	std::vector<std::shared_ptr<SparseLinearSolver>> solvers_;
	std::vector<unsigned int> sweeps_;

	void resize_rhs_resids();
	void verify_inputs() const;
};

#endif
