#ifndef MGTWOGRIDCORRECTION_H_
#define MGTWOGIRDCORRECTION_H_

#include<memory>
#include<vector>

#include<SparseMatrix.h>
#include<LinearIterativeSolver.h>
#include<GridOperators.h>

class TwoGridCorrection 
{
public:
	TwoGridCorrection(unsigned int numRelaxSweeps,
		const std::shared_ptr<SparseLinearSolver>& coarseGridSolver,
		const std::shared_ptr<SparseIterativeLinearSolver> fineGridSolver);

	bool solve(bool interiorPointsOnly, std::vector<double>& sol) const;

	void set_coarse_grid_solver(const std::shared_ptr<SparseLinearSolver>& coarseGridSolver);
	void set_fine_grid_solver(const std::shared_ptr<SparseIterativeLinearSolver>& fineGridSolver);
	void set_A_fine(const std::shared_ptr<SparseMatrix>& A);
	void set_A_coarse(const std::shared_ptr<SparseMatrix>& A);
	void set_rhs(const std::shared_ptr<std::vector<double>>& rhs);
	void set_restriction(const std::shared_ptr<RestrictionOperator>& R);
	void set_interpolation(const std::shared_ptr<InterpolationOperator>& I);

	const std::shared_ptr<SparseMatrix>& A_fine() const;
	const std::shared_ptr<SparseMatrix>& A_coarse() const;
	const std::shared_ptr<std::vector<double>>& rhs() const;

private:
	unsigned int numRelaxSweeps_;
	mutable std::shared_ptr<SparseLinearSolver> coarseGridSolver_;
	mutable std::shared_ptr<SparseIterativeLinearSolver> fineGridSolver_;
	std::shared_ptr<SparseMatrix> fineA_, coarseA_;
	std::shared_ptr<std::vector<double>> rhs_;
	std::shared_ptr<RestrictionOperator> R_;
	std::shared_ptr<InterpolationOperator> I_;

	void verify_inputs() const;
};

#endif