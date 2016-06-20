#ifndef MGVCYCLE_H_
#define MGVCYCLE_H_

#include <MGSolver.h>

class MGSolverVCycle : public MultiGridSolver
{
public:
	MGSolverVCycle(unsigned int numLevels, double tolerance,
		const std::vector<unsigned int>& sweeps,
		const std::vector<std::shared_ptr<SparseLinearSolver>>& solvers,
		const std::shared_ptr<InterpolationOperator>& I,
		const std::shared_ptr<RestrictionOperator>& R);

	virtual bool solve(bool interiorPointsOnly, std::vector<double>& sol) const;
};

#endif
