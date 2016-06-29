#ifndef FMG_H_
#define FMG_H_

#include <MGSolver.h>

class FMGSolver : public MultiGridSolver
{
public:
	FMGSolver(unsigned int numLevels, double tolerance,
		const std::vector<unsigned int>& sweeps,
		const std::vector<std::shared_ptr<SparseLinearSolver>>& solvers,
		const std::shared_ptr<InterpolationOperator>& I,
		const std::shared_ptr<RestrictionOperator>& R);

	virtual bool solve(bool interiorPointsOnly, std::vector<double>& sol) const;

private:
	void create_rhs() const;
	void v_cycle(int startLevel) const;
};

#endif

