#ifndef GRID_OPERATORS_H_
#define GRID_OPERATORS_H_

#include <vector>

class RestrictionOperator
{
public:
	virtual void apply_operator(std::vector<double>& fineGrid, 
		std::vector<double>& coarseGrid, bool interiorPointsOnly) const = 0;
};

class Injecction : public RestrictionOperator
{
public:
	void apply_operator(std::vector<double>& fineGrid,
		std::vector<double>& coarseGrid, bool interiorPointsOnly) const;
};

class FullWeighing : public RestrictionOperator
{
public:
	void apply_operator(std::vector<double>& fineGrid,
		std::vector<double>& coarseGrid, bool interiorPointsOnly) const;
};

#endif
