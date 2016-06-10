#ifndef GRID_OPERATORS_H_
#define GRID_OPERATORS_H_

#include <vector>

class RestrictionOperator
{
public:
	virtual void apply_operator(std::vector<double>& fineGrid, 
		std::vector<double>& coarseGrid) const = 0;
	virtual unsigned int dim() const;
};

class InterpolationOperator
{
public:
	virtual void apply_operator(std::vector<double>& fineGrid, 
		std::vector<double>& coarseGrid) const = 0;
	virtual unsigned int dim() const;
};

class Injecction1D : public RestrictionOperator
{
public:
	void apply_operator(std::vector<double>& fineGrid,
		std::vector<double>& coarseGrid) const;
	virtual unsigned int dim() const;
};

class FullWeighing1D : public RestrictionOperator
{
public:
	void apply_operator(std::vector<double>& fineGrid,
		std::vector<double>& coarseGrid) const;
	virtual unsigned int dim() const;
};

class LinearInterpolation1D : public InterpolationOperator
{
public:
	virtual void apply_operator(std::vector<double>& fineGrid,
		std::vector<double>& coarseGrid) const;
	virtual unsigned int dim() const;
};

#endif
