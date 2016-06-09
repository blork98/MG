#include <GridOperators.h>

void Injecction::apply_operator(std::vector<double>& fineGrid,
	std::vector<double>& coarseGrid, bool interiorPointsOnly) const
{
	size_t h = fineGrid.size()-1;
	size_t h2 = coarseGrid.size()-1;

	if (h2 != h / 2)
		coarseGrid.resize( (h/2) , 0.0);

	for (size_t i = 0; i < coarseGrid.size(); ++i)
		coarseGrid[i] = fineGrid[2 * i];

}

void FullWeighing::apply_operator(std::vector<double>& fineGrid,
	std::vector<double>& coarseGrid, bool interiorPointsOnly) const
{
	size_t h = fineGrid.size() - 1;
	size_t h2 = coarseGrid.size() - 1;

	if (h2 != h / 2)
		coarseGrid.resize((h / 2), 0.0);

	//copy endpoints
	coarseGrid[0] = fineGrid[0];
	*(coarseGrid.end()-1) = *(fineGrid.end()-1);
	for (size_t i = 1; i < coarseGrid.size(); ++i)
		coarseGrid[i] = 0.25*(fineGrid[2*i-1] + 2*fineGrid[2*i] + fineGrid[2*i+1]);
}