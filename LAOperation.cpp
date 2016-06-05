#include<LAOperation.h>

#include<cassert>
#include<cmath>
#include<cstddef>

double vector_dot_product(const std::vector<double>& vec1,
	const std::vector<double>& vec2)
{
	double result = 0.0;

	assert(vec1.size() == vec2.size());

	for (size_t i = 0; i < vec1.size(); ++i)
		result += vec1[i] * vec2[i];

	return result;
}

double vector_norm2(const std::vector<double>& vec1)
{
	double result = 0.0;

	for (size_t i = 0; i < vec1.size(); ++i)
		result += vec1[i] * vec1[i];

	result = std::sqrt(result);
	return result;
}

void vector_add(const std::vector<double>&& vec1,
	const std::vector<double>& vec2, std::vector<double>& result)
{
	assert(vec1.size() == vec2.size());

	if (result.size() != vec1.size())
		result.resize(vec1.size());

	for (size_t i = 0; i < vec1.size(); ++i)
		result[i] = vec1[i] + vec2[i];
}

void vector_subtract(const std::vector<double>&& vec1,
	const std::vector<double>& vec2, std::vector<double>& result)
{
	assert(vec1.size() == vec2.size());

	if (result.size() != vec1.size())
		result.resize(vec1.size());

	for (size_t i = 0; i < vec1.size(); ++i)
		result[i] = vec1[i] - vec2[i];
}

void vector_scalar_mult(double scalar, const std::vector<double>& vec1,
	std::vector<double>& res)
{
	for (size_t i = 0; i < vec1.size(); ++i)
		res[i] = scalar*vec1[i];
}

void vector_axpy(double scalar, const std::vector<double>& vecMult,
	const std::vector<double>& vecPlus, std::vector<double>& result)
{
	assert(vecMult.size() == vecPlus.size());

	if (result.size() != vecPlus.size())
		result.resize(vecPlus.size());

	for (size_t i = 0; i < vecMult.size(); ++i)
		result[i] = scalar*vecMult[i] + vecPlus[i];
}

void vector_axmy(double scalar, const std::vector<double>& vecMult,
	const std::vector<double>& vecMinus, std::vector<double>& result)
{
	assert(vecMult.size() == vecMinus.size());

	if (result.size() != vecMinus.size())
		result.resize(vecMinus.size());

	for (size_t i = 0; i < vecMult.size(); ++i)
		result[i] = scalar*vecMult[i] - vecMinus[i];
}