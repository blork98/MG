#ifndef LAOPERATION_H_
#define LAOPERATION_H_

#include<vector>

double vector_dot_product(const std::vector<double>& vec1,
	const std::vector<double>& vec2);
double vector_norm_2(const std::vector<double>& vec1);

void vector_add(const std::vector<double>& vec1,
	const std::vector<double>& vec2, std::vector<double>& result);
void vector_subtract(const std::vector<double>& vec1,
	const std::vector<double>& vec2, std::vector<double>& result);
void vector_scalar_mult(double scalar, const std::vector<double>& vec1,
	std::vector<double>& res);
void vector_axpy(double scalar, const std::vector<double>& vecMult,
	const std::vector<double>& vecPlus, std::vector<double>& result);
void vector_axmy(double scalar, const std::vector<double>& vecMult,
	const std::vector<double>& vecMinus, std::vector<double>& result);
void vector_xmay(double scalar, const std::vector<double>& vecMult,
	const std::vector<double>& vecMinus, std::vector<double>& result);

#endif