#ifndef SPARSEMATRIX_H_
#define SPARSEMATRIX_H_

#include <vector>

class SparseMatrix
{
public:
	SparseMatrix( size_t numRows, size_t numCols, 
		const std::vector<double>& val, 
		const std::vector<unsigned int>& colInd, 
		const std::vector<unsigned int>& rowPtr);

	SparseMatrix& operator += ( double a );
	SparseMatrix& operator *= ( double a );
	SparseMatrix& operator -= ( double a );
	SparseMatrix& operator /= ( double a );

	void vec_multiply(const std::vector<double>& in, std::vector<double>& out) const;
	void get_diagonals ( std::vector<double>& in ) const;

	double vec_row_mult( size_t row, const std::vector<double>& in) const;

	size_t rows() const;
	size_t cols() const;
	size_t nonzero_elems() const;

	double entry_norm( int p ) const;

private:
	std::vector<double> val_; 
	std::vector<unsigned int> colInd_, rowPtr_;
	size_t numRows_, numCols_;
};

#endif