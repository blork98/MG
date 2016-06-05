#include <SparseMatrix.h>

#include <assert.h>
#include <cmath>

SparseMatrix::SparseMatrix(size_t numRows, size_t numCols,  
	const std::vector<double>& val, 
	const std::vector<unsigned int>& colInd, 
	const std::vector<unsigned int>& rowPtr)
	: val_(val), colInd_(colInd), rowPtr_(rowPtr),
	numRows_(numRows), numCols_(numCols)
{
};

SparseMatrix& SparseMatrix::operator += ( double a )
{
	for ( unsigned int i = 0; i < val_.size(); ++i )
		val_[i] += a;

	return *this;
};

SparseMatrix& SparseMatrix::operator *= ( double a )
{
	for ( unsigned int i = 0; i < val_.size(); ++i )
		val_[i] *= a;

	return *this;
};

SparseMatrix& SparseMatrix::operator -= ( double a )
{
	for ( unsigned int i = 0; i < val_.size(); ++i )
		val_[i] -= a;

	return *this;
};

SparseMatrix& SparseMatrix::operator /= ( double a )
{
	for ( unsigned int i = 0; i < val_.size(); ++i )
		val_[i] /= a;

	return *this;
};

double SparseMatrix::vec_row_mult(size_t row, const std::vector<double>& in) const
{
	assert(numCols_ == in.size());

	double result = 0.0;

	//get range of elements of each row
	unsigned int startIndex = rowPtr_[row];
	unsigned int endIndex = rowPtr_[row + 1];

	for (unsigned int colCtr = startIndex; colCtr < endIndex; ++colCtr)
	{
		result += val_[colCtr] * in[colInd_[colCtr]];
	};

	return result;
}

void SparseMatrix::vec_multiply(const std::vector<double>& in, 
	std::vector<double>& out) const
{
	
	assert(numCols_ == in.size());

	if( in.size() != out.size() )
		out.resize(in.size());

	for( unsigned int rowCtr = 0; rowCtr < numRows_; ++rowCtr )
	{
		double result = 0.0;

		//get range of elements of each row
		unsigned int startIndex = rowPtr_[rowCtr];
		unsigned int endIndex = rowPtr_[rowCtr+1];

		for( unsigned int colCtr = startIndex; colCtr < endIndex; ++colCtr)
		{
			result += val_[colCtr]*in[colInd_[colCtr]];
		};

		out[rowCtr] = result;
	};

};

void SparseMatrix::get_diagonals ( std::vector<double>& diags ) const
{
	if( diags.size() != rows() )
		diags.resize(cols(),0.0);

	for( size_t rowCtr = 0; rowCtr < numRows_; ++rowCtr )
	{
		//get range of elements of each row
		unsigned int startIndex = rowPtr_[rowCtr];
		unsigned int endIndex = rowPtr_[rowCtr+1];

		diags[rowCtr] = 0.0;

		for( unsigned int colCtr = startIndex; colCtr < endIndex; ++colCtr)
		{
			if( colInd_[colCtr] == rowCtr )
			{
				diags[rowCtr] = val_[colCtr];
				break;
			}
		}
	}
}

size_t SparseMatrix::rows() const
{
	return numRows_;
}

size_t SparseMatrix::cols() const
{
	return numCols_;
}

size_t SparseMatrix::nonzero_elems() const
{
	return val_.size();
};

double SparseMatrix::entry_norm( int p ) const
{
	double result = 0.0;

	for( auto it = val_.begin(); it != val_.end(); ++it)
		result += std::pow(*it,p);

	result = std::pow(result,1/p);

	return result;
};