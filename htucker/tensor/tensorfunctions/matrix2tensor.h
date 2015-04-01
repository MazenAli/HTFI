#ifndef HTUCKER_TENSOR_TENSORFUNCTIONS_MATRIX2TENSOR_H
#define HTUCKER_TENSOR_TENSORFUNCTIONS_MATRIX2TENSOR_H 1

#include <htucker/dimensionindex/dimensionindex.h>
#include <htucker/indextools/m2vindexconverter.h>

namespace htucker{

class Matrix2Tensor{
	int d;
	DimensionIndex minval,maxval;
	bool matrixset;
	flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > * M;
	m2vindexconverter conv;

	public:
	typedef double type;
	
	Matrix2Tensor(const int _min, const int _max, const int _d);

	Matrix2Tensor(const DimensionIndex &_min, const DimensionIndex &_max);

	Matrix2Tensor(const Matrix2Tensor & copy); // hier auch noch den converter, dann copy...

	void 
	setMatrix(flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > &_M);

	void 
	setMatrix(flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > &_M, const DimensionIndex &_minidxrow, const DimensionIndex &_maxidxrow);

	double operator()(const DimensionIndex &vals) const;

	bool vecEval() const;

	int dim() const;

	const DimensionIndex &
	getmin() const;

	const DimensionIndex &
	getmax() const;

	void  
	vec(const DimensionIndex & vals, const int dim, flens::DenseVector<flens::Array<type> > & vec) const;

	flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > &
	getMatrix() const;

	bool 
	getmatrixset() const;

	const m2vindexconverter &
	getconverter() const;

	Matrix2Tensor &
	operator=(const Matrix2Tensor &copy); // den hier implementieren, dazu noch den converter kopieren....
};

#include<htucker/tensor/tensorfunctions/matrix2tensor.tcc>


}


#endif // HTUCKER_TENSOR_TENSORFUNCTIONS_MATRIX2TENSOR_H
