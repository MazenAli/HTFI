#ifndef HTUCKER_TENSOR_TENSORFUNCTIONS_VECTOR2TENSOR_H
#define HTUCKER_TENSOR_TENSORFUNCTIONS_VECTOR2TENSOR_H 1

#include <htucker/dimensionindex/dimensionindex.h>

namespace htucker{

class Vector2Tensor{
	int d;
	DimensionIndex minval,maxval;
	bool vectorset;
	flens::DenseVector<flens::Array<double > > * v;
	
	public:
	typedef double type;
	
	Vector2Tensor(const int _min, const int _max, const int _d);

	Vector2Tensor(const DimensionIndex &_min, const DimensionIndex &_max);

	Vector2Tensor(const Vector2Tensor & copy);

	void 
	setVector(flens::DenseVector<flens::Array<double > > &_v);

	double operator()(const DimensionIndex &vals) const;

	bool vecEval() const;

	bool VectorSet() const;

	int dim() const;

	void
	vec() const;

	const DimensionIndex &
	getmin() const;

	const DimensionIndex &
	getmax() const;

	const flens::DenseVector<flens::Array<double > > &
	getVector() const;
	
	Vector2Tensor&
	operator=(const Vector2Tensor &copy);
};

#include<htucker/tensor/tensorfunctions/vector2tensor.tcc>


}


#endif // HTUCKER_TENSOR_TENSORFUNCTIONS_VECTOR2TENSOR_H
