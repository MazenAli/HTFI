#ifndef HTUCKER_TENSOR_TENSORFUNCTIONS_INTEGERTENSOR_H
#define HTUCKER_TENSOR_TENSORFUNCTIONS_INTEGERTENSOR_H 1

#include <htucker/dimensionindex/dimensionindex.h>

namespace htucker{

class IntegerTensor{	
	int d;
	DimensionIndex minval,maxval;
	int offset;
public:
	
	typedef double type;
	
	IntegerTensor(const int _min, const int _max, const int _d);

	IntegerTensor(const DimensionIndex &_min, const DimensionIndex &_max);

	IntegerTensor(const IntegerTensor &copy);

	double operator()(const DimensionIndex &vals) const;

	void setoffset(const int _offset);

	bool vecEval() const;

	int getoffset() const;

	const DimensionIndex &
	getmax() const;

	const DimensionIndex & 
	getmin() const;

	int dim() const;

	void  
	vec(const DimensionIndex & vals, const int dim, flens::DenseVector<flens::Array<type> > & vec) const;

	IntegerTensor&
	operator=(const IntegerTensor &copy);

};

#include<htucker/tensor/tensorfunctions/integertensor.tcc>


}


#endif // HTUCKER_TENSOR_TENSORFUNCTIONS_INTEGERTENSOR_H
