#ifndef HTUCKER_TENSOR_TENSORFUNCTIONS_EXPONENTIALSUM_H
#define HTUCKER_TENSOR_TENSORFUNCTIONS_EXPONENTIALSUM_H 1

#include <htucker/dimensionindex/dimensionindex.h>

namespace htucker{


class exponentialSum{
	int d;
	DimensionIndex minval,maxval; 

public:
	typedef double type;

	exponentialSum(const int _min, const int _max, const int _d); 

	exponentialSum(const DimensionIndex &_min, const DimensionIndex &_max);

	exponentialSum(const exponentialSum &copy);

	double operator() (const DimensionIndex &vals) const;

	int dim() const;

	const DimensionIndex & 
	getmin() const;

	const DimensionIndex &
	getmax() const;
	
	bool vecEval() const;

	void
	vec() const;


	exponentialSum&
	operator= (const exponentialSum & copy);
};


#include<htucker/tensor/tensorfunctions/exponentialsum.tcc>
} // htucker


#endif // HTUCKER_TENSOR_TENSORFUNCTIONS_EXPONENTIALSUM
