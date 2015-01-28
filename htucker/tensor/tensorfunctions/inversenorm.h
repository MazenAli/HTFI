#ifndef HTUCKER_TENSOR_TENSORFUNCTIONS_INVERSENORM_H
#define HTUCKER_TENSOR_TENSORFUNCTIONS_INVERSENORM_H 1

#include <htucker/dimensionindex/dimensionindex.h>

namespace htucker{

class inverseNorm{	
public:
	typedef double type;
	
	inverseNorm(const int _min, const int _max, const int _d);

	inverseNorm(const DimensionIndex &_min, const DimensionIndex &_max);

	inverseNorm(const inverseNorm & copy);

	double operator()(const DimensionIndex &vals) const;

	int dim() const;

	const DimensionIndex &
	getmin() const;

	const DimensionIndex &
	getmax() const;

	bool vecEval() const;

	void  
	vec(const DimensionIndex & vals, const int dim, flens::DenseVector<flens::Array<type> > & vec) const;

	inverseNorm & 
	operator=(const inverseNorm & copy);

private:
	int d;
	DimensionIndex minval,maxval;
};

#include<htucker/tensor/tensorfunctions/inversenorm.tcc>


}


#endif // HTUCKER_TENSOR_TENSORFUNCTIONS_INVERSENORM
