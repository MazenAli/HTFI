#ifndef HTUCKER_TENSOR_TENSORFUNCTIONS_EXPONENTIALSUMAPPROXIMATION_H
#define HTUCKER_TENSOR_TENSORFUNCTIONS_EXPONENTIALSUMAPPROXIMATION_H 1

#include <htucker/dimensionindex/dimensionindex.h>

namespace htucker{


class exponentialSumApproximation{
	int d;
	DimensionIndex minval,maxval; 
	flens::DenseVector<flens::Array<double> > alpha,weights;

public:
	typedef double type;

	exponentialSumApproximation(const int _min, const int _max, const int _d); 

	exponentialSumApproximation(const DimensionIndex &_min, const DimensionIndex &_max);

	exponentialSumApproximation(const exponentialSumApproximation &copy);

	double operator() (const DimensionIndex &vals) const;

	int dim() const;

	const DimensionIndex & 
	getmin() const;

	const DimensionIndex &
	getmax() const;


	DenseVectorList<double> 
	getVectors() const;

	bool vecEval() const;

	void  
	vec(const DimensionIndex & vals, const int dim, flens::DenseVector<flens::Array<type> > & vec) const;



	exponentialSumApproximation &
	operator=(const exponentialSumApproximation & copy);
};


#include<htucker/tensor/tensorfunctions/exponentialsumapproximation.tcc>
} // htucker


#endif // HTUCKER_TENSOR_TENSORFUNCTIONS_EXPONENTIALSUMAPPROXIMATION
