#ifndef HTUCKER_TENSOR_TENSORFUNCTIONS_INVERSENORMAPPROXIMATION_H
#define HTUCKER_TENSOR_TENSORFUNCTIONS_INVERSENORMAPPROXIMATION_H 1

#include <htucker/dimensionindex/dimensionindex.h>

namespace htucker{

class inverseNormApproximation{
	int d;
	DimensionIndex minval,maxval;
	flens::DenseVector<flens::Array<double> > alpha,weights;

	public:
	
	typedef double type;

	inverseNormApproximation(const int _min, const int _max, const int _d);

	inverseNormApproximation(const DimensionIndex &_min, const DimensionIndex &_max);

	inverseNormApproximation(const inverseNormApproximation &copy);

	double operator()(const DimensionIndex &vals) const;

	DenseVectorList<double>
	getVectors() const;

	int dim() const;

	const DimensionIndex &
	getmin() const;

	const DimensionIndex &
	getmax() const;
	
	bool vecEval() const;

	void
	vec() const;

	inverseNormApproximation&
	operator=(const inverseNormApproximation & copy);
};

#include <htucker/tensor/tensorfunctions/inversenormapproximation.tcc>


} // namespace htucker


#endif // HTUCKER_TENSOR_TENSORFUNCTIONS_INVERSENORMAPPROXIMATION
