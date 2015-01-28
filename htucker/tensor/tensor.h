#ifndef HTUCKER_TENSOR_TENSOR_H
#define HTUCKER_TENSOR_TENSOR_H 1

//wrapps a discrete function into a tensor, i.e.
#include <htucker/dimensionindex/dimensionindex.h>

namespace htucker{

template <typename TensorFunction>
class tensor{
	int dimension;
	DimensionIndex minval,maxval;
	TensorFunction tf;

	public:
		tensor();

		tensor(const int _min, const int _max, const int d);

		tensor(const DimensionIndex &_min, const DimensionIndex &_max);

		tensor(const tensor<TensorFunction> &copy);

		typename TensorFunction::type
		LinfNorm(const int n) const;

		//Evaluation operator
		typename TensorFunction::type
        operator()(const DimensionIndex &vals) const;

		int dim() const;

		DimensionIndex & 
		getmax() const;

		DimensionIndex &
		getmin() const;

		TensorFunction &
		getTensorFunction() const;

		bool 
		vecEval() const;

		flens::DenseVector<flens::Array<typename TensorFunction::type> >
		vec(const DimensionIndex & vals, const int dim) const;

		tensor<TensorFunction> &
		operator=(const tensor<TensorFunction> &copy);
};

#include <htucker/tensor/tensor.tcc>

} // namespace htucker


#include <htucker/tensor/tensorfunctions/tensorfunctions.h>
#endif // HTUCKER_TENSOR_TENSOR_H 
