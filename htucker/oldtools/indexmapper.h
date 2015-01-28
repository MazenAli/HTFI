#ifndef LAWA_HTUCKER_TENSOR_H
#define LAWA_HTUCKER_TENSOR_H 1

//wrapps a discrete function into a tensor, i.e.

#include <fima/HTucker/dimensionindex.h>

template <typename T>
class tensor{
    public:
	typedef T (*funpoint)(DimensionIndex);

	funpoint f;
	int dimension;
	DimensionIndex minval,maxval;

	//constructors
	tensor():dimension(0){};

	tensor(funpoint g, int d):f(g),dimension(d){};

	tensor(funpoint g, int d, int _min, int _max);

	tensor(funpoint g, int d, DimensionIndex _min, DimensionIndex _max);

	T
	LinfNorm(int n);

	//Evaluation operator
	T operator()(DimensionIndex para);
};

#include <fima/HTucker/tensor.tcc>

#endif // LAWA_HTUCKER_TENSOR_H