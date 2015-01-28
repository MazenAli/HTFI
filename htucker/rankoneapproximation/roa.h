#ifndef HTUCKER_RANKONEAPPROXIMATION_ROA_H
#define HTUCKER_RANKONEAPPROXIMATION_ROA_H 1

//Class for a sucessive rank one approximation

#include <htucker/dimensionindex/dimensionindex.h>
#include <htucker/tensor/tensor.h>
#include <iostream>


namespace htucker{




template <typename T, typename TensorFunction>
class ROA{
	private:
		int d;
		
		DimensionIndex t,tbar,tcomplement;
	public:
		
		tensor<TensorFunction> A;
		flens::GeMatrix<flens::FullStorage<T,flens::ColMajor> > Atensor,U,V;
		flens::DenseVector<flens::Array<T> > s,inverse_s;
		DimensionIndexList pivots;
		double SVD_epsilon;
	
	//constructors

		ROA(const double eps = 0.000000000001);

		ROA(const tensor<TensorFunction> &_A, const DimensionIndex &_t, const DimensionIndex &_tbar, const int _dim, const double _eps = 0.000000000001 );

		ROA(const ROA<T,TensorFunction> &copy);

		DimensionIndex 
		GreedyPivotSearch(const int lmax, const DimensionIndexList & fatherpivots, double &error) const;

		flens::DenseVector<flens::Array<T> >
		evaluate(const DimensionIndex &initialIndex, const DimensionIndex &_t, const  int dim, const int _min, const int _max) const;

		flens::DenseVector<flens::Array<T> >
		evaluate(const DimensionIndex &initialIndex, const DimensionIndex & _t, const DimensionIndexList &list, const DimensionIndex & list_dims) const;

		//approximates the given tensor with a rank one approximation of rank "rank"
		void approximate(const int rank, const DimensionIndexList &fatherpivots, const int l);
	
		//approximate the given tensor with a rank one approximation in blackbox-fashion (until error is smaller than a given boundary
		void approximate(const double epsilon, const DimensionIndexList &fatherpivots, const int l);

		// adds another pivot element to the ROA
		void addpivot(const DimensionIndex & pivot);

		DimensionIndex 
		getTbar() const;

		DimensionIndex
		getT() const;

		int dim() const;

		void print() const;

		// Evaluation operator
		T 
		operator()(const DimensionIndex &vals) const;

		ROA<T,TensorFunction> &
		operator=(const ROA<T,TensorFunction> &rhs);

};

} // namespace htucker


#include <htucker/rankoneapproximation/roa.tcc>

#endif // HTUCKER_RANKONEAPPROXIMATION_ROA_H
