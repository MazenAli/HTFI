#ifndef HTUCKER_LINEARSOLVER_MAXSOLVER_H
#define HTUCKER_LINEARSOLVER_MAXSOLVER_H 1

namespace  htucker{

	
	template <typename T>
	HTuckerTree<T> positivepart(const HTuckerTree<T> & tree,const int maxrank, const T tol = std::numeric_limits<T>::epsilon(),const T newtontol = 0.000001, long maxIterations = LONG_MAX);
	
	template <typename T>
	HTuckerTree<T> absolutevalue(const HTuckerTree<T> & tree,const int maxrank, const T tol = std::numeric_limits<T>::epsilon(),const T newtontol = 0.000001, long maxIterations = LONG_MAX);

} //namespace htucker

#include <htucker/linearsolver/maxsolver.tcc>
#endif // HTUCKER_LINEARSOLVER_MAXSOLVER_H
