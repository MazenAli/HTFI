#ifndef HTUCKER_LINEARSOLVER_LINEARSOLVER_H
#define HTUCKER_LINEARSOLVER_LINEARSOLVER_H 1


#include <climits>
#include <limits>
#include <flens/flens.h>
#include <htucker/matrixtensors/matrixtensor.h>


using namespace std;

namespace  htucker{




//--------------------------------------------------------------
//CG
//--------------------------------------------------------------

template <typename T>
int
cg(const HTuckerTree<T> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const T tol = std::numeric_limits<T>::epsilon(), const long maxIterations = LONG_MAX);


//solver versions with an additional rank parameter use truncation!
template <typename T>
int
cg(const HTuckerTree<T> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const int maxrank, const T tol = std::numeric_limits<T>::epsilon(), const long maxIterations = LONG_MAX);



//--------------------------------------------------------------
//Bicgstab
//--------------------------------------------------------------

template <typename T>
int
bicgstab(const HTuckerTree<T> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const T tol = std::numeric_limits<T>::epsilon(), const long maxIterations = LONG_MAX);

//solver versions with an additional rank parameter use truncation!
template <typename T>
int
bicgstab(const HTuckerTree<T> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const int maxrank, const T tol  = std::numeric_limits<T>::epsilon(), const long maxIterations = LONG_MAX);

template <typename T, typename Aop, class A1, class A2>
int
bicgstab(const HTuckerClosure<Aop,A1,A2> &A, HTuckerTree<T> &x,const HTuckerTree<T> &b, const T tol = std::numeric_limits<T>::epsilon(), const long maxIterations = LONG_MAX);

template <typename T, typename Aop, class A1, class A2>
int
bicgstab(const HTuckerClosure<Aop,A1,A2> &A, HTuckerTree<T> &x,const HTuckerTree<T> &b, const int maxrank, const T tol = std::numeric_limits<T>::epsilon(), const long maxIterations = LONG_MAX);

template <typename T, typename Aop, class A1, class A2>
int
bicgstablinf(const HTuckerClosure<Aop,A1,A2> &A, HTuckerTree<T> &x,const HTuckerTree<T> &b, const T tol = std::numeric_limits<T>::epsilon(), const long maxIterations = LONG_MAX);


//--------------------------------------------------------------
//gmres + restarted gmres
//--------------------------------------------------------------


template <typename T>
int
gmres(const HTuckerTree<T> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const T tol = std::numeric_limits<T>::epsilon(), long maxIterations = -1);

template <typename T>
int
gmres(const HTuckerTree<T> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const int maxrank, const T tol = std::numeric_limits<T>::epsilon(), long maxIterations = -1);


template <typename T>
int
gmresm(const HTuckerTree<T> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const T tol=std::numeric_limits<T>::epsilon(),
       const int restart=20, long maxIterations=-1);

template <typename T>
int
gmresm(const HTuckerTree<T> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const int maxrank, const T tol=std::numeric_limits<T>::epsilon(),
       const int restart=20, long maxIterations=-1);

template <typename T, typename Aop, class A1, class A2>
int
gmres(const HTuckerClosure<Aop,A1,A2> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const T tol = std::numeric_limits<T>::epsilon(), long maxIterations = -1);


//solver versions with an additional rank parameter use truncation!
template <typename T, typename Aop, class A1, class A2>
int
gmres(const HTuckerClosure<Aop,A1,A2> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const int maxrank, const T tol = std::numeric_limits<T>::epsilon(), long maxIterations = -1);


template <typename T, typename Aop, class A1, class A2>
int
gmresm(const HTuckerClosure<Aop,A1,A2> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const T tol=std::numeric_limits<T>::epsilon(),
       const int restart=20, long maxIterations=-1);


//solver versions with an additional rank parameter use truncation!
template <typename T, typename Aop, class A1, class A2>
int
gmresm(const HTuckerClosure<Aop,A1,A2> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const int maxrank, const T tol=std::numeric_limits<T>::epsilon(),
       const int restart=20, long maxIterations=-1);






} //namespace htucker

#include <htucker/linearsolver/linearsolver.tcc>
#include <htucker/linearsolver/maxsolver.h>

#endif // HTUCKER_LINEARSOLVER_LINEARSOLVER_H
