#ifndef HTUCKER_HTUCKERTREE_HTUCKERTREE_H
#define HTUCKER_HTUCKERTREE_HTUCKERTREE_H 1

#include <cdo/andisspecialtools/matrix_tools.h>
#include <htucker/tensor/tensor.h>
#include <htucker/htuckertree/htuckerelement.h>
#include <htucker/htuckertree/htuckertreenode.h>
#include <cdo/andisspecialtools/gnuplot/gnuplot_i.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <htucker/rankoneapproximation/roa.h>

namespace htucker{


template <typename T>
class HTuckerTree{
	
	GeneralTree<HTuckerTreeNode<T> > tree;
	int d;
	
	void init(const int _d,const double _split);

	template <typename TensorFunction>
	void setUorB(const GeneralTree<ROA<T,TensorFunction> > & gt);

public:

	HTuckerTree(){};

	HTuckerTree(const int _d);

	HTuckerTree(const int _d, const double _split);

    void
    set_tree(const HTuckerTree<T>& X);

	template <typename NType>
	GeneralTree<NType>
	copy() const;

	void print() const;

	void print_info() const;

	void print_w_UorB() const;

    void
    print_svs(bool isorth=false);

	//please call this function only for vectors of small length!!!
	void print_values() const;

	const GeneralTree<HTuckerTreeNode<T> > &
	getGeneralTree() const;

	GeneralTree<HTuckerTreeNode<T> > &
	getGeneralTree();

	T evaluate(const DimensionIndex &index) const; 
	
	flens::DenseVector<flens::Array<T> >
	vec_evaluate(const DimensionIndex &index, const int vardim) const;

    void orthogonalize();

    void orthogonalize_svd(std::vector
                           <flens::DenseVector
                           <flens::Array<T> > >& sigmas,
                           bool isorth = false);

	T L2norm() const;

	T L2normorthogonal() const;

	T ScalarProduct(const HTuckerTree<T> & anothertree) const;

	//Generate an H-Tucker version of the Tensor \sum_{i=1}^k \otimes_\mu^d a_{i,\mu} where a_{i,\mu} is a vector
	//The dense vector list contains these Vectors in the order a11, a21,..., ak1,a12,a22,..,ak2,....,akd
	void generateTofElementary(DenseVectorList<T> &list,const int k,const int d);

	T average_rank() const;

	int max_rank() const;

	int max_n() const;

	T effective_rank() const;

	template <typename TensorFunction>
	void approximate(const tensor<TensorFunction> &tf, const int rank, const int l = 3);

    template <typename TensorFunction>
	void approximate(const tensor<TensorFunction> &tf, const double epsilon, const int l = 3);

    void
    truncate(const int rank, bool isorth = false);

    void
    truncate(double eps, bool isorth = false);

    void
    truncate_hsvd(const T eps);

    void
    truncate(const HTuckerTree<T>& gram, const int rank);

    void
    truncate(const HTuckerTree<T>& gram,
             const flens::DenseVector<flens::Array<T> >& eps);
    void
    truncate_elementary(const HTuckerTree<T>& gram,
                        const flens::DenseVector<flens::Array<T> >& eps);

	T Linfnorm(HTuckerTree<T> & anothertree,const DimensionIndex & minval,const DimensionIndex & maxval,const int n); //not const because evaluation of this and anothertree is needed!

	T Linfnorm(const DimensionIndex & minval, const DimensionIndex & maxval, const int n);  //not const because of evaluate!

	template <typename TensorFunction>
	T Linfnorm(const tensor<TensorFunction> & tens,const DimensionIndex & minval, const DimensionIndex & maxval,const int n);

	void spy(const char* filename, const char* terminal, const char* options, const double width = 1.0, const double height = 1.0, const double pointsize = 1.0) const;

	void spy(const DimensionIndex &subknoten,const char* filename,const char* terminal, const char* options, const double width = 1.0, const double height = 1.0, const double pointsize = 1.0) const;

	void plot_sv(const char* filename, const char* terminal, const char* options, const double width = 1.0, const double height = 1.0, const double pointsize = 1.0) const;

	void plot_sv(const DimensionIndex &subknoten,const char* filename,const char* terminal, const char* options, const double width = 1.0, const double height = 1.0, const double pointsize = 1.0) const;


	DimensionIndex
	getmax() const;

	int dim() const;

	void setdim(const int dim);

    void
    scal(T alpha);

    int
    depth() const;

    HTuckerTree<T>&
    operator=(const HTuckerTree<T>& copy);

};


template <typename T>
HTuckerTree<T>
reapproximate(HTuckerTree<T> & tree1, const double epsilon, const int l = 3);

template <typename T>
HTuckerTree<T>
reapproximate(HTuckerTree<T> & tree1, const int rank,const int l = 3);

template <typename T>
HTuckerTree<T>   add_truncate(const HTuckerTree<T> & tree1,
                              const HTuckerTree<T> & tree2,
                              const double eps);

template <typename T>
HTuckerTree<T>   operator+(const HTuckerTree<T> & tree1, const HTuckerTree<T> & tree2);

template <typename T>
HTuckerTree<T>  operator*(const double number, const HTuckerTree<T> &tree);

template <typename T>
HTuckerTree<T>  operator*(const HTuckerTree<T> &tree, const double number);

template <typename T>
HTuckerTree<T> operator-(const HTuckerTree<T> & tree1, const HTuckerTree<T> & tree2);

template <typename T>
HTuckerTree<T> operator*(const HTuckerTree<T> & tree1, const HTuckerTree<T> & tree2);

template <typename T>
HTuckerTree<T> ones(const DimensionIndex &dims);

template <typename T>
HTuckerTree<T> identity(const DimensionIndex &dims);


template <typename T> 
HTuckerTree<T>
transpose(const HTuckerTree<T> & tree);

template <typename T>
HTuckerTree<T> 
vec2diag(const HTuckerTree<T> &tree);

template <typename T>
HTuckerTree<T> concat(const HTuckerTree<T> &tree1, const HTuckerTree<T> &tree2);

template <typename T>
HTuckerTree<T> subtree(const HTuckerTree<T> &tree, const int mindim, const int maxdim);

template <typename T>
HTuckerTree<T> gramians_orthogonal(const HTuckerTree<T> & tree);

template <typename T>
HTuckerTree<T> gramians_nonorthogonal(const HTuckerTree<T> & tree);

template <typename T>
HTuckerTree<T> gramians_elementary(const HTuckerTree<T>& tree);

template <typename T>
HTuckerTree<T> truncate(const HTuckerTree<T> & tree, const int rank);

template <typename T>
HTuckerTree<T> truncate_nonorthogonal(const HTuckerTree<T> & tree, const int rank);






template <typename T>
class HTuckerTree2Tensor{
	bool isset;
	public:
	typedef T type;
	const HTuckerTree<T> *tree;

	HTuckerTree2Tensor();

	HTuckerTree2Tensor(const int min, const int max, const int d); //dummy for tensor

	HTuckerTree2Tensor(const DimensionIndex & _min, const DimensionIndex &_max); //dummy

	void setTree(const HTuckerTree<T> &_tree);

	T operator()(const DimensionIndex &vals) const;

	bool vecEval() const;

	void  
	vec(const DimensionIndex & vals, const int dim, flens::DenseVector<flens::Array<type> > & vec) const;
};



} // namespace htucker
#include <htucker/htuckertree/htuckertree.tcc>



#endif // HTUCKER_HTUCKERTREE_HTUCKERTREE_H
