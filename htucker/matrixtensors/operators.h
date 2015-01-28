#ifndef HTUCKER_MATRIXTENSORS_OPERATORS_H
#define HTUCKER_MATRIXTENSORS_OPERATORS_H 1

namespace flens{
	
	struct OpMat{};

	struct OpVec{};

	struct OpTensor{};

};

using namespace flens;



namespace htucker{


//The following operators build up the operator tree structures
	
template <typename T>
HTuckerClosure<OpMat,HTuckerTree<T>, HTuckerTree<T> >
mat(const HTuckerTree<T> & tree1);

template <typename T>
HTuckerClosure<OpVec,HTuckerTree<T>, HTuckerTree<T> >
vec(const HTuckerTree<T> & tree1);

//OpAdd

template <typename op1, class A1, class B1, typename op2, class A2, class B2>
HTuckerClosure<OpAdd,HTuckerClosure<op1,A1,B1>,HTuckerClosure<op2,A2,B2> >
operator+ (const HTuckerClosure<op1,A1,B1> & htc1, const HTuckerClosure<op2,A2,B2> & htc2);

template <typename op, class A, class B>
HTuckerClosure<OpAdd, HTuckerClosure<op,A,B>, IdentityTensor>
operator+ (const HTuckerClosure<op,A,B> & htc, const IdentityTensor &it);

template <typename op, class A, class B, typename mattype>
HTuckerClosure<OpAdd, HTuckerClosure<op,A,B>, MatrixTensor<mattype> >
operator+ (const HTuckerClosure<op,A,B> & htc, const MatrixTensor<mattype> & mt);

template <typename op, class A, class B>
HTuckerClosure<OpAdd, IdentityTensor, HTuckerClosure<op,A,B> >
operator+ (const IdentityTensor &it, const HTuckerClosure<op,A,B> & htc);

template <typename op, class A, class B, typename mattype>
HTuckerClosure<OpAdd, MatrixTensor<mattype>,  HTuckerClosure<op,A,B> >
operator+ (const MatrixTensor<mattype> & mt, const HTuckerClosure<op,A,B> & htc);

template <typename mattype1, typename mattype2>
HTuckerClosure<OpAdd, MatrixTensor<mattype1>, MatrixTensor<mattype2> >
operator + (const MatrixTensor<mattype1> & mt1, const MatrixTensor<mattype2> & mt2);

HTuckerClosure<OpAdd,IdentityTensor, IdentityTensor>
operator+ (const IdentityTensor &it1, const IdentityTensor &it2);

template <typename mattype>
HTuckerClosure<OpAdd,IdentityTensor, MatrixTensor<mattype> >
operator+ (const IdentityTensor &it, const MatrixTensor<mattype> & mt);

template <typename mattype>
HTuckerClosure<OpAdd, MatrixTensor<mattype>,IdentityTensor >
operator+ (const MatrixTensor<mattype> & mt, const IdentityTensor &it);

//OpSub

template <typename op1, class A1, class B1, typename op2, class A2, class B2>
HTuckerClosure<OpSub,HTuckerClosure<op1,A1,B1>,HTuckerClosure<op2,A2,B2> >
operator- (const HTuckerClosure<op1,A1,B1> & htc1, const HTuckerClosure<op2,A2,B2> & htc2);

template <typename op, class A, class B>
HTuckerClosure<OpSub, HTuckerClosure<op,A,B>, IdentityTensor>
operator- (const HTuckerClosure<op,A,B> & htc, const IdentityTensor &it);

template <typename op, class A, class B, typename mattype>
HTuckerClosure<OpSub, HTuckerClosure<op,A,B>, MatrixTensor<mattype> >
operator- (const HTuckerClosure<op,A,B> & htc, const MatrixTensor<mattype> & mt);

template <typename op, class A, class B>
HTuckerClosure<OpSub, IdentityTensor, HTuckerClosure<op,A,B> >
operator- (const IdentityTensor &it, const HTuckerClosure<op,A,B> & htc);

template <typename op, class A, class B, typename mattype>
HTuckerClosure<OpSub, MatrixTensor<mattype>,  HTuckerClosure<op,A,B> >
operator- (const MatrixTensor<mattype> & mt, const HTuckerClosure<op,A,B> & htc);

template <typename mattype1, typename mattype2>
HTuckerClosure<OpSub, MatrixTensor<mattype1>, MatrixTensor<mattype2> >
operator - (const MatrixTensor<mattype1> & mt1, const MatrixTensor<mattype2> & mt2);

HTuckerClosure<OpSub,IdentityTensor, IdentityTensor>
operator- (const IdentityTensor &it1, const IdentityTensor &it2);

template <typename mattype>
HTuckerClosure<OpSub,IdentityTensor, MatrixTensor<mattype> >
operator- (const IdentityTensor &it, const MatrixTensor<mattype> & mt);

template <typename mattype>
HTuckerClosure<OpSub, MatrixTensor<mattype>,IdentityTensor >
operator- (const MatrixTensor<mattype> & mt, const IdentityTensor &it);


//OPTensor


template <typename op1, class A1, class B1, typename op2, class A2, class B2>
HTuckerClosure<OpTensor,HTuckerClosure<op1,A1,B1>,HTuckerClosure<op2,A2,B2> >
operator* (const HTuckerClosure<op1,A1,B1> & htc1, const HTuckerClosure<op2,A2,B2> & htc2);

template <typename op, class A, class B>
HTuckerClosure<OpTensor, HTuckerClosure<op,A,B>, IdentityTensor>
operator* (const HTuckerClosure<op,A,B> & htc, const IdentityTensor &it);

template <typename op, class A, class B, typename mattype>
HTuckerClosure<OpTensor, HTuckerClosure<op,A,B>, MatrixTensor<mattype> >
operator* (const HTuckerClosure<op,A,B> & htc,const  MatrixTensor<mattype> & mt);

template <typename op, class A, class B>
HTuckerClosure<OpTensor, IdentityTensor, HTuckerClosure<op,A,B> >
operator* (const IdentityTensor &it, const HTuckerClosure<op,A,B> & htc);

template <typename op, class A, class B, typename mattype>
HTuckerClosure<OpTensor, MatrixTensor<mattype>,  HTuckerClosure<op,A,B> >
operator* (const MatrixTensor<mattype> & mt, const HTuckerClosure<op,A,B> & htc);

template <typename mattype1, typename mattype2>
HTuckerClosure<OpTensor, MatrixTensor<mattype1>, MatrixTensor<mattype2> >
operator * (const MatrixTensor<mattype1> & mt1, const MatrixTensor<mattype2> & mt2);

HTuckerClosure<OpTensor,IdentityTensor, IdentityTensor>
operator* (const IdentityTensor &it1, const IdentityTensor &it2);

template <typename mattype>
HTuckerClosure<OpTensor,IdentityTensor, MatrixTensor<mattype> >
operator* (const IdentityTensor &it, const MatrixTensor<mattype> & mt);

template <typename mattype>
HTuckerClosure<OpTensor, MatrixTensor<mattype>,IdentityTensor >
operator* (const MatrixTensor<mattype> & mt, const IdentityTensor &it);

template <typename T, typename op, class A, class B>
HTuckerClosure<OpTensor, HTuckerClosure<op,A,B>, VectorTensor<T> >
operator* (const HTuckerClosure<op,A,B> & htc, const VectorTensor<T> &vt);
//-------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------


//This operators evaluate the tree structure

//addition and subtraction makes not much sense here, here we baisically implement the matrix-vector mutliplication

//OpMult

template <typename T>
HTuckerTree<T>
operator*(const HTuckerClosure<OpMat,HTuckerTree<T>,HTuckerTree<T> > & htc, const HTuckerTree<T> & tree);

template <class A, class B,typename T> 
HTuckerTree<T>
operator*(const HTuckerClosure<OpAdd,A,B> & htc, const HTuckerTree<T> & tree);

template <class A, class B,typename T> 
HTuckerTree<T>
operator*(const HTuckerClosure<OpSub,A,B> & htc, const HTuckerTree<T> & tree);


template <class A, class B,typename T> 
HTuckerTree<T>
operator*(const HTuckerClosure<OpTensor,A,B> & htc, const HTuckerTree<T> & tree);

template <typename T>
HTuckerTree<T>
operator*(const HTuckerClosure<OpMat,HTuckerTree<T>,HTuckerTree<T> > & htc, const HTuckerTreePart<T> & treepart);


//die nächsten 4 sollen ausmultiplizieren, wir wollen keine additionen auf teilbäumen
template <class A1, class A2, class B, typename T> 
HTuckerTree<T>
operator* (const HTuckerClosure<OpTensor,HTuckerClosure<OpAdd,A1,A2> ,B> & htc, const HTuckerTreePart<T> & treepart);

template <class A1, class A2, class B, typename T> 
HTuckerTree<T>
operator* (const HTuckerClosure<OpTensor,HTuckerClosure<OpSub,A1,A2>,B> & htc,const HTuckerTreePart<T> & treepart);

template <class A, class B1, class B2, typename T> 
HTuckerTree<T>
operator* (const HTuckerClosure<OpTensor,A,HTuckerClosure<OpAdd,B1,B2> > & htc, const HTuckerTreePart<T> & treepart);

template <class A, class B1, class B2, typename T> 
HTuckerTree<T>
operator* (const HTuckerClosure<OpTensor,A,HTuckerClosure<OpSub,B1,B2> > & htc, const HTuckerTreePart<T> & treepart);

//ausmultiplizieren ende

//matrix-matrix product
template <typename T>
HTuckerTree<T>
operator*(const HTuckerClosure<OpMat,HTuckerTree<T>,HTuckerTree<T> > & htc, const HTuckerClosure<OpMat,HTuckerTree<T>,HTuckerTree<T> > & htc2);




template <class A, class B,typename T> 
HTuckerTree<T>
operator*(const HTuckerClosure<OpTensor,A,B> & htc, const HTuckerTreePart<T> & treepart);

template <typename mattype, typename T>
HTuckerTree<T>
operator*(const MatrixTensor<mattype> & mat, const HTuckerTree<T> & tree);

template <typename mattype, typename T>
HTuckerTree<T>
operator*(const MatrixTensor<mattype> & mat, const HTuckerTreePart<T> & treepart);

template <typename T>
HTuckerTree<T>
operator*(const IdentityTensor &id, const HTuckerTree<T> & tree);

template <typename T>
HTuckerTree<T>
operator*(const IdentityTensor &id, const HTuckerTreePart<T> & treepart);

template <typename T>
HTuckerTree<T>
operator*(const HTuckerTree<T> &tree, const HTuckerTreePart<T> &treepart);





};

#include <htucker/matrixtensors/operators.tcc>

#endif // HTUCKER_MATRIXTENSORS_OPERATORS_H
