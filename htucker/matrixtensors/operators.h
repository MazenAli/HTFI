#ifndef HTUCKER_MATRIXTENSORS_OPERATORS_H
#define HTUCKER_MATRIXTENSORS_OPERATORS_H 1

namespace flens{
	
	struct OpMat{};

	struct OpVec{};

	struct OpTensor{};

};


namespace htucker{


//The following operators build up the operator tree structures
	
template <typename T>
HTuckerClosure<flens::OpMat,HTuckerTree<T>, HTuckerTree<T> >
mat(const HTuckerTree<T> & tree1);

template <typename T>
HTuckerClosure<flens::OpVec,HTuckerTree<T>, HTuckerTree<T> >
vec(const HTuckerTree<T> & tree1);

//flens::OpAdd

template <typename op1, class A1, class B1, typename op2, class A2, class B2>
HTuckerClosure<flens::OpAdd,HTuckerClosure<op1,A1,B1>,HTuckerClosure<op2,A2,B2> >
operator+ (const HTuckerClosure<op1,A1,B1> & htc1, const HTuckerClosure<op2,A2,B2> & htc2);

template <typename op, class A, class B>
HTuckerClosure<flens::OpAdd, HTuckerClosure<op,A,B>, IdentityTensor>
operator+ (const HTuckerClosure<op,A,B> & htc, const IdentityTensor &it);

template <typename op, class A, class B, typename mattype>
HTuckerClosure<flens::OpAdd, HTuckerClosure<op,A,B>, MatrixTensor<mattype> >
operator+ (const HTuckerClosure<op,A,B> & htc, const MatrixTensor<mattype> & mt);

template <typename op, class A, class B>
HTuckerClosure<flens::OpAdd, IdentityTensor, HTuckerClosure<op,A,B> >
operator+ (const IdentityTensor &it, const HTuckerClosure<op,A,B> & htc);

template <typename op, class A, class B, typename mattype>
HTuckerClosure<flens::OpAdd, MatrixTensor<mattype>,  HTuckerClosure<op,A,B> >
operator+ (const MatrixTensor<mattype> & mt, const HTuckerClosure<op,A,B> & htc);

template <typename mattype1, typename mattype2>
HTuckerClosure<flens::OpAdd, MatrixTensor<mattype1>, MatrixTensor<mattype2> >
operator + (const MatrixTensor<mattype1> & mt1, const MatrixTensor<mattype2> & mt2);

HTuckerClosure<flens::OpAdd,IdentityTensor, IdentityTensor>
operator+ (const IdentityTensor &it1, const IdentityTensor &it2);

template <typename mattype>
HTuckerClosure<flens::OpAdd,IdentityTensor, MatrixTensor<mattype> >
operator+ (const IdentityTensor &it, const MatrixTensor<mattype> & mt);

template <typename mattype>
HTuckerClosure<flens::OpAdd, MatrixTensor<mattype>,IdentityTensor >
operator+ (const MatrixTensor<mattype> & mt, const IdentityTensor &it);

//flens::OpSub

template <typename op1, class A1, class B1, typename op2, class A2, class B2>
HTuckerClosure<flens::OpSub,HTuckerClosure<op1,A1,B1>,HTuckerClosure<op2,A2,B2> >
operator- (const HTuckerClosure<op1,A1,B1> & htc1, const HTuckerClosure<op2,A2,B2> & htc2);

template <typename op, class A, class B>
HTuckerClosure<flens::OpSub, HTuckerClosure<op,A,B>, IdentityTensor>
operator- (const HTuckerClosure<op,A,B> & htc, const IdentityTensor &it);

template <typename op, class A, class B, typename mattype>
HTuckerClosure<flens::OpSub, HTuckerClosure<op,A,B>, MatrixTensor<mattype> >
operator- (const HTuckerClosure<op,A,B> & htc, const MatrixTensor<mattype> & mt);

template <typename op, class A, class B>
HTuckerClosure<flens::OpSub, IdentityTensor, HTuckerClosure<op,A,B> >
operator- (const IdentityTensor &it, const HTuckerClosure<op,A,B> & htc);

template <typename op, class A, class B, typename mattype>
HTuckerClosure<flens::OpSub, MatrixTensor<mattype>,  HTuckerClosure<op,A,B> >
operator- (const MatrixTensor<mattype> & mt, const HTuckerClosure<op,A,B> & htc);

template <typename mattype1, typename mattype2>
HTuckerClosure<flens::OpSub, MatrixTensor<mattype1>, MatrixTensor<mattype2> >
operator - (const MatrixTensor<mattype1> & mt1, const MatrixTensor<mattype2> & mt2);

HTuckerClosure<flens::OpSub,IdentityTensor, IdentityTensor>
operator- (const IdentityTensor &it1, const IdentityTensor &it2);

template <typename mattype>
HTuckerClosure<flens::OpSub,IdentityTensor, MatrixTensor<mattype> >
operator- (const IdentityTensor &it, const MatrixTensor<mattype> & mt);

template <typename mattype>
HTuckerClosure<flens::OpSub, MatrixTensor<mattype>,IdentityTensor >
operator- (const MatrixTensor<mattype> & mt, const IdentityTensor &it);


//OPTensor


template <typename op1, class A1, class B1, typename op2, class A2, class B2>
HTuckerClosure<flens::OpTensor,HTuckerClosure<op1,A1,B1>,HTuckerClosure<op2,A2,B2> >
operator* (const HTuckerClosure<op1,A1,B1> & htc1, const HTuckerClosure<op2,A2,B2> & htc2);

template <typename op, class A, class B>
HTuckerClosure<flens::OpTensor, HTuckerClosure<op,A,B>, IdentityTensor>
operator* (const HTuckerClosure<op,A,B> & htc, const IdentityTensor &it);

template <typename op, class A, class B, typename mattype>
HTuckerClosure<flens::OpTensor, HTuckerClosure<op,A,B>, MatrixTensor<mattype> >
operator* (const HTuckerClosure<op,A,B> & htc,const  MatrixTensor<mattype> & mt);

template <typename op, class A, class B>
HTuckerClosure<flens::OpTensor, IdentityTensor, HTuckerClosure<op,A,B> >
operator* (const IdentityTensor &it, const HTuckerClosure<op,A,B> & htc);

template <typename op, class A, class B, typename mattype>
HTuckerClosure<flens::OpTensor, MatrixTensor<mattype>,  HTuckerClosure<op,A,B> >
operator* (const MatrixTensor<mattype> & mt, const HTuckerClosure<op,A,B> & htc);

template <typename mattype1, typename mattype2>
HTuckerClosure<flens::OpTensor, MatrixTensor<mattype1>, MatrixTensor<mattype2> >
operator * (const MatrixTensor<mattype1> & mt1, const MatrixTensor<mattype2> & mt2);

HTuckerClosure<flens::OpTensor,IdentityTensor, IdentityTensor>
operator* (const IdentityTensor &it1, const IdentityTensor &it2);

template <typename mattype>
HTuckerClosure<flens::OpTensor,IdentityTensor, MatrixTensor<mattype> >
operator* (const IdentityTensor &it, const MatrixTensor<mattype> & mt);

template <typename mattype>
HTuckerClosure<flens::OpTensor, MatrixTensor<mattype>,IdentityTensor >
operator* (const MatrixTensor<mattype> & mt, const IdentityTensor &it);

template <typename T, typename op, class A, class B>
HTuckerClosure<flens::OpTensor, HTuckerClosure<op,A,B>, VectorTensor<T> >
operator* (const HTuckerClosure<op,A,B> & htc, const VectorTensor<T> &vt);
//-------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------


//This operators evaluate the tree structure

//addition and subtraction makes not much sense here, here we baisically implement the matrix-vector mutliplication

//OpMult

template <typename T>
HTuckerTree<T>
operator*(const HTuckerClosure<flens::OpMat,HTuckerTree<T>,HTuckerTree<T> > & htc, const HTuckerTree<T> & tree);

template <class A, class B,typename T> 
HTuckerTree<T>
operator*(const HTuckerClosure<flens::OpAdd,A,B> & htc, const HTuckerTree<T> & tree);

template <class A, class B,typename T> 
HTuckerTree<T>
operator*(const HTuckerClosure<flens::OpSub,A,B> & htc, const HTuckerTree<T> & tree);


template <class A, class B,typename T> 
HTuckerTree<T>
operator*(const HTuckerClosure<flens::OpTensor,A,B> & htc, const HTuckerTree<T> & tree);

template <typename T>
HTuckerTree<T>
operator*(const HTuckerClosure<flens::OpMat,HTuckerTree<T>,HTuckerTree<T> > & htc, const HTuckerTreePart<T> & treepart);


//die nächsten 4 sollen ausmultiplizieren, wir wollen keine additionen auf teilbäumen
template <class A1, class A2, class B, typename T> 
HTuckerTree<T>
operator* (const HTuckerClosure<flens::OpTensor,HTuckerClosure<flens::OpAdd,A1,A2> ,B> & htc, const HTuckerTreePart<T> & treepart);

template <class A1, class A2, class B, typename T> 
HTuckerTree<T>
operator* (const HTuckerClosure<flens::OpTensor,HTuckerClosure<flens::OpSub,A1,A2>,B> & htc,const HTuckerTreePart<T> & treepart);

template <class A, class B1, class B2, typename T> 
HTuckerTree<T>
operator* (const HTuckerClosure<flens::OpTensor,A,HTuckerClosure<flens::OpAdd,B1,B2> > & htc, const HTuckerTreePart<T> & treepart);

template <class A, class B1, class B2, typename T> 
HTuckerTree<T>
operator* (const HTuckerClosure<flens::OpTensor,A,HTuckerClosure<flens::OpSub,B1,B2> > & htc, const HTuckerTreePart<T> & treepart);

//ausmultiplizieren ende

//matrix-matrix product
template <typename T>
HTuckerTree<T>
operator*(const HTuckerClosure<flens::OpMat,HTuckerTree<T>,HTuckerTree<T> > & htc, const HTuckerClosure<flens::OpMat,HTuckerTree<T>,HTuckerTree<T> > & htc2);




template <class A, class B,typename T> 
HTuckerTree<T>
operator*(const HTuckerClosure<flens::OpTensor,A,B> & htc, const HTuckerTreePart<T> & treepart);

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
