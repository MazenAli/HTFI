
namespace htucker{


//The following operators build up the operator tree structures
	
template <typename T>
HTuckerClosure<flens::OpMat,HTuckerTree<T>, HTuckerTree<T> >
mat(const HTuckerTree<T> & tree1){
	return HTuckerClosure<flens::OpMat,HTuckerTree<T>, HTuckerTree<T> >(tree1,tree1,tree1.dim());
};

template <typename T>
HTuckerClosure<flens::OpVec,HTuckerTree<T>, HTuckerTree<T> >
vec(const HTuckerTree<T> & tree1){
	return HTuckerClosure<flens::OpVec,HTuckerTree<T>, HTuckerTree<T> >(tree1,tree1,tree1.dim());
};


//flens::OpAdd

template <typename op1, class A1, class B1, typename op2, class A2, class B2>
HTuckerClosure<flens::OpAdd,HTuckerClosure<op1,A1,B1>,HTuckerClosure<op2,A2,B2> >
operator+ (const HTuckerClosure<op1,A1,B1> & htc1, const HTuckerClosure<op2,A2,B2> & htc2){
	assert(htc1.dim() == htc2.dim());
	return HTuckerClosure<flens::OpAdd,HTuckerClosure<op1,A1,B1>,HTuckerClosure<op2,A2,B2> >(htc1,htc2,htc1.dim());
};

template <typename op, class A, class B>
HTuckerClosure<flens::OpAdd, HTuckerClosure<op,A,B>, IdentityTensor>
operator+ (const HTuckerClosure<op,A,B> & htc, const IdentityTensor &it){
	assert(it.maxvals.length() ==  htc.dim);
	return HTuckerClosure<flens::OpAdd, HTuckerClosure<op,A,B>, IdentityTensor>(htc,it,htc.dim());
};

template <typename op, class A, class B, typename mattype>
HTuckerClosure<flens::OpAdd, HTuckerClosure<op,A,B>, MatrixTensor<mattype> >
operator+ (const HTuckerClosure<op,A,B> & htc, const MatrixTensor<mattype> & mt){
	assert(mt.dim() == htc.dim);
	return HTuckerClosure<flens::OpAdd, HTuckerClosure<op,A,B>, MatrixTensor<mattype> >(htc,mt,htc.dim());
};

template <typename op, class A, class B>
HTuckerClosure<flens::OpAdd, IdentityTensor, HTuckerClosure<op,A,B> >
operator+ (const IdentityTensor &it, const HTuckerClosure<op,A,B> & htc){
	assert(it.maxvals.length() == htc.dim);
	return HTuckerClosure<flens::OpAdd, IdentityTensor, HTuckerClosure<op,A,B> >(it, htc,htc.dim());
};

template <typename op, class A, class B, typename mattype>
HTuckerClosure<flens::OpAdd, MatrixTensor<mattype>,  HTuckerClosure<op,A,B> >
operator+ ( const MatrixTensor<mattype> & mt, const HTuckerClosure<op,A,B> & htc){
	assert(mt.dim() == htc.dim);
	return HTuckerClosure<flens::OpAdd, MatrixTensor<mattype>,  HTuckerClosure<op,A,B> >(mt,htc,htc.dim());
};

template <typename mattype1, typename mattype2>
HTuckerClosure<flens::OpAdd, MatrixTensor<mattype1>, MatrixTensor<mattype2> >
operator + ( const MatrixTensor<mattype1> & mt1, const MatrixTensor<mattype2> & mt2){
	assert(mt1.dim() == mt2.dim());
	return HTuckerClosure<flens::OpAdd, MatrixTensor<mattype1>, MatrixTensor<mattype2> >(mt1,mt2,mt1.dim());
};

HTuckerClosure<flens::OpAdd,IdentityTensor, IdentityTensor>
operator+ (const IdentityTensor &it1, const IdentityTensor &it2){
	assert(it1.maxvals.length()==it2.maxvals.length());
	return HTuckerClosure<flens::OpAdd,IdentityTensor, IdentityTensor>(it1,it2,it1.dim());
};

template <typename mattype>
HTuckerClosure<flens::OpAdd,IdentityTensor, MatrixTensor<mattype> >
operator+ (const IdentityTensor &it, const MatrixTensor<mattype> & mt){
	assert(it.maxvals.length() == mt.dim());
	return HTuckerClosure<flens::OpAdd,IdentityTensor, MatrixTensor<mattype> >(it,mt,mt.dim());
};

template <typename mattype>
HTuckerClosure<flens::OpAdd, MatrixTensor<mattype>,IdentityTensor >
operator+ (const MatrixTensor<mattype> & mt, const IdentityTensor &it){
	assert(mt.dim() == it.maxvals.length());
	return HTuckerClosure<flens::OpAdd, MatrixTensor<mattype>,IdentityTensor >(mt,it,mt.dim());
};

//flens::OpSub

template <typename op1, class A1, class B1, typename op2, class A2, class B2>
HTuckerClosure<flens::OpSub,HTuckerClosure<op1,A1,B1>,HTuckerClosure<op2,A2,B2> >
operator- (const HTuckerClosure<op1,A1,B1> & htc1, const HTuckerClosure<op2,A2,B2> & htc2){
	assert(htc1.dim == htc2.dim);
	return HTuckerClosure<flens::OpSub,HTuckerClosure<op1,A1,B1>,HTuckerClosure<op2,A2,B2> >(htc1,htc2,htc1.dim());
};

template <typename op, class A, class B>
HTuckerClosure<flens::OpSub, HTuckerClosure<op,A,B>, IdentityTensor>
operator- (const HTuckerClosure<op,A,B> & htc, const IdentityTensor &it){
	assert(htc.dim == it.maxvals.length());
	return HTuckerClosure<flens::OpSub, HTuckerClosure<op,A,B>, IdentityTensor>(htc,it,htc.dim());
};

template <typename op, class A, class B, typename mattype>
HTuckerClosure<flens::OpSub, HTuckerClosure<op,A,B>, MatrixTensor<mattype> >
operator- (const HTuckerClosure<op,A,B> & htc, const MatrixTensor<mattype> & mt){
	assert(htc.dim == mt.dim());
	return HTuckerClosure<flens::OpSub, HTuckerClosure<op,A,B>, MatrixTensor<mattype> >(htc,mt,htc.dim());
};

template <typename op, class A, class B>
HTuckerClosure<flens::OpSub, IdentityTensor, HTuckerClosure<op,A,B> >
operator- (const IdentityTensor &it, const HTuckerClosure<op,A,B> & htc){
	assert(it.maxvals.length() == htc.dim);
	return HTuckerClosure<flens::OpSub, IdentityTensor, HTuckerClosure<op,A,B> >(it, htc,htc.dim());
};

template <typename op, class A, class B, typename mattype>
HTuckerClosure<flens::OpSub, MatrixTensor<mattype>,  HTuckerClosure<op,A,B> >
operator- (const MatrixTensor<mattype> & mt, const HTuckerClosure<op,A,B> & htc){
	assert(mt.dim() == htc.dim);
	return HTuckerClosure<flens::OpSub, MatrixTensor<mattype>,  HTuckerClosure<op,A,B> >(mt,htc,htc.dim());
};

template <typename mattype1, typename mattype2>
HTuckerClosure<flens::OpSub, MatrixTensor<mattype1>, MatrixTensor<mattype2> >
operator - (const MatrixTensor<mattype1> & mt1, const MatrixTensor<mattype2> & mt2){
	assert(mt1.dim() == mt2.dim());
	return HTuckerClosure<flens::OpSub, MatrixTensor<mattype1>, MatrixTensor<mattype2> >(mt1,mt2,mt1.dim());
};

HTuckerClosure<flens::OpSub,IdentityTensor, IdentityTensor>
operator- (const IdentityTensor &it1, const IdentityTensor &it2){
	assert(it1.maxvals.length() == it2.maxvals.length());
	return HTuckerClosure<flens::OpSub,IdentityTensor, IdentityTensor>(it1,it2,it1.maxvals.length());
};

template <typename mattype>
HTuckerClosure<flens::OpSub,IdentityTensor, MatrixTensor<mattype> >
operator- (const IdentityTensor &it, const MatrixTensor<mattype> & mt){
	assert(it.maxvals.length() == mt.dim());
	return HTuckerClosure<flens::OpSub,IdentityTensor, MatrixTensor<mattype> >(it,mt,mt.dim());
};

template <typename mattype>
HTuckerClosure<flens::OpSub, MatrixTensor<mattype>,IdentityTensor >
operator- (const MatrixTensor<mattype> & mt, const IdentityTensor &it){
	assert(mt.dim() == it.maxvals.length());
	return HTuckerClosure<flens::OpSub, MatrixTensor<mattype>,IdentityTensor >(mt,it,mt.dim());
};


//OPTensor


template <typename op1, class A1, class B1, typename op2, class A2, class B2>
HTuckerClosure<flens::OpTensor,HTuckerClosure<op1,A1,B1>,HTuckerClosure<op2,A2,B2> >
operator* (const HTuckerClosure<op1,A1,B1> & htc1, const HTuckerClosure<op2,A2,B2> & htc2){
	return HTuckerClosure<flens::OpTensor,HTuckerClosure<op1,A1,B1>,HTuckerClosure<op2,A2,B2> >(htc1,htc2,htc1.dim()+htc2.dim());
};

template <typename op, class A, class B>
HTuckerClosure<flens::OpTensor, HTuckerClosure<op,A,B>, IdentityTensor>
operator* (const HTuckerClosure<op,A,B> & htc, const IdentityTensor &it){
	return HTuckerClosure<flens::OpTensor, HTuckerClosure<op,A,B>, IdentityTensor>(htc,it,htc.dim() + it.dim());
};

template <typename op, class A, class B, typename mattype>
HTuckerClosure<flens::OpTensor, HTuckerClosure<op,A,B>, MatrixTensor<mattype> >
operator* (const HTuckerClosure<op,A,B> & htc, const MatrixTensor<mattype> & mt){
	return HTuckerClosure<flens::OpTensor, HTuckerClosure<op,A,B>, MatrixTensor<mattype> >(htc,mt,htc.dim() + mt.dim());
};

template <typename op, class A, class B>
HTuckerClosure<flens::OpTensor, IdentityTensor, HTuckerClosure<op,A,B> >
operator* (const IdentityTensor &it, const HTuckerClosure<op,A,B> & htc){
	return HTuckerClosure<flens::OpTensor, IdentityTensor, HTuckerClosure<op,A,B> >(it, htc,it.maxvals.length() + htc.dim());
};

template <typename op, class A, class B, typename mattype>
HTuckerClosure<flens::OpTensor, MatrixTensor<mattype>,  HTuckerClosure<op,A,B> >
operator* (const MatrixTensor<mattype> & mt, const HTuckerClosure<op,A,B> & htc){
	return HTuckerClosure<flens::OpTensor, MatrixTensor<mattype>,  HTuckerClosure<op,A,B> >(mt,htc,mt.dim() + htc.dim());
};

template <typename mattype1, typename mattype2>
HTuckerClosure<flens::OpTensor, MatrixTensor<mattype1>, MatrixTensor<mattype2> >
operator * (const MatrixTensor<mattype1> & mt1, const MatrixTensor<mattype2> & mt2){
	return HTuckerClosure<flens::OpTensor, MatrixTensor<mattype1>, MatrixTensor<mattype2> >(mt1,mt2,mt1.dim() + mt2.dim());
};

HTuckerClosure<flens::OpTensor,IdentityTensor, IdentityTensor>
operator* (const IdentityTensor &it1, const IdentityTensor &it2){
	return HTuckerClosure<flens::OpTensor,IdentityTensor, IdentityTensor>(it1,it2,it1.dim()+ it2.dim());
};

template <typename mattype>
HTuckerClosure<flens::OpTensor,IdentityTensor, MatrixTensor<mattype> >
operator* (const IdentityTensor &it, const MatrixTensor<mattype> & mt){
	return HTuckerClosure<flens::OpTensor,IdentityTensor, MatrixTensor<mattype> >(it,mt,it.dim() + mt.dim());
};

template <typename mattype>
HTuckerClosure<flens::OpTensor, MatrixTensor<mattype>,IdentityTensor >
operator* (const MatrixTensor<mattype> & mt, const IdentityTensor &it){
	return HTuckerClosure<flens::OpTensor, MatrixTensor<mattype>,IdentityTensor >(mt,it,mt.dim()+it.dim());
};

template <typename T, typename op, class A, class B>
HTuckerClosure<flens::OpTensor, HTuckerClosure<op,A,B>, VectorTensor<T> >
operator* (const HTuckerClosure<op,A,B> & htc, const VectorTensor<T> &vt){
	return HTuckerClosure<flens::OpTensor, HTuckerClosure<op,A,B>, VectorTensor<T> >(htc,vt,htc.dim()+vt.dim());
};

//-------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------


//This operators evaluate the tree structure

//addition and subtraction makes not much sense here, here we baisically implement the matrix-vector mutliplication

//OpMult

template <typename T>
HTuckerTree<T>
operator*(const HTuckerClosure<flens::OpMat,HTuckerTree<T>,HTuckerTree<T> > & htc, const HTuckerTree<T> & tree){
	return htc.getLeft()*tree;
};

template <class A, class B,typename T> 
HTuckerTree<T>
operator*(const HTuckerClosure<flens::OpAdd,A,B> & htc, const HTuckerTree<T> & tree){
	//std::cout << "operator*(HTC<Add,A,B>, HT>) return htc.left*ht+ htc.right*ht" <<std::endl;
	return htc.getLeft()*tree + htc.getRight()*tree;
};

template <class A, class B,typename T> 
HTuckerTree<T>
operator*(const HTuckerClosure<flens::OpSub,A,B> & htc, const HTuckerTree<T> & tree){
	return htc.getLeft()*tree - htc.getRight()*tree;
};


template <class A, class B,typename T> 
HTuckerTree<T>
operator*(const HTuckerClosure<flens::OpTensor,A,B> & htc, const HTuckerTree<T> & tree){
	//std::cout << "operator*(HTC<Tens>,HT)" << std::endl;
	HTuckerTreePart<T> p1(tree,1,htc.getLeft().dim());
	HTuckerTreePart<T> p2(tree,htc.getLeft().dim()+1,htc.getLeft().dim()+htc.getRight().dim());
	HTuckerTree<T> ret = concat(htc.getLeft()*p1,htc.getRight()*p2);
	ret.getGeneralTree().root->getContent()->setUorB(tree.getGeneralTree().root->getContent()->getUorB());
	ret.getGeneralTree().root->getContent()->setNumRows(tree.getGeneralTree().root->getContent()->getNumRows());
	ret.getGeneralTree().root->getContent()->setLeftChildNumRows(tree.getGeneralTree().root->getContent()->getLeftChildNumRows());
	ret.getGeneralTree().root->getContent()->setRightChildNumRows(tree.getGeneralTree().root->getContent()->getRightChildNumRows());
	return ret;
};

template <typename T>
HTuckerTree<T>
operator*(const HTuckerClosure<flens::OpMat,HTuckerTree<T>,HTuckerTree<T> > & htc, const HTuckerTreePart<T> & treepart){
	 assert(htc.dim() == treepart.maxdim - treepart.mindim + 1);
	 //std::cout << "operator*(HTC<Mat>,HTPart)" << std::endl;
	 return htc.getLeft()*treepart;
};


//die nächsten 4 sollen ausmultiplizieren, wir wollen keine additionen auf teilbäumen
template <class A1, class A2, class B, typename T> 
HTuckerTree<T>
operator* (const HTuckerClosure<flens::OpTensor,HTuckerClosure<flens::OpAdd,A1,A2> ,B> & htc, const HTuckerTreePart<T> & treepart){
	return htc.getLeft().getLeft()*treepart + htc.getLeft().getRight()*treepart;
};

template <class A1, class A2, class B, typename T> 
HTuckerTree<T>
operator* (const HTuckerClosure<flens::OpTensor,HTuckerClosure<flens::OpSub,A1,A2>,B> & htc, const HTuckerTreePart<T> & treepart){
	return htc.getLeft().getLeft()*treepart - htc.getLeft().getRight()*treepart;
};

template <class A, class B1, class B2, typename T> 
HTuckerTree<T>
operator* (const HTuckerClosure<flens::OpTensor,A,HTuckerClosure<flens::OpAdd,B1,B2> > & htc, const HTuckerTreePart<T> & treepart){
	return htc.getRight().getLeft()*treepart + htc.getRight().getRight()*treepart;
};

template <class A, class B1, class B2, typename T> 
HTuckerTree<T>
operator* (const HTuckerClosure<flens::OpTensor,A,HTuckerClosure<flens::OpSub,B1,B2> > & htc, const HTuckerTreePart<T> & treepart){
	return htc.getRight().getLeft()*treepart - htc.getRight().getRight()*treepart;
};

//ausmultiplizieren ende

//matrix-matrix product

template <typename T>
HTuckerTree<T>
operator*(const HTuckerClosure<flens::OpMat,HTuckerTree<T>,HTuckerTree<T> > & htc, const HTuckerClosure<flens::OpMat,HTuckerTree<T>,HTuckerTree<T> > & htc2){
    using flens::_;

	assert(htc.getLeft().dim() == htc2.getLeft().dim());
	std::cout << "Hallo " << std::endl;
	HTuckerTree<T>  tmp = htc.getLeft();
	GeneralTreeNode<HTuckerTreeNode<T> > *node1,*node2,*nodetmp;
	
    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT1 = htc.getLeft().getGeneralTree().end();
	GeneralTreeIterator<HTuckerTreeNode<T> > TIT2 = htc2.getLeft().getGeneralTree().end();
	GeneralTreeIterator<HTuckerTreeNode<T> > TITtmp = tmp.getGeneralTree().end();
    for(; TIT1 >= htc.getLeft().getGeneralTree().begin(); TIT1--,TIT2--,TITtmp--){
		node1 = TIT1.getNode();
		node2 = TIT2.getNode();
		nodetmp = TITtmp.getNode();
		if(TIT1.getNode()->isLeaf()){
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > & node1UorB = node1->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > & node2UorB = node2->getContent()->getUorB();
			
			int veclen = sqrt(node2UorB.numRows() + 0.0); // Achtung, die erfordert quadratische Matrizen!
			std::cout << "Veclen = " << veclen << std::endl;
			std::cout << node2UorB.numRows() << std::endl;
			int len = node1UorB.numRows() / veclen;
			assert(node1UorB.numRows() == veclen * len);
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Rneu(len*veclen,node1UorB.numCols() * node2UorB.numCols());
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > matversion(len,veclen);
			int node1nc = node1->getContent()->getUorB().numCols();
			int node2nr = node2->getContent()->getUorB().numRows();
			
			for(int i = 1; i <= node1nc; ++i){
				for(int j = 1; j <= veclen; ++j){
					matversion(_,j) = node1UorB(_((j-1)*len + 1, j * len),i);
				}
				for(int l = 1; l <= veclen; ++l){
					Rneu(_((l-1)*veclen+1,l*veclen),_((i-1)*node2UorB.numCols()+1,i * node2UorB.numCols())) = matversion * node2UorB(_((l-1)*veclen+1,l*veclen),_);
				}
			}
			nodetmp->getContent()->setUorB(Rneu);
		} else  {
			int node1numel = node1->getContent()->getNumRows();
			int node1lcnumel = node1->getContent()->getLeftChildNumRows();
			int node1rcnumel = node1->getContent()->getRightChildNumRows();
			int node2numel = node2->getContent()->getNumRows();
			int node2lcnumel= node2->getContent()->getLeftChildNumRows();
			int node2rcnumel = node2->getContent()->getRightChildNumRows();

			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Rneu(node1rcnumel * node1numel *node2rcnumel * node2numel,node1lcnumel * node2lcnumel);
			int newblock_numrows = node1rcnumel * node2rcnumel;
			int newblock_numcols = node1lcnumel * node2lcnumel;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &node1UorB = node1->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &node2UorB = node2->getContent()->getUorB();
			
			for(int i = 1; i <= node1numel; ++i){
				for(int j = 1; j<= node2numel; ++j){
					for(int k = 1; k <= node1rcnumel; ++k){
						for(int l = 1; l <= node1lcnumel; ++l){
						
							Rneu(_((i-1)*(newblock_numrows *node2numel) + (j-1) * newblock_numrows + (k-1) * node2rcnumel + 1,(i-1)*(newblock_numrows * node2numel) + (j-1) * newblock_numrows + k * node2rcnumel),_((l-1) * node2lcnumel + 1,l*node2lcnumel)) = node1UorB((i-1)*node1rcnumel + k, l)*node2UorB(_((j-1) * node2rcnumel + 1,j * node2rcnumel),_);
						}
					}
				}
			}

			nodetmp->getContent()->setUorB(Rneu);
			nodetmp->getContent()->setRightChildNumRows(newblock_numrows);
			nodetmp->getContent()->setLeftChildNumRows(newblock_numcols);
			nodetmp->getContent()->setNumRows(node1numel * node2numel);

		}

	}
    }
	return tmp;
	//return htc.getLeft();
};






template <class A, class B,typename T> 
HTuckerTree<T>
operator*(const HTuckerClosure<flens::OpTensor,A,B> & htc, const HTuckerTreePart<T> & treepart){
	assert(treepart.maxdim-treepart.mindim + 1 == htc.dim());

	HTuckerTreePart<T> p1(treepart.httree,treepart.mindim,treepart.mindim + htc.getLeft().dim() - 1);
	HTuckerTreePart<T> p2(treepart.httree,treepart.mindim + htc.getLeft().dim(),treepart.maxdim);

	HTuckerTree<T> rettree = concat(htc.getLeft()*p1,htc.getRight()*p2);
	
	bool found = false;
	GeneralTreeNode<HTuckerTreeNode<T> > * node;
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = treepart.httree.getGeneralTree().end(); TIT >= treepart.httree.getGeneralTree().begin(); TIT--){
		DimensionIndex& index = TIT.getNode()->getContent()->getIndex();
		if(index[0] == treepart.mindim && index[index.length()-1] == treepart.maxdim){
			found = true;
			node = TIT.getNode();
			break;
		}
	}

	if(found){
		rettree.getGeneralTree.root->getContent().setUorB(node->getContent()->getUorB());
		rettree.getGeneralTree.root->getContent()->UorB_numel = node->getContent()->UorB_numel;
		rettree.getGeneralTree.root->getContent()->UorB_lcnumel = node->getContent()->UorB_lcnumel;
		rettree.getGeneralTree.root->getContent()->UorB_rcnumel = node->getContent()->UorB_rcnumel;
	} else {
		assert(0);
	}

	return rettree;
};

template <typename mattype, typename T>
HTuckerTree<T>
operator*(const MatrixTensor<mattype> & mat, const HTuckerTree<T> & tree){
	assert(mat.dim() == tree.dim());
	HTuckerTree<T> rettree = tree;
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = rettree.getGeneralTree().end(); TIT >= rettree.getGeneralTree().begin(); TIT--){
		if(TIT.getNode()->isLeaf()){
			TIT.getNode()->getContent()->setUorB(mat.getMatrix((TIT.getNode()->getContent()->getIndex())[0])*TIT.getNode()->getContent()->getUorB());
		}
	}
	return rettree;
};

template <typename mattype, typename T>
HTuckerTree<T>
operator*(const MatrixTensor<mattype> & mat, const HTuckerTreePart<T> & treepart){
	assert(treepart.maxdim - treepart.mindim + 1 == mat.dim());
	for(int i = 1; i <= mat.dim(); ++i){
		//std::cout << "i: " << i << "   " << mat.getMatrix(i).numRows() << " x " << mat.getMatrix(i).numCols() << std::endl;
	}
	//std::cout << "treepart: " << std::endl;
	//treepart.httree.print();
	//std::cout << "mindim = " << treepart.mindim << std::endl;
	//std::cout << "maxdim = " << treepart.maxdim << std::endl;
	HTuckerTree<T> rettree = treepart.httree;
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = rettree.getGeneralTree().end(); TIT >= rettree.getGeneralTree().begin(); TIT--){
		if(TIT.getNode()->isLeaf()){
			int dimind = TIT.getNode()->getContent()->getIndex()[0];
			if(dimind >= treepart.mindim && dimind <= treepart.maxdim){
				//std::cout << "dimind = " << dimind << std::endl;
				//std::cout << ((TIT.getNode()->getContent()->getIndex())[0]-treepart.mindim + 1) << std::endl;
				//std::cout << mat.getMatrix((TIT.getNode()->getContent()->getIndex())[0]-treepart.mindim + 1) << std::endl;
				//std::cout << TIT.getNode()->getContent()->getUorB() << std::endl;
				TIT.getNode()->getContent()->setUorB(mat.getMatrix((TIT.getNode()->getContent()->getIndex())[0] - treepart.mindim + 1)*TIT.getNode()->getContent()->getUorB());
			}
		}
	}
	//std::cout << "MatrixTensor* HTP" << std::endl;
	//rettree.print_w_UorB();
	//std::cout << "mindim = " << treepart.mindim << "   maxdim = " << treepart.maxdim << std::endl;
	return subtree(rettree,treepart.mindim,treepart.maxdim);
};

template <typename T>
HTuckerTree<T>
operator*(const IdentityTensor &id, const HTuckerTree<T> & tree){
	assert(id.dim() == tree.dim());
	return tree;
};

template <typename T>
HTuckerTree<T>
operator*(const IdentityTensor &id, const HTuckerTreePart<T> & treepart){
	assert(id.dim() == treepart.maxdim - treepart.mindim + 1);
	return subtree(treepart.httree,treepart.mindim,treepart.maxdim);
};

template <typename T>
HTuckerTree<T>
operator*(const HTuckerTree<T> &tree, const HTuckerTreePart<T> &treepart){
	//std::cout << "operator*(HT,HTP)" << std::endl;
	return tree*subtree(treepart.httree,treepart.mindim,treepart.maxdim);
};









}; //namespace htucker
