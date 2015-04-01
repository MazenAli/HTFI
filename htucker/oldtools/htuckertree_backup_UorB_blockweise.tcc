
template <typename T, class EVALTYPE>
void 
HTuckerTree<T,EVALTYPE>::init(int _d,double _split,tensor<T> & _A){
	assert(0 <= _split && _split <= 1);
	
	DimensionIndex all(_d);
	for(int i = 0; i < _d; ++i){
		all[i] = i+1;
	}
	
	HTuckerTreeNode<T,EVALTYPE> rootnode(all,all,_A,_d);
	
	tree = GeneralTree<HTuckerTreeNode<T,EVALTYPE> >(rootnode);
	
	GeneralTreeNode<HTuckerTreeNode<T,EVALTYPE> > * node = tree.root;
	GeneralTreeNode<HTuckerTreeNode<T,EVALTYPE> > * nextlev = NULL;

	while(1){
		int len = (node->getContent()->getIndex()).length();
		if(len > 1){
			int splitpoint = ceil(len * _split);
			if(len <= splitpoint){
				splitpoint = len-1;
			} else if(len <= 0){
				splitpoint = 1;
			}
			
			DimensionIndex fi(splitpoint);
			DimensionIndex la(len - splitpoint);
			for(int i = 0; i< splitpoint; ++i){
				fi[i] = (node->getContent()->getIndex())[i];
			}
			for(int i = 0; i < len - splitpoint; ++i){
				la[i] = (node->getContent()->getIndex())[splitpoint + i];
			}

			HTuckerTreeNode<T,EVALTYPE> fino(fi,la,_A,_d);
			HTuckerTreeNode<T,EVALTYPE> lano(la,fi,_A,_d);
			node->appendChild(fino);
			node->appendChild(lano);
			if(nextlev == 0){
				nextlev = node->getfirstChild();
			}
		}
		if(node->getlevelright() != NULL){
			node = node->getlevelright();
		} else if(nextlev != NULL) {
			node = nextlev;
			nextlev = NULL;
		} else {
			break;
		}
	}
};

template <typename T, class EVALTYPE>
HTuckerTree<T,EVALTYPE>::HTuckerTree(int _d, tensor<T> & _A):d(_d){
	init(_d,0.5,_A);
};

template <typename T, class EVALTYPE>
HTuckerTree<T,EVALTYPE>::HTuckerTree(int _d,double _split, tensor<T> & _A):d(_d){
	init(_d,_split,_A);
};

template <typename T, class EVALTYPE>
GeneralTree<HTuckerTreeNode<T,EVALTYPE> > &
HTuckerTree<T,EVALTYPE>::getGeneralTree(){
	return tree;
};

template <typename T, class EVALTYPE>
void 
HTuckerTree<T,EVALTYPE>::approximate(int rank){
	GeneralTreeNode<HTuckerTreeNode<T,EVALTYPE> > *node;
	for(GeneralTreeIterator<HTuckerTreeNode<T,EVALTYPE> > TIT = tree.begin(); TIT <= tree.end(); TIT++){
		node = TIT.getNode();
		if(!(node->isRoot())){
			std::cout << "HTUCKERTREE.tcc: approximate " << node->getContent()->getIndex() << std::endl;
			std::cout << "HTUCKERTREE.tcc: roa " << std::endl;
			node->getContent()->getROA().print();
			node->getContent()->getROA().approximate(rank,node->getParent()->getContent()->getROA().pivots);
		}
	}
};


template <typename T, class EVALTYPE>
void 
HTuckerTree<T,EVALTYPE>::approximate(double epsilon){
	GeneralTreeNode<HTuckerTreeNode<T,EVALTYPE> > *node;
	for(GeneralTreeIterator<HTuckerTreeNode<T,EVALTYPE> > TIT = tree.begin(); TIT <= tree.end(); TIT++){
		node = TIT.getNode();
		if(!(node->isRoot())){
			std::cout << "HTUCKERTREE.tcc: approximate " << node->getContent()->getIndex() << std::endl;
			std::cout << "HTUCKERTREE.tcc: roa " << std::endl;
			node->getContent()->getROA().print();
			node->getContent()->getROA().approximate(epsilon,node->getParent()->getContent()->getROA().pivots);
		}
	}
};




template <typename T, class EVALTYPE>
void 
HTuckerTree<T,EVALTYPE>::print(){
	tree.print();
};

template <typename T, class EVALTYPE>
void 
HTuckerTree<T,EVALTYPE>::print_w_UorB(){
	tree.print();
	for(GeneralTreeIterator<HTuckerTreeNode<T,EVALTYPE> > TIT = tree.begin(); TIT <= tree.end(); TIT++){
		std::cout << "Node " << TIT.getNode()->getContent()->getIndex() << std::endl;
		std::cout << "UorB " << TIT.getNode()->getContent()->getUorB() << std::endl;
	}

};


template <typename T, class EVALTYPE>
T 
HTuckerTree<T,EVALTYPE>::evaluate(DimensionIndex index){
	T ret = 0.0;
	GeneralTreeNode<HTuckerTreeNode<T,EVALTYPE> > *node;
	flens::DenseVector<flens::Array<T> > evaluate;
	for(GeneralTreeIterator<HTuckerTreeNode<T,EVALTYPE> > TIT = tree.end(); TIT >= tree.begin(); TIT--){
		node = TIT.getNode();
		//std::cout << "htucker.tcc Evaluate at node " << node->getContent()->getIndex() << std::endl;
		//std::cout << "htucker.tcc UorB " << node->getContent()->getUorB() << std::endl;
		if(node->isRoot()){
			evaluate = flens::DenseVector<flens::Array<T> >(1);
			flens::DenseVector<flens::Array<T> > intermediateresult;
			intermediateresult =node->getContent()->getUorB() * node->getfirstChild()->getContent()->getEvaluate();
			evaluate(1) = 0;
			flens::DenseVector<flens::Array<T> > evaluate_right = node->getlastChild()->getContent()->getEvaluate();
			for(int j = 1; j<= intermediateresult.length(); ++j){
				evaluate(1) += intermediateresult(j)*evaluate_right(j);
			}
			return evaluate(1);		
		} else if(node->isLeaf()){
			assert(node->getContent()->getIndex().length() == 1);
			int numcols = node->getContent()->getUorB().numCols(); // in den Blättern geht das so, nicht aber in den inneren Knoten!
			evaluate = flens::DenseVector<flens::Array<T> >(numcols);
			int val = index[(node->getContent()->getIndex())[0] - 1] - node->getContent()->getROA().A.minval[(node->getContent()->getIndex())[0] - 1];
			evaluate = (node->getContent()->getUorB())(val+1,_);
			node->getContent()->setEvaluate(evaluate);
		} else{
			int numblocks = node->getContent()->UorB_numel;
			evaluate = flens::DenseVector<flens::Array<T> >(numblocks);
			int lcnumel = node->getContent()->UorB_lcnumel;
			int rcnumel = node->getContent()->UorB_rcnumel;

			flens::DenseVector<flens::Array <T> > intermediateresult;
			
			for(int i= 1; i<= numblocks; ++i){
				intermediateresult = (node->getContent()->getUorB())(_((i-1)*rcnumel + 1,i*rcnumel),_)  * node->getfirstChild()->getContent()->getEvaluate();
				evaluate(i) = 0;
				for(int j = 1; j<= intermediateresult.length(); ++j){
					evaluate(i) += intermediateresult(j)*(node->getlastChild()->getContent()->getEvaluate())(j);
				}
			}
			node->getContent()->setEvaluate(evaluate);
		}
	}
};



template <typename T, class EVALTYPE>
void
HTuckerTree<T,EVALTYPE>::orthogonalize(){
	GeneralTreeNode<HTuckerTreeNode<T,EVALTYPE> > *node;
	for(GeneralTreeIterator<HTuckerTreeNode<T,EVALTYPE> > TIT = tree.end(); TIT >= tree.begin(); TIT--){
		node = TIT.getNode();
		//std::cout << "node = " << node->getContent()->getIndex() << std::endl;
		if(node->isLeaf()){
			//nothing happens
		} else {
			//compute QR of left child
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Rl,Rr, Ql,Qr;
			flens::DenseVector<flens::Array<T> > taul,taur;
			int numRows,numCols;
			if(node->getfirstChild()->isLeaf()){
				Rl = node->getfirstChild()->getContent()->getUorB();
			} else {
				int numel = node->getfirstChild()->getContent()->UorB_numel;
				int lcnumel = node->getfirstChild()->getContent()->UorB_lcnumel;
				int rcnumel = node->getfirstChild()->getContent()->UorB_rcnumel;
				
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Usave;
				Usave = node->getfirstChild()->getContent()->getUorB();
				Rl = flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> >(lcnumel * rcnumel,numel);
				for(int i = 1; i <= lcnumel; ++i){
					for(int j = 1; j <= numel; ++j){
						Rl(_((i-1)*rcnumel + 1,i*rcnumel),j) = Usave(_((j-1)*rcnumel + 1,j*rcnumel),i);
					}
				}
			}
			
			numCols = Rl.numCols();
			numRows = Rl.numRows();
			qrf(Rl,taul);
			Ql = Rl(_(1,numRows),_(1,min(numRows,numCols)));
			orgqr(Ql,taul);

			//compute QR of right child

			if(node->getlastChild()->isLeaf()){
				Rr = node->getlastChild()->getContent()->getUorB();
			} else {
				int numel = node->getlastChild()->getContent()->UorB_numel;
				int lcnumel = node->getlastChild()->getContent()->UorB_lcnumel;
				int rcnumel = node->getlastChild()->getContent()->UorB_rcnumel;
				
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Usave;
				Usave = node->getlastChild()->getContent()->getUorB();
				Rr = flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> >(lcnumel * rcnumel,numel);
				for(int i = 1; i <= lcnumel; ++i){
					for(int j = 1; j <= numel; ++j){
						Rr(_((i-1)*rcnumel + 1,i*rcnumel),j) = Usave(_((j-1)*rcnumel + 1,j*rcnumel),i);
					}
				}
			}
			
			numCols = Rr.numCols();
			numRows = Rr.numRows();
			qrf(Rr,taur);
			Qr = Rr(_(1,numRows),_(1,min(numRows,numCols)));
			orgqr(Qr,taur);

			//save orthogonalized childs and adapt this node (if the child is an inner node, we have to retransform into block-format)

			if(node->getfirstChild()->isLeaf()){
				node->getfirstChild()->getContent()->setUorB(Ql);
			} else {
				int numel = node->getfirstChild()->getContent()->UorB_numel;
				int lcnumel = node->getfirstChild()->getContent()->UorB_lcnumel;
				int rcnumel = node->getfirstChild()->getContent()->UorB_rcnumel;
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Qneu(Ql.numCols()*rcnumel,lcnumel);
				for(int i = 1; i <= Ql.numCols(); ++i){
					for(int j = 1; j <= lcnumel; ++j){
						Qneu(_((i-1)*rcnumel + 1, i*rcnumel),j) = Ql(_((j-1)*rcnumel +1, j*rcnumel),i);
					}
				}
				node->getfirstChild()->getContent()->UorB_numel = Ql.numCols();
				node->getfirstChild()->getContent()->setUorB(Qneu);
			}

			if(node->getlastChild()->isLeaf()){
				node->getlastChild()->getContent()->setUorB(Ql);
			} else {
				int numel = node->getlastChild()->getContent()->UorB_numel;
				int lcnumel = node->getlastChild()->getContent()->UorB_lcnumel;
				int rcnumel = node->getlastChild()->getContent()->UorB_rcnumel;
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Qneu(Qr.numCols()*rcnumel,lcnumel);
				for(int i = 1; i <= Qr.numCols(); ++i){
					for(int j = 1; j <= lcnumel; ++j){
						Qneu(_((i-1)*rcnumel + 1, i*rcnumel),j) = Qr(_((j-1)*rcnumel +1, j*rcnumel),i);
					}
				}
				node->getlastChild()->getContent()->UorB_numel = Qr.numCols();
				node->getlastChild()->getContent()->setUorB(Qneu);
			}

			//transform this node with the triangular-matrices
			int numel = node->getContent()->UorB_numel;
			int lcnumel = node->getContent()->UorB_lcnumel;
			int rcnumel = node->getContent()->UorB_rcnumel;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Rneu(numel*Rr.numRows(),Rl.numRows());
			for(int i = 1; i<= numel; ++i){
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > temp1(Rr.numRows(),Rl.numRows());
				for(int j = 1; j<= Rr.numRows(); j++){
					for(int k = 1; k <= Rl.numRows(); ++k){
						temp1(j,k) = 0;
						for(int mu = j; mu <= Rr.numCols(); ++mu){
							for(int nu = k; nu <= Rl.numCols(); ++nu){
								temp1(j,k) += Rr(j,mu)*Rl(k,nu)*(node->getContent()->getUorB())((i-1)*rcnumel + mu,nu);
							}
						}
					}
				}
				// die vierfachschleife oben soll die folgenden zwei Zeilen ersetzten, weil für Rr.upper() offensichtlich kein mm definiert ist.
				//temp1 = Rr.upper() * (node->getContent()->getUorB())(_((i-1)*rcnumel + 1, i*rcnumel),_);
				//temp2 = temp1*transpose(Rl.upper());
				Rneu(_((i-1)*Rr.numRows() + 1, i*Rr.numRows()),_) = temp1;
			}
			node->getContent()->UorB_lcnumel = Rl.numRows();
			node->getContent()->UorB_rcnumel = Rr.numRows();
			node->getContent()->setUorB(Rneu);

		}
	}
};


template <typename T, class EVALTYPE>
T
HTuckerTree<T,EVALTYPE>::L2norm(){
	T erg = ScalarProduct(*this);
	std::cout << "htuckertree.tcc in L2norm: ScalarProduct: " << erg << std::endl;
	return sqrt(erg);
};

template <typename T, class EVALTYPE>
T
HTuckerTree<T,EVALTYPE>::L2normorthogonal(){
	T erg;
	/*flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > U,VT;
	flens::DenseVector<flens::Array<T> > s;
	svd(tree.root->getContent()->getUorB(),s,U,VT);
	for(int i = 1; i <= s.length(); ++i){
		erg += s(i);
	}
	return sqrt(erg);*/
	for(int i = 1; i <= tree.root->getContent()->getUorB().numRows(); ++i){
		for(int j = 1; j<= tree.root->getContent()->getUorB().numCols(); ++j){
		erg += ( tree.root->getContent()->getUorB())(i,j)*( tree.root->getContent()->getUorB())(i,j);
		}
	}
	return sqrt(erg);
};

template <typename T, class EVALTYPE>
T
HTuckerTree<T,EVALTYPE>::ScalarProduct(HTuckerTree<T, EVALTYPE> & anothertree){
	assert(this->getGeneralTree().root->getContent()->getIndex().length() == anothertree.getGeneralTree().root->getContent()->getIndex().length());
	GeneralTree<HTuckerTreeNode<T, EVALTYPE> > tmp = anothertree.getGeneralTree();
	GeneralTreeNode<HTuckerTreeNode<T, EVALTYPE> > *nodethis, *nodeanother, *nodetmp;
	for(GeneralTreeIterator<HTuckerTreeNode<T,EVALTYPE> > TITthis = tree.end(),
		GeneralTreeIterator<HTuckerTreeNode<T,EVALTYPE> > TITanother = anothertree.getGeneralTree().end(),
		GeneralTreeIterator<HTuckerTreeNode<T,EVALTYPE> > TITtmp = tmp.end(); TITthis >= tree.begin(); TITthis--,TITanother--,TITtmp--){
	
		nodethis = TITthis.getNode();
		nodeanother = TITanother.getNode();
		nodetmp = TITtmp.getNode();
		//std::cout << nodethis->getContent()->getIndex() << std::endl;
		if(nodethis->isLeaf()){
			//std::cout << " " << nodethis->getContent()->getUorB().numRows()  << "== " <<  nodeanother->getContent()->getUorB().numRows()  << std::endl;
			assert(nodethis->getContent()->getUorB().numRows() == nodeanother->getContent()->getUorB().numRows());
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmpM;
			
			tmpM = transpose(nodethis->getContent()->getUorB())*nodeanother->getContent()->getUorB();
			nodetmp->getContent()->setUorB(tmpM);
		} else{
			// the resalt has to be saved in non-blockwise format, because we can than multiply more easily
			// i.e. at this point we can assume that nodeanother and notethis have blockwise format and nodetemp non-blockwise format
			int numrows_left = nodetmp->getfirstChild()->getContent()->getUorB().numRows();
			int numrows_right = nodetmp->getlastChild()->getContent()->getUorB().numRows();
			int numcols_left = nodetmp->getfirstChild()->getContent()->getUorB().numCols();
			int numcols_right = nodetmp->getlastChild()->getContent()->getUorB().numCols();
			int anothernumcols = nodeanother->getContent()->getUorB().numCols();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tensorresult(nodethis->getContent()->UorB_numel,nodeanother->getContent()->UorB_numel);
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > intermediateresult,intermediateresult2;
			for(int i = 1; i<= nodeanother->getContent()->UorB_numel;++i){
				intermediateresult = nodetmp->getlastChild()->getContent()->getUorB() * (nodeanother->getContent()->getUorB())(_((i-1)*nodeanother->getContent()->UorB_rcnumel + 1 ,i * nodeanother->getContent()->UorB_rcnumel),_);
				intermediateresult2 = intermediateresult * transpose(nodetmp->getfirstChild()->getContent()->getUorB());
				flens::DenseVector<flens::Array<T> > rhresult(intermediateresult2.numRows() * intermediateresult2.numCols());
				for(int j = 1; j <= intermediateresult.numCols(); ++j){
					rhresult(_((j-1)*intermediateresult2.numRows() + 1, j*intermediateresult2.numRows())) = intermediateresult2(_,j);
				}

				for(int j = 1; j<= nodethis->getContent()->UorB_numel; ++j){
					tensorresult(j,i) = 0;
					for(int k = 1; k <= rhresult.length(); ++k){
						tensorresult(j,i) += rhresult(k)*(nodethis->getContent()->getUorB())((j-1)*nodethis->getContent()->UorB_numel + ((k-1) % nodethis->getContent()->UorB_numel) + 1, (int)(((k-1) - (k-1) % nodethis->getContent()->UorB_numel)/nodethis->getContent()->UorB_numel) + 1);
					}
				}
			}
			
			nodetmp->getContent()->setUorB(tensorresult);
		}
	}
	return (tmp.root->getContent()->getUorB())(1,1);
};


template <typename T, class EVALTYPE>
void 
HTuckerTree<T,EVALTYPE>::generateTofElementary(DenseVectorList<T> list,int k,int d){
	GeneralTreeNode<HTuckerTreeNode<T,EVALTYPE> > *node;
	for(GeneralTreeIterator<HTuckerTreeNode<T,EVALTYPE> > TIT = tree.end(); TIT >= tree.begin(); TIT--){
		node = TIT.getNode();
		if(node->isLeaf()){
			int pos = (node->getContent()->getIndex())[0];
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > mat3((*list((pos-1)*k+1)).length(),k);
			for(int i = 1; i<= k; ++i){
				mat3(_,i) = *list((pos-1)*k+i);
			}
			node->getContent()->setUorB(mat3);
		} else if(node->isRoot()) {
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > mat(k,k);
			for(int i = 1; i<=k; ++i){
				mat(i,i) = 1.0;
			}
			node->getContent()->setUorB(mat);
			node->getContent()->UorB_numel = 1;
			node->getContent()->UorB_rcnumel = k;
			node->getContent()->UorB_lcnumel = k;
		} else {
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > mat2(k*k,k);
			for(int i = 1; i<= k; ++i){
				mat2((i-1)*k+i,i) = 1.0;
			}
			node->getContent()->setUorB(mat2);
			node->getContent()->UorB_numel = k;
			node->getContent()->UorB_rcnumel = k;
			node->getContent()->UorB_lcnumel = k;
		}
	}

};

template <typename T, class EVALTYPE>
T
HTuckerTree<T,EVALTYPE>::Linfnorm(HTuckerTree<T,EVALTYPE > & anothertree,DimensionIndex minval,DimensionIndex maxval,int n){
	T maxi = -1;
	int r;
	T res;
	DimensionIndex eval(minval.length());
	for(int i = 1; i<= n; ++i){
		std::cout << "i = " << i << std::endl;
		eval.setRandom(minval,maxval);
		
		res = abs(evaluate(eval) -anothertree.evaluate(eval));
		if(res > maxi){
			maxi = res;
		}
	}
	return maxi;
};

template <typename T, class EVALTYPE>
T
HTuckerTree<T,EVALTYPE>::Linfnorm(DimensionIndex minval,DimensionIndex maxval,int n){
	T maxi = -1;
	int r;
	T res;
	DimensionIndex eval(minval.length());
	for(int i = 1; i<= n; ++i){
		
		eval.setRandom(minval,maxval);
		
		res = abs(evaluate(eval));
		if(res > maxi){
			maxi = res;
		}
	}
	return maxi;

};

template <typename T, class EVALTYPE>
T
HTuckerTree<T,EVALTYPE>::Linfnorm(tensor<T> tens,DimensionIndex minval, DimensionIndex maxval,int n){
	T maxi = -1;
	int r;
	T res;
	stopwatch s1;
	DimensionIndex eval(minval.length());
	s1.start();
	for(int i = 1; i<= n; ++i){
		eval.setRandom(minval,maxval);
		if(i % 1000 == 0){
			s1.stop();
			std::cout << "Linfnorm " << i << "  time = " << s1.getTime() << std::endl;
			s1.start();
		}
		res = abs(evaluate(eval) - tens(eval));
		if(res > maxi){
			maxi = res;
		}
	}
	s1.stop();
	return maxi;
};

template <typename T, class EVALTYPE>
T
HTuckerTree<T,EVALTYPE>::Linfnorm(funpoint fun, DimensionIndex minval, DimensionIndex maxval,int n){
	T maxi = -1;
	int r;
	T res;
	DimensionIndex eval(minval.length());
	for(int i = 1; i<= n; ++i){
		eval.setRandom(minval,maxval);
		
		res = abs(evaluate(eval) - fun(eval));
		if(res > maxi){
			maxi = res;
		}
	}
	return maxi;
};


template <typename T>
void 
setUorB(HTuckerTree<T,SVD > & tree1){
	typedef double T;
	double eps = 0.00000000001;
	int d = tree1.getGeneralTree().root->getContent()->getIndex().length();
	GeneralTreeNode<HTuckerTreeNode<double,SVD> > *node;
	for(GeneralTreeIterator<HTuckerTreeNode<double,SVD> > TIT = tree1.getGeneralTree().end(); TIT >= tree1.getGeneralTree().begin(); TIT--){
		
		node = TIT.getNode();
		//std::cout << "htuckertree.tcc setUorB(HTuckerTree<T,SVD > & tree1), " << node->getContent()->getIndex() << std::endl;
		if(node->isLeaf()){
			int minval = node->getContent()->getROA().A.minval[(node->getContent()->getIndex())[0]-1];
			int maxval = node->getContent()->getROA().A.maxval[(node->getContent()->getIndex())[0]-1];
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > U(maxval - minval + 1,node->getContent()->getROA().pivots.length());
			for(int i = 0; i < node->getContent()->getROA().pivots.length(); ++i){
				U(_,i) = node->getContent()->getROA().evaluate(node->getContent()->getROA().pivots[i],node->getContent()->getIndex(),(node->getContent()->getIndex())[0],minval,maxval);
			}
			flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > Utilde;
			Utilde = U*transpose(node->getContent()->getROA().V);
			node->getContent()->setUorB(Utilde);
		} else {
			int lcnumel = node->getfirstChild()->getContent()->getROA().pivots.length();
			int rcnumel = node->getlastChild()->getContent()->getROA().pivots.length();
			int numel;
			node->getContent()->UorB_numel = numel;
			node->getContent()->UorB_lcnumel = lcnumel;
			node->getContent()->UorB_rcnumel = rcnumel;
			if(node->isRoot()){
				numel = 1;
			} else {
				numel = node->getContent()->getROA().pivots.length();
			}

			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Tinner(numel * rcnumel,lcnumel);
			DimensionIndex idx(tree1.getGeneralTree().root->getContent()->getIndex().length());
			for(int i = 0; i < numel; ++i){
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Amat(rcnumel,lcnumel);
				if(!(node->isRoot())){
					idx = (node->getContent()->getROA().pivots)[i];
				} 
				for(int j = 0; j < lcnumel; ++j){
					for(int k = 0; k < rcnumel; ++k){
						idx.setValue(node->getfirstChild()->getContent()->getROA().pivots[j],node->getfirstChild()->getContent()->getIndex());
						idx.setValue(node->getlastChild()->getContent()->getROA().pivots[k],node->getlastChild()->getContent()->getIndex());
						Amat(k+1,j+1) = (node->getContent()->getROA())(idx);
					}
				}
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > AU, AVT;
				flens::DenseVector<flens::Array<T> > As;
				svd(Amat,As,AU,AVT);
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > O1,O2;
				O1 = transpose(node->getlastChild()->getContent()->getROA().U)*AU;
				O2 = AVT*node->getfirstChild()->getContent()->getROA().U;
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > o1help(O1.numRows(),O1.numCols());
				for(int j = 1; j <= O1.numRows(); ++j){
					o1help(j,_) = O1(j,_)*As(j);
				}
				O1 = o1help * O2;
				O2 = O1;
				flens::DenseVector<flens::Array<T > > sl,sr;
				sl = node->getfirstChild()->getContent()->getROA().s;
				sr = node->getlastChild()->getContent()->getROA().s;
				for(int j = 1; j<= O2.numRows(); ++j){
					for(int k = 1; k<= O2.numCols(); ++k){
						if(abs(sl(k)) > eps && abs(sr(j)) > eps){
							O2(j,k) = O1(j,k) / sl(k)/sr(k);
						} else {
							O2(j,k) = 0;
						}
					}
				}
				Tinner(_(i*rcnumel+1, (i+1)*rcnumel),_) = O2;
			} //end for(i = 0; i < numel; ++i){
			
			if(!(node->isRoot())){
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > V;
				V = transpose(node->getContent()->getROA().V);
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Ttilde(V.numCols() * rcnumel,lcnumel);
				for(int i = 1; i<= V.numCols() * rcnumel; ++i){
					for(int j = 1; j<= lcnumel; ++j){
						Ttilde(i,j) = 0;
						for(int k = 1; k <= numel; ++k){
							Ttilde(i,j) += Tinner( ((i-1) % rcnumel) + 1 + (k-1)*rcnumel,j) * V(k, (int)((i-1) - ((i-1) % rcnumel)) + 1);
						} 
					}
				}
				node->getContent()->setUorB(Ttilde);
			} else {
				node->getContent()->setUorB(Tinner);
			}

			
		
		} //end else {
	} // end for
};

template <typename T, class EVALTYPE>
HTuckerTree<T,EVALTYPE >  operator+(HTuckerTree<T,EVALTYPE>  & tree1, HTuckerTree<T,EVALTYPE > & tree2){
	assert(tree1.getGeneralTree().root->getContent()->getIndex().length() == tree2.getGeneralTree().root->getContent()->getIndex().length());
	HTuckerTree<T,EVALTYPE>  tmp = tree1;
	GeneralTreeNode<HTuckerTreeNode<T,EVALTYPE> > *node1,*node2,*nodetmp;
	int kt1,kt2,kt1s,kt2s,kt1ss,kt2ss,kt,kts,ktss;
	for(GeneralTreeIterator<HTuckerTreeNode<T,EVALTYPE> > TIT1 = tree1.getGeneralTree().end(),
		GeneralTreeIterator<HTuckerTreeNode<T,EVALTYPE> > TIT2 = tree2.getGeneralTree().end(),
		GeneralTreeIterator<HTuckerTreeNode<T,EVALTYPE> > TITtmp = tmp.getGeneralTree().end(); TIT1 >= tree1.getGeneralTree().begin(); TIT1--,TIT2--,TITtmp--){
		
		node1 = TIT1.getNode();
		node2 = TIT2.getNode();
		nodetmp = TITtmp.getNode();
		//std::cout << "htuckertree.tcc operator+ " << node1->getContent()->getIndex() << " " << node2->getContent()->getIndex() << " " << nodetmp->getContent()->getIndex() << std::endl;
		if(TIT1.getNode()->isLeaf()){
			//std::cout << node1->getContent()->getUorB() << "  " << node2->getContent()->getUorB() << std::endl;
			flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > tmpM(max(node1->getContent()->getUorB().numRows(),node2->getContent()->getUorB().numRows()),node1->getContent()->getUorB().numCols()+node2->getContent()->getUorB().numCols());
			tmpM(_(1,node1->getContent()->getUorB().numRows()),_(1,node1->getContent()->getUorB().numCols())) = node1->getContent()->getUorB();
			tmpM(_(1,node2->getContent()->getUorB().numRows()),_(node1->getContent()->getUorB().numCols()+1,node1->getContent()->getUorB().numCols()+ node2->getContent()->getUorB().numCols())) = node2->getContent()->getUorB();
			nodetmp->getContent()->setUorB(tmpM);
		} else if(TIT1.getNode()->isRoot()){
			int rcnumel1 = node1->getContent()->UorB_rcnumel;
			int rcnumel2 = node2->getContent()->UorB_rcnumel;
			int lcnumel1 = node1->getContent()->UorB_lcnumel;
			int lcnumel2 = node2->getContent()->UorB_lcnumel;
			// the variables numel1 and numel2 should be = 1;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > UorBsum(rcnumel1+rcnumel2,lcnumel1+lcnumel2);
			UorBsum(_(1,rcnumel1),_(1,lcnumel1)) = node1->getContent()->getUorB();
			UorBsum(_(rcnumel1+1, rcnumel1+rcnumel2),_(lcnumel1+1,lcnumel1+lcnumel2)) = node2->getContent()->getUorB();
			nodetmp->getContent()->setUorB(UorBsum);
			nodetmp->getContent()->UorB_numel = 1;
			nodetmp->getContent()->UorB_lcnumel = lcnumel1 + lcnumel2;
			nodetmp->getContent()->UorB_rcnumel = rcnumel1 + rcnumel2;
			
		} else {
			int rcnumel1 = node1->getContent()->UorB_rcnumel;
			int rcnumel2 = node2->getContent()->UorB_rcnumel;
			int lcnumel1 = node1->getContent()->UorB_lcnumel;
			int lcnumel2 = node2->getContent()->UorB_lcnumel;
			int numel1 = node1->getContent()->UorB_numel;
			int numel2 = node2->getContent()->UorB_numel;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > UorBsum((rcnumel1+rcnumel2)*(numel1+numel2),lcnumel1 + lcnumel2);
			for(int i = 1; i <= numel1; ++i){
				UorBsum(_((i-1)*(numel1 + numel2) + 1,(i-1)*(numel1 + numel2) + numel1),_(1,lcnumel1)) = (node1->getContent()->getUorB())(_((i-1)*rcnumel1+1,i*rcnumel1),_);
			}
			for(int i = 1; i<= numel2; ++i){
				UorBsum(_(numel1*(numel1 + numel2) + (i-1) * (numel1 + numel2) + numel1 + 1,numel1*(numel1 + numel2) + (i-1) * (numel1 + numel2)),_(lcnumel1+1,lcnumel1 + lcnumel2)) = (node2->getContent()->getUorB())(_((i-1)*rcnumel2+1,i*rcnumel2),_);
			}
			nodetmp->getContent()->setUorB(UorBsum);
			nodetmp->getContent()->UorB_numel = numel1 + numel2;
			nodetmp->getContent()->UorB_lcnumel = lcnumel1 + lcnumel2;
			nodetmp->getContent()->UorB_rcnumel = rcnumel1 + rcnumel2;
		}

	}

	return tmp;
};

template <typename T, class EVALTYPE>
T 
HTuckerTree<T,EVALTYPE>::average_rank(){
	int count = 0;
	T rank = 0.0;
	for(GeneralTreeIterator<HTuckerTreeNode<T,EVALTYPE> > TIT = tree.end(); TIT >= tree.begin(); TIT--){
		if(TIT.getNode()->isLeaf()){
			rank += TIT.getNode()->getContent()->getUorB().numCols();
		} else  {
			rank += TIT.getNode()->getContent()->UorB_numel;
		}
		count++;
	}

	return rank/(count +0.0);
};

template <typename T, class EVALTYPE>
T 
HTuckerTree<T,EVALTYPE>::max_rank(){
	T rank = 0.0;
	for(GeneralTreeIterator<HTuckerTreeNode<T,EVALTYPE> > TIT = tree.end(); TIT >= tree.begin(); TIT--){
		int newrank; 
		if(TIT.getNode()->isLeaf()){
			 newrank = TIT.getNode()->getContent()->getUorB().numCols();
		} else  {
			 newrank = TIT.getNode()->getContent()->UorB_numel;
		}
		if(newrank > rank){
			rank = newrank;
		}
	}

	return rank +0.0;
};

template <typename T, class EVALTYPE>
T 
HTuckerTree<T,EVALTYPE>::effective_rank(){
	int storage = 0;
	int n = 0;
	int count = 0;
	for(GeneralTreeIterator<HTuckerTreeNode<T,EVALTYPE> > TIT = tree.end(); TIT >= tree.begin(); TIT--){
		if(TIT.getNode()->isLeaf()){
			storage = storage + TIT.getNode()->getContent()->getUorB().numRows()*TIT.getNode()->getContent()->getUorB().numCols();
			n = n + TIT.getNode()->getContent()->getUorB().numRows();
			count ++;
		} else {
			storage = storage + TIT.getNode()->getContent()->UorB_numel*TIT.getNode()->getContent()->UorB_lcnumel*TIT.getNode()->getContent()->UorB_rcnumel;
		}
	}
	
	double n_av = n/(count + 0.0);
	double eps = 0.0001;
	double k = 0;
	double K = 10;
	double middle;
	while(abs((d-1)*k*k*k + d * n_av * k - storage) > eps){
		if((d-1)*K*K*K + d * n_av * K - storage < 0){
			K *= 2;
		} else {
			middle = (K+k)/2;
			if((d-1)*middle*middle*middle + d * n_av * middle - storage < 0){
				k = middle;
			} else {
				K = middle;
			}
		}
	}
	return T(k);
};

template <typename T, class EVALTYPE>
HTuckerTree<T,EVALTYPE>   operator*(double number, HTuckerTree<T,EVALTYPE>  & tree){
	HTuckerTree<T,EVALTYPE>  tmp = tree;
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmpmat;
	tmpmat = number * tmp.getGeneralTree().root->getContent()->getUorB();
	tmp.getGeneralTree().root->getContent()->setUorB(tmpmat);
	return tmp;
};

template <typename T, class EVALTYPE>
HTuckerTree<T,EVALTYPE>   operator*(HTuckerTree<T,EVALTYPE> & tree,double number){
	return number * tree;
};

template <typename T, class EVALTYPE>
HTuckerTree<T,EVALTYPE>   operator-(HTuckerTree<T,EVALTYPE>  & tree1, HTuckerTree<T,EVALTYPE>  & tree2){
	return tree1 + (-1)*tree2;
};
	


