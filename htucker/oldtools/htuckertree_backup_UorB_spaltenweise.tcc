
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
			
			int numcols_left = node->getfirstChild()->getContent()->getUorB().numCols();
			int numcols_right = node->getlastChild()->getContent()->getUorB().numCols();
			flens::DenseVector<flens::Array <T> > bpart(numcols_left*numcols_right);
			flens::DenseVector<flens::Array <T> > intermediateresult;
			
			bpart = (node->getContent()->getUorB())(_,1);
		
			intermediateresult = Vec2Mat(bpart,numcols_right,numcols_left)*node->getfirstChild()->getContent()->getEvaluate();
			evaluate(1) = 0;
			flens::DenseVector<flens::Array<T> > evaluate_right = node->getlastChild()->getContent()->getEvaluate();
			for(int j = 1; j<= numcols_right; ++j){
				evaluate(1) += intermediateresult(j)*evaluate_right(j);
			}
			return evaluate(1);		
		} else if(node->isLeaf()){
			assert(node->getContent()->getIndex().length() == 1);
			int numcols = node->getContent()->getUorB().numCols();
			evaluate = flens::DenseVector<flens::Array<T> >(numcols);
			int val = index[(node->getContent()->getIndex())[0] - 1];
			evaluate = (node->getContent()->getUorB())(val,_);
			node->getContent()->setEvaluate(evaluate);
		} else{
			int numcols = node->getContent()->getUorB().numCols();
			evaluate = flens::DenseVector<flens::Array<T> >(numcols);
			
			int numcols_left = node->getfirstChild()->getContent()->getUorB().numCols();
			int numcols_right = node->getlastChild()->getContent()->getUorB().numCols();

			flens::DenseVector<flens::Array <T> > bpart(numcols_left*numcols_right);
			flens::DenseVector<flens::Array <T> > intermediateresult;
			
			for(int i= 1; i<= numcols; ++i){
				for(int j = 1; j<= numcols_left*numcols_right ; j++){
					bpart(j) = (node->getContent()->getUorB())(j,i);
				}
				
				intermediateresult = Vec2Mat(bpart,numcols_right,numcols_left)*node->getfirstChild()->getContent()->getEvaluate();
				evaluate(i) = 0;
				for(int j = 1; j<= numcols_right; ++j){
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
			//compute QR-decomposition of left child
			//std::cout << "QR of left child " << std::endl;
			int left_numrows = node->getfirstChild()->getContent()->getUorB().numRows();
			int left_numcols = node->getfirstChild()->getContent()->getUorB().numCols();
			//std::cout << "left_numrows = " << left_numrows<< "   left_numcols =  " << left_numcols << std::endl;
			flens::DenseVector<flens::Array<double> > taul(min(left_numrows,left_numcols));
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Rl;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Ql;
			Rl = node->getfirstChild()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Rlneu(min(Rl.numRows(),Rl.numCols()),Rl.numCols());
			//std::cout << "qrf " << std::endl;
			qrf(Rl,taul);
			//std::cout << "taul = " << taul << std::endl;
			if(left_numcols > left_numrows){
				Ql = Rl(_(1,left_numrows),_(1,left_numrows));
			} else {
				Ql = Rl;
			}
			
			//std::cout << "orgqr " << std::endl;
			orgqr(Ql,taul);
			//std::cout << "taul " << taul << std::endl;
			
			for(int i = Rlneu.firstRow(); i <= Rlneu.lastRow(); ++i){
				for(int j = Rlneu.firstCol(); j <= Rlneu.lastCol(); ++j){
					if(i>j){
						Rlneu(i,j) = 0;
					} else {
						Rlneu(i,j) = Rl(i,j);
					}
				}
			}
				
			//compute QR-decomposition of right child
			//std::cout << "QR of right child " << std::endl;
			int right_numrows = node->getlastChild()->getContent()->getUorB().numRows();
			int right_numcols = node->getlastChild()->getContent()->getUorB().numCols();
			flens::DenseVector<flens::Array<double> > taur(min(right_numrows,right_numcols));
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Rr = node->getlastChild()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Qr;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Rrneu(min(Rr.numRows(),Rr.numCols()),Rr.numCols());
			qrf(Rr,taur);
			if(right_numcols > right_numrows){
				Qr = Rr(_(1,right_numrows),_(1,right_numrows));
			} else {
				Qr = Rr;
			}
			orgqr(Qr,taur);
			for(int i = Rrneu.firstRow(); i <= Rrneu.lastRow(); ++i){
				for(int j = Rrneu.firstCol(); j <= Rrneu.lastCol(); ++j){
					if(i>j){
						Rrneu(i,j) = 0;
					} else {
						Rrneu(i,j) = Rr(i,j);
					}
				}
			}
				

			//std::cout << " Set Q matrices " << std::endl;
			//set orthogonal Matrices to childnodes
			node->getfirstChild()->getContent()->setUorB(Ql);
			node->getlastChild()->getContent()->setUorB(Qr);

			//compute new transfer Tensor:
			//std::cout << "compute new transfer tensor " << std::endl;
			int numcols = node->getContent()->getUorB().numCols();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Bneu(Rlneu.numRows()*Rrneu.numRows(),numcols);
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > intermediateresult,endresult;
			
			for(int i = node->getContent()->getUorB().firstCol(); i<= node->getContent()->getUorB().lastCol(); ++i){
				intermediateresult = Vec2Mat((node->getContent()->getUorB())(_,i),Rrneu.numCols(),Rlneu.numCols())*transpose(Rlneu);
				endresult = Rrneu*intermediateresult;
				Bneu(_,i) = Mat2Vec(endresult);
			}
			node->getContent()->setUorB(Bneu);
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
		erg += ( tree.root->getContent()->getUorB())(i,1)*( tree.root->getContent()->getUorB())(i,1);
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
			int numrows_left = nodetmp->getfirstChild()->getContent()->getUorB().numRows();
			int numrows_right = nodetmp->getlastChild()->getContent()->getUorB().numRows();
			int numcols_left = nodetmp->getfirstChild()->getContent()->getUorB().numCols();
			int numcols_right = nodetmp->getlastChild()->getContent()->getUorB().numCols();
			int anothernumcols = nodeanother->getContent()->getUorB().numCols();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tensorresult(numrows_left * numrows_right,anothernumcols);
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > intermediateresult,intermediateresult2;
			
			for(int i = nodeanother->getContent()->getUorB().firstCol(); i <= nodeanother->getContent()->getUorB().lastCol(); ++i){
				intermediateresult = Vec2Mat((nodeanother->getContent()->getUorB())(_,i),numcols_right,numcols_left)*transpose(nodetmp->getfirstChild()->getContent()->getUorB());
				intermediateresult2 = nodetmp->getlastChild()->getContent()->getUorB() * intermediateresult;
				tensorresult(_,i) = Mat2Vec(intermediateresult2);
			}
			
			intermediateresult = transpose(nodethis->getContent()->getUorB()) * tensorresult;
			nodetmp->getContent()->setUorB(intermediateresult);
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
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > mat(k*k,1);
			for(int i = 1; i<=k; ++i){
				mat((i-1)*k+i,1) = 1.0;
			}
			node->getContent()->setUorB(mat);
		} else {
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > mat2(k*k,k);
			for(int i = 1; i<= k; ++i){
				mat2((i-1)*k+i,i) = 1.0;
			}
			node->getContent()->setUorB(mat2);
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
setUorB(HTuckerTree<T,INVERSE > & tree1){
	typedef double T;
	int d = tree1.getGeneralTree().root->getContent()->getIndex().length();
	GeneralTreeNode<HTuckerTreeNode<double,INVERSE> > *node;
	for(GeneralTreeIterator<HTuckerTreeNode<double,INVERSE> > TIT = tree1.getGeneralTree().end(); TIT >= tree1.getGeneralTree().begin(); TIT--){
		node = TIT.getNode();
		if(TIT.getNode()->isLeaf()){
			//here we should have only one dimension in the t-index, otherwise it is not a leaf.
			assert(node->getContent()->getIndex().length() == 1);
			int maxval = node->getContent()->getROA().A.maxval[(node->getContent()->getIndex())[0] -1 ];
			int minval = node->getContent()->getROA().A.minval[(node->getContent()->getIndex())[0] -1 ]; 
			
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > UorBnew(maxval - minval + 1,node->getContent()->getROA().pivots.length() ,1,1);
			DimensionIndex idx(d);
			for(int j = 0; j <= node->getContent()->getROA().pivots.length(); ++j){
				idx = (node->getContent()->getROA().pivots)[j];
				for(DimensionIndexIterator<Iterator1D> it = idx.getIterator((node->getContent()->getIndex())[0],minval,maxval),int i = 1; it.inRange(); it++,i++){
					UorBnew(i,j+1) = (node->getContent()->getROA())(idx);
				}
			}
			node->getContent()->setUorB(UorBnew);
		} else if(TIT.getNode()->isRoot()) {
				
			flens::DenseVector<flens::Array<T> > tensorprod;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tensorres;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tensorres2;
			DimensionIndex idx;

			int right_pivots_length  = node->getlastChild()->getContent()->getROA().pivots.length();
			int left_pivots_length  = node->getlastChild()->getContent()->getROA().pivots.length();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > UorBnew(left_pivots_length * right_pivots_length,1,1,1);
			
			flens::DenseVector<flens::Array<T> > input(left_pivots_length * right_pivots_length);
			for(int j = 0; j < left_pivots_length ; ++j){
				idx = (node->getfirstChild()->getContent()->getROA().pivots)[j];
				for(int k = 0; k < right_pivots_length; ++k){
					idx.setValue(node->getlastChild()->getContent()->getROA().pivots[k],node->getlastChild()->getContent()->getIndex());
					input(j * right_pivots_length + k + 1) = node->getContent()->getROA().A(idx);
				}
			}
			
			tensorres = node->getlastChild()->getContent()->getROA().Ainverse * Vec2Mat(input,right_pivots_length,left_pivots_length);
			tensorres2 = tensorres * transpose(node->getfirstChild()->getContent()->getROA().Ainverse);
			tensorprod = Mat2Vec(tensorres);
			
			for(int j = tensorprod.firstIndex(); j <= tensorprod.lastIndex(); ++j){
				UorBnew(j,1) = tensorprod(j);  // Spaltenweise eintragen... fertig!
			}
			
			node->getContent()->setUorB(UorBnew);
		} else {
				
			flens::DenseVector<flens::Array<T> > tensorprod;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tensorres;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tensorres2;

			int right_pivots_length  = node->getlastChild()->getContent()->getROA().pivots.length();
			int left_pivots_length  = node->getlastChild()->getContent()->getROA().pivots.length();
			int pivots_length = node->getContent()->getROA().pivots.length();
			
			DimensionIndex idx;

			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > UorBnew(left_pivots_length * right_pivots_length,pivots_length,1,1);
			for(int i = 0; i < pivots_length; ++i){
				flens::DenseVector<flens::Array<T> > input(left_pivots_length * right_pivots_length);
				idx = (node->getContent()->getROA().pivots)[i];
				for(int j = 0; j < left_pivots_length; ++j){
					for(int k = 0; k < right_pivots_length; ++k){
						idx.setValue(node->getfirstChild()->getContent()->getROA().pivots[j],node->getfirstChild()->getContent()->getIndex());
						idx.setValue(node->getlastChild()->getContent()->getROA().pivots[k],node->getlastChild()->getContent()->getIndex());
					
						input(j*right_pivots_length + k + 1) = (node->getContent()->getROA())(idx);
					}
				}
				
				tensorres = node->getlastChild()->getContent()->getROA().Ainverse * Vec2Mat(input,right_pivots_length,left_pivots_length);
				tensorres2 = transpose(tensorres);
				tensorres = tensorres2 * transpose(node->getlastChild()->getContent()->getROA().Ainverse);

				tensorprod = Mat2Vec(tensorres);
				for(int j = tensorprod.firstIndex(); j <= tensorprod.lastIndex(); ++j){
					UorBnew(j,i + 1) = tensorprod(j);  // Spaltenweise eintragen... fertig!
				}

			}
				
			node->getContent()->setUorB(UorBnew);	
		}
	}
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
		if(TIT.getNode()->isLeaf()){
			//here we should have only one dimension in the t-index, otherwise it is not a leaf.
			assert(node->getContent()->getIndex().length() == 1);
			int maxval = node->getContent()->getROA().A.maxval[(node->getContent()->getIndex())[0] -1 ];
			int minval = node->getContent()->getROA().A.minval[(node->getContent()->getIndex())[0] -1 ]; 
			
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > UorBnew(maxval - minval + 1,node->getContent()->getROA().pivots.length() ,1,1);
			DimensionIndex idx(d);
			for(int j = 0; j < node->getContent()->getROA().pivots.length(); ++j){
				idx = (node->getContent()->getROA().pivots)[j];
				
				
				for(DimensionIndexIterator<Iterator1D> it = idx.getIterator((node->getContent()->getIndex())[0],minval,maxval),int i = 1; it.inRange(); it++,i++){
					//std::cout << "htuckertree.tcc setUorB, i = " << i << "  idx = " << it.getIndex() << std::endl;
					UorBnew(i,j+1) = (node->getContent()->getROA())(it.getIndex());
				}
			}
			//std::cout << "htuckertree.tcc, setUorB(HTuckerTree<T,SVD > & tree1), " << UorBnew << std::endl;
			node->getContent()->setUorB(UorBnew);
		} else if(TIT.getNode()->isRoot()) {
				
			flens::DenseVector<flens::Array<T> > tensorprod;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tensorres;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tensorres2;
			DimensionIndex idx;

			int right_pivots_length  = node->getlastChild()->getContent()->getROA().pivots.length();
			int left_pivots_length  = node->getfirstChild()->getContent()->getROA().pivots.length();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > UorBnew(left_pivots_length * right_pivots_length,1,1,1);
			
			flens::DenseVector<flens::Array<T> > input(left_pivots_length * right_pivots_length);
			for(int j = 0; j < left_pivots_length ; ++j){
				idx = (node->getfirstChild()->getContent()->getROA().pivots)[j];
				for(int k = 0; k < right_pivots_length; ++k){
					idx.setValue(node->getlastChild()->getContent()->getROA().pivots[k],node->getlastChild()->getContent()->getIndex());
					input(j * right_pivots_length + k + 1) = node->getContent()->getROA().A(idx);
				}
			}
			
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > newdiagr(node->getlastChild()->getContent()->getROA().s.length(),node->getlastChild()->getContent()->getROA().s.length());
			for(int j = node->getlastChild()->getContent()->getROA().s.firstIndex(); j<= node->getlastChild()->getContent()->getROA().s.lastIndex(); ++j){
				if(abs(node->getlastChild()->getContent()->getROA().s(j)) > eps){
					newdiagr(j,j) = 1/node->getlastChild()->getContent()->getROA().s(j);
				} else {
					newdiagr(j,j) = 0;
				}
			}

			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > newdiagl(node->getfirstChild()->getContent()->getROA().s.length(),node->getfirstChild()->getContent()->getROA().s.length());
			for(int j = node->getfirstChild()->getContent()->getROA().s.firstIndex(); j<= node->getfirstChild()->getContent()->getROA().s.lastIndex(); ++j){
				if(abs(node->getfirstChild()->getContent()->getROA().s(j)) > eps){
					newdiagl(j,j) = 1/node->getfirstChild()->getContent()->getROA().s(j);
				} else {
					newdiagl(j,j) = 0;
				}
			}
			
			tensorres2 = transpose(node->getlastChild()->getContent()->getROA().U) * Vec2Mat(input,right_pivots_length,left_pivots_length);
			tensorres = newdiagr * tensorres2;
			tensorres2 = transpose(node->getlastChild()->getContent()->getROA().V) * tensorres;
			tensorres = tensorres2 * node->getfirstChild()->getContent()->getROA().U;
			tensorres2 = tensorres * newdiagl;
			tensorres = tensorres2 * node->getfirstChild()->getContent()->getROA().V;
			
			tensorprod = Mat2Vec(tensorres);
			for(int j = tensorprod.firstIndex(); j <= tensorprod.lastIndex(); ++j){
				UorBnew(j,1) = tensorprod(j);  // Spaltenweise eintragen... fertig!
			}
			
			node->getContent()->setUorB(UorBnew);
		} else {
				
			flens::DenseVector<flens::Array<T> > tensorprod;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tensorres;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tensorres2;

			int right_pivots_length  = node->getlastChild()->getContent()->getROA().pivots.length();
			int left_pivots_length  = node->getfirstChild()->getContent()->getROA().pivots.length();
			int pivots_length = node->getContent()->getROA().pivots.length();
			
			//std::cout << "htuckertree.tcc: setUorB: right: " << right_pivots_length << " left : " << left_pivots_length << " this = " << pivots_length << std::endl;
			DimensionIndex idx;

			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > UorBnew(left_pivots_length * right_pivots_length,pivots_length,1,1);
			for(int i = 0; i < pivots_length; ++i){
				flens::DenseVector<flens::Array<T> > input(left_pivots_length * right_pivots_length);
				idx = (node->getContent()->getROA().pivots)[i];
				for(int j = 0; j < left_pivots_length; ++j){
					for(int k = 0; k < right_pivots_length; ++k){
						idx.setValue(node->getfirstChild()->getContent()->getROA().pivots[j],node->getfirstChild()->getContent()->getIndex());
						idx.setValue(node->getlastChild()->getContent()->getROA().pivots[k],node->getlastChild()->getContent()->getIndex());
					
						input(j*right_pivots_length + k + 1) = (node->getContent()->getROA())(idx);
					}
				}

				//std::cout << "htuckertree.tcc. setUorB  input = " << input << std::endl;
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > newdiagr(node->getlastChild()->getContent()->getROA().s.length(),node->getlastChild()->getContent()->getROA().s.length());
				for(int j = node->getlastChild()->getContent()->getROA().s.firstIndex(); j<= node->getlastChild()->getContent()->getROA().s.lastIndex(); ++j){
					if(abs(node->getlastChild()->getContent()->getROA().s(j)) > eps){
						newdiagr(j,j) = 1/node->getlastChild()->getContent()->getROA().s(j);
					} else {
						newdiagr(j,j) = 0;
					}
				}
				
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > newdiagl(node->getfirstChild()->getContent()->getROA().s.length(),node->getfirstChild()->getContent()->getROA().s.length());
				for(int j = node->getfirstChild()->getContent()->getROA().s.firstIndex(); j<= node->getfirstChild()->getContent()->getROA().s.lastIndex(); ++j){
					if(abs(node->getfirstChild()->getContent()->getROA().s(j)) > eps){
						newdiagl(j,j) = 1/node->getfirstChild()->getContent()->getROA().s(j);
					} else {
						newdiagl(j,j) = 0;
					}
				}
				
				

				tensorres2 = transpose(node->getlastChild()->getContent()->getROA().U) * Vec2Mat(input,right_pivots_length,left_pivots_length);
				tensorres = newdiagr * tensorres2;
				tensorres2 = transpose(node->getlastChild()->getContent()->getROA().V) * tensorres;
				tensorres = tensorres2 * node->getfirstChild()->getContent()->getROA().U ;
				tensorres2 = tensorres * newdiagl;
				tensorres = tensorres2 * node->getfirstChild()->getContent()->getROA().V ;
				
				tensorprod = Mat2Vec(tensorres);

				for(int j = tensorprod.firstIndex(); j <= tensorprod.lastIndex(); ++j){
					UorBnew(j,i + 1) = tensorprod(j);  // Spaltenweise eintragen... fertig!
				}

			}
				
			node->getContent()->setUorB(UorBnew);	
		}
	}
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
			kt1 = node1->getfirstChild()->getContent()->getUorB().numCols();
			kt2 = node1->getlastChild()->getContent()->getUorB().numCols();
			kt1s =	node2->getfirstChild()->getContent()->getUorB().numCols();
			kt2s = node2->getlastChild()->getContent()->getUorB().numCols();
			kt1ss = kt1+kt1s;
			kt2ss = kt2 + kt2s;
			kt = node1->getContent()->getUorB().numCols();
			kts = node2->getContent()->getUorB().numCols();
			ktss = kt + kts;
			flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > mat(kt1ss*kt2ss,1);
			//std::cout << "hier vor der ersten for" << std::endl;
			for(int i = 1; i <= kt1; ++i){
				//std::cout << "in +: root: mat(_(" << (i-1)*(kt2ss)+1 << ", " << (i-1)*(kt2ss) + kt2 << ") , " << 1 << ") = lnode1->rankoneapprox.UorB(_(" << (i-1)*kt2 + 1 << ", " <<  i*kt2 << "),1)" << std::endl;
				mat(_((i-1)*(kt2ss)+1,(i-1)*(kt2ss) + kt2),1) = (node1->getContent()->getUorB())(_((i-1)*kt2 + 1, i*kt2),1);
			}
			//std::cout << "hier vor der zweiten for" << std::endl;
			for(int i = 1; i<= kt1s; ++i){
				//std::cout << "in +: root: mat(_( " << kt1*kt2ss + (i-1)*kt2ss + kt2+1 << ", " << kt1*kt2ss + i*kt2ss << "),1) = lnode2->rankoneapprox.UorB(_(" << (i-1)*kt2s + 1 << ", " << i*kt2s << "),1);" << std::endl;
				mat(_(kt1*kt2ss + (i-1)*kt2ss + kt2+1,kt1*kt2ss + i*kt2ss),1) = (node2->getContent()->getUorB())(_((i-1)*kt2s + 1,i*kt2s),1);
			}
			nodetmp->getContent()->setUorB(mat);
		} else {
			kt1 = node1->getfirstChild()->getContent()->getUorB().numCols();
			kt2 = node1->getlastChild()->getContent()->getUorB().numCols();
			kt1s = node2->getfirstChild()->getContent()->getUorB().numCols();
			kt2s = node2->getlastChild()->getContent()->getUorB().numCols();
			kt1ss = kt1+kt1s;
			kt2ss = kt2 + kt2s;
			kt = node1->getContent()->getUorB().numCols();
			kts = node2->getContent()->getUorB().numCols();
			ktss = kt + kts;
			//std::cout << "kt1 = " << kt1 << "  kt2 = " << kt2 << "  kt1s = " << kt1s << "  kt2s = " << kt2s << "  kt1ss = " << kt1ss << "  kt2ss = " << kt2ss << " kt = " << kt << " kts = " << kts << "  ktss = " << ktss << std::endl;
			flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > tmpMat(kt1ss* kt2ss,ktss);
			for(int i = 1; i <= kt1; ++i){
				//std::cout << "in + : tmpMat(_(" << (i-1)*(kt2ss)+1<< " , " << (i-1)*(kt2ss) + kt2 << "),_(" << 1 << ", " << kt << ")) = lnode1->rankoneapprox.UorB(_(" << (i-1)*kt2 + 1 << ", " << i*kt2 << "),_);" << std::endl;
				tmpMat(_((i-1)*(kt2ss)+1,(i-1)*(kt2ss) + kt2),_(1,kt)) = (node1->getContent()->getUorB())(_((i-1)*kt2 + 1, i*kt2),_);
			}
			
			for(int i = 1; i<= kt1s; i++){
				//std::cout << node2->getContent()->getUorB().numRows() << "  " << node2->getContent()->getUorB().numCols() << std::endl;
				//std::cout << "in +2 : tmpMat(_(" << kt1*kt2ss + (i-1)*kt2ss + kt2 +1 << "," << kt1*kt2ss + i*kt2ss << "),_(" << kt+1 << "," << ktss << ")) = lnode2->rankoneapprox.UorB(_(" << (i-1)*kt2s + 1 << "," << i*kt2s << "),_)" << std::endl;
				tmpMat(_(kt1*kt2ss + (i-1)*kt2ss + kt2 + 1,kt1*kt2ss + i*kt2ss),_(kt+1,ktss)) = (node2->getContent()->getUorB())(_((i-1)*kt2s + 1,i*kt2s),_);
			}
			//std::cout << "vor set " << std::endl;
			nodetmp->getContent()->setUorB(tmpMat);
			//std::cout << "nach set " << std::endl;
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
		rank += TIT.getNode()->getContent()->getUorB().numCols();
		count++;
	}

	return rank/(count +0.0);
};

template <typename T, class EVALTYPE>
T 
HTuckerTree<T,EVALTYPE>::max_rank(){
	T rank = 0.0;
	for(GeneralTreeIterator<HTuckerTreeNode<T,EVALTYPE> > TIT = tree.end(); TIT >= tree.begin(); TIT--){
		int newrank =  TIT.getNode()->getContent()->getUorB().numCols();
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
		storage = storage + TIT.getNode()->getContent()->getUorB().numRows()*TIT.getNode()->getContent()->getUorB().numCols();
		if(TIT.getNode()->isLeaf()){
			n = n + TIT.getNode()->getContent()->getUorB().numRows();
			count ++;
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
	


