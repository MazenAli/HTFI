#include <flens/flens.cxx>
#include <vector>
#include <cstddef>
#include <cassert>
#include <cmath>
#include <chrono>

namespace htucker{

template <typename T>
void
HTuckerTree<T>::init(const int _d, const double _split){
	assert(0 <= _split && _split <= 1);
	
	DimensionIndex all(_d);
	all.setValueAsc();
	
	HTuckerTreeNode<T> rootnode(all);
	
    rootnode.setNumRows(0);
    rootnode.setLeftChildNumRows(0);
    rootnode.setRightChildNumRows(0);

	tree = GeneralTree<HTuckerTreeNode<T> >(rootnode);
	
	GeneralTreeNode<HTuckerTreeNode<T> > * node = tree.root;
	GeneralTreeNode<HTuckerTreeNode<T> > * nextlev = NULL;

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

			HTuckerTreeNode<T> fino(fi);
			HTuckerTreeNode<T> lano(la);

            fino.setNumRows(0);
            fino.setLeftChildNumRows(0);
            fino.setRightChildNumRows(0);

            lano.setNumRows(0);
            lano.setLeftChildNumRows(0);
            lano.setRightChildNumRows(0);

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
}



template <typename T>
HTuckerTree<T>::HTuckerTree(const int _d):d(_d){
	init(_d,0.5);
}

template <typename T>
HTuckerTree<T>::HTuckerTree(const int _d,const double _split):d(_d){
	init(_d,_split);
}


template <typename T>
void
HTuckerTree<T>::set_tree(const HTuckerTree<T>& X)
{
    typedef     HTuckerTreeNode<T>          NType;

    tree.root   = new GeneralTreeNode<NType>(X.getGeneralTree().root->
                                            getContent()->getIndex());
    d           = X.dim();
    GeneralTreeNode<NType> * newbase = tree.root;
    GeneralTreeNode<NType> * base = X.getGeneralTree().root;
    GeneralTreeNode<NType> * save;
    GeneralTreeNode<NType> * nextlev = NULL;
    GeneralTreeNode<NType> * newnextlev = NULL;

    while(1){
        save = base->getfirstChild();
        if(save != NULL){
            do{
                NType   temp(save->getContent()->getIndex());
                newbase->appendChild(temp);
                if(nextlev == NULL && newnextlev == NULL){
                    nextlev = save;
                    newnextlev = newbase->getfirstChild();
                }
                save = save->getnextSibling();
            } while(save != base->getfirstChild());
        }
        if(base->getlevelright() != NULL){
            base = base->getlevelright();
            newbase = newbase->getlevelright();
        } else if(nextlev != NULL && newnextlev != NULL){
            base = nextlev;
            newbase = newnextlev;
            nextlev = NULL;
            newnextlev = NULL;
        } else {
            break;
        }
    }
}


template <typename T>
template <typename NType>
GeneralTree<NType>
HTuckerTree<T>::copy() const{
	GeneralTree<NType> ret(tree);
	return ret;
}

template <typename T>
void 
HTuckerTree<T>::print() const{
	tree.print();
}


template <typename T>
void 
HTuckerTree<T>::print_info() const{
	tree.print();
	HTuckerTreeNode<T> * save;
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.begin(); TIT <= tree.end(); TIT++){
		save = TIT.getNode()->getContent();
		std::cout << "Node " << save->getIndex() << std::endl;
		std::cout << "UorB " << save->getUorB().numRows() << " x " << save->getUorB().numCols() << std::endl;
		std::cout << " numel = " << save->getNumRows() << "  lcnumel = " << save->getLeftChildNumRows() << "  rcnumel = " << save->getRightChildNumRows() << std::endl;
	}

}

template <typename T>
void
HTuckerTree<T>::print_w_UorB() const{
	tree.print();
	HTuckerTreeNode<T> * save;
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.begin(); TIT <= tree.end(); TIT++){
		save = TIT.getNode()->getContent();
		std::cout << "Node " << save->getIndex() << std::endl;
		std::cout << "UorB " << save->getUorB() << std::endl;
		std::cout << "UorB " << save->getUorB().numRows() << " x " << save->getUorB().numCols() << std::endl;
		std::cout << " numel = " << save->getNumRows() << "  lcnumel = " << save->getLeftChildNumRows() << "  rcnumel = " << save->getRightChildNumRows() << std::endl;
	}

}


template <typename T>
void
HTuckerTree<T>::print_svs(bool isorth)
{
    using flens::_;

    HTuckerTree<T> Stree;
    Stree.set_tree(*this);
    if (!isorth) {
       this->orthogonalize();
    }

    HTuckerTree<T> gram = gramians_orthogonal(*this);
    flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > U;
    flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > Vt;
    flens::DenseVector<flens::Array<T> > s;

    GeneralTreeIterator<HTuckerTreeNode<T> > TITgram = gram.getGeneralTree().begin();
    std::vector<flens::DenseVector<flens::Array<T> > >  svals;
    int count = 0;
    for(; TITgram <= gram.getGeneralTree().end(); TITgram ++){
        flens::svd(TITgram.getNode()->getContent()->getUorB(),s,U,Vt);
        svals.push_back(s);
        count += s.length();
    }

    flens::DenseVector<flens::Array<T> > s_all(count);
    count = 1;
    for (std::size_t i=0; i<svals.size(); ++i) {
        s_all(_(count,count+svals[i].length()-1)) = svals[i];
        count += svals[i].length();
    }


    flens::DenseVector<flens::Array<T> > rho(s_all.length());
    flens::sort(s_all, rho);
    std::cout << "***The singular vals are***\n" << s_all
              << std::endl;
}


template <typename T>
void 
HTuckerTree<T>::print_values() const{
	std::cout << "d = " << d << std::endl;
	
	DimensionIndex start(1,d);
	DimensionIndex _min(1,d);
	DimensionIndex _max = this->getmax();
	int * dims = new int[d];
	for(int i = 0; i < d; ++i){
		dims[i] = i+1;
	}
	DimensionIndex activeDims(dims,d);
	//std::cout << "_max = " << _max << std::endl;
	std::cout << "[";
	for(DimensionIndexIterator<IteratorXD> iter = start.getIterator(activeDims,_min,_max); iter.inRange(); iter++){
		if(iter.getIndex().equals(_max,activeDims)){
		    std::cout << this->evaluate(iter.getIndex()) << "]" << std::endl;
		} else {
			std::cout << this->evaluate(iter.getIndex()) << ", ";
		}
	}
}

template <typename T>
const GeneralTree<HTuckerTreeNode<T> > &
HTuckerTree<T>::getGeneralTree() const{
	return tree;
}


template <typename T>
GeneralTree<HTuckerTreeNode<T> > &
HTuckerTree<T>::getGeneralTree(){
	return tree;
}


template <typename T>
T
HTuckerTree<T>::evaluate(const DimensionIndex &index) const{
    using flens::_;

	flens::DenseVector<flens::Array<T> > * evalarrays = new flens::DenseVector<flens::Array<T> >[d];
	GeneralTreeNode<HTuckerTreeNode<T> > *node;
	flens::DenseVector<flens::Array<T> > evaluate;
	HTuckerTreeNode<T> *save;
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.end(); TIT >= tree.begin(); TIT--){
		node = TIT.getNode();
		save = node->getContent();
		//std::cout << "htucker.tcc Evaluate at node " << node->getContent()->getIndex() << std::endl;
		//std::cout << "htucker.tcc UorB " << node->getContent()->getUorB() << std::endl;
		int evaluatepos = save->getIndex()[0]-1;
		if(node->isRoot()){
			evaluate.resize(1);
			flens::DenseVector<flens::Array<T> > intermediateresult(save->getUorB().numRows());
			intermediateresult = save->getUorB() * evalarrays[node->getfirstChild()->getContent()->getIndex()[0]-1];
			evaluate(1) = 0;
			flens::DenseVector<flens::Array<T> > &evaluate_right = evalarrays[node->getlastChild()->getContent()->getIndex()[0]-1];
			//std::cout << "intermediateresult = " << intermediateresult << std::endl;
			//std::cout << "intermediate: " << intermediateresult.firstIndex() << " " << intermediateresult.lastIndex() << std::endl;
			//std::cout << "evaluate_right" << evaluate_right << std::endl;
			//std::cout << "evaluate_right" << evaluate_right.firstIndex() << " " << evaluate_right.lastIndex() << std::endl;
			int interlen = intermediateresult.length();
			for(int j = 1; j<= interlen ; ++j){
				//std::cout << "j = " << j << std::endl;
				evaluate(1) += intermediateresult(j)*evaluate_right(j);
			}
			//std::cout << "root ende " << std::endl;
			delete [] evalarrays;
			return evaluate(1);		
		} else if(node->isLeaf()){
			assert(save->getIndex().length() == 1);
			int numcols = save->getUorB().numCols(); // in den Blättern geht das so, nicht aber in den inneren Knoten!
			evalarrays[evaluatepos] = flens::DenseVector<flens::Array<T> >(numcols,1);
			//std::cout << node->getContent()->getROA().A.minval << std::endl;
			//std::cout << (node->getContent()->getIndex()) << std::endl;
			//std::cout << (node->getContent()->getIndex())[0] - 1 << std::endl;
			int val = index[(save->getIndex())[0] - 1];
			// std::cout << "here" << std::endl;
			evalarrays[evaluatepos] = (save->getUorB())(val,_);
		} else{
			//std::cout << "inner " << std::endl;
			//std::cout << node->getContent()->getIndex() << std::endl;
			int numblocks = save->getNumRows();
            evaluate.resize(0);
			evaluate = flens::DenseVector<flens::Array<T> >(numblocks,1);
			int rcnumel = save->getRightChildNumRows();

			flens::DenseVector<flens::Array <T> > intermediateresult(rcnumel);
			flens::DenseVector<flens::Array <T> > &evalVectorlc = evalarrays[node->getlastChild()->getContent()->getIndex()[0]-1];
			flens::DenseVector<flens::Array <T> > &evalVectorfc = evalarrays[node->getfirstChild()->getContent()->getIndex()[0]-1];
			for(int i= 1; i<= numblocks; ++i){
				//std::cout << "     i = " << i << "  intermediateblock: " << rcnumel << " x " << node->getContent()->getUorB().numCols() << "  evaluate: " << (node->getfirstChild()->getContent()->getEvaluate()).length() << std::endl;
				intermediateresult = (save->getUorB())(_((i-1)*rcnumel + 1,i*rcnumel),_)  * evalVectorfc;
				evaluate(i) = 0;
				//std::cout << "here: " << intermediateresult.firstIndex() << "  " << intermediateresult.lastIndex() << "  " << node->getlastChild()->getContent()->getEvaluate().firstIndex() << "  " << node->getlastChild()->getContent()->getEvaluate().lastIndex() << std::endl;
				int interli = intermediateresult.lastIndex();
				for(int j = intermediateresult.firstIndex(); j<= interli ; ++j){
					evaluate(i) += intermediateresult(j)*evalVectorlc(j);
				}
			}
            evalarrays[evaluatepos].resize(evaluate.length());
			evalarrays[evaluatepos] = evaluate;
			//std::cout << "inner done" << std::endl;
		}
		//std::cout << "Evaluate: " << node->getContent()->getEvaluate() << std::endl;
	}
    return evaluate(1);
}

template <typename T>
flens::DenseVector<flens::Array<T> >
HTuckerTree<T>::vec_evaluate(const DimensionIndex &index, const int vardim) const{
    using flens::_;

	//std::cout << "vec_evaluate " << index << " vardim = " << vardim << std::endl;
	flens::DenseVector<flens::Array<T> > ret;
	flens::DenseVector<flens::Array<T> > * evalarrays = new flens::DenseVector<flens::Array<T> >[d];
	flens::GeMatrix<flens::FullStorage<T,cxxblas::RowMajor> > specialdim, tmpmat;
	GeneralTreeNode<HTuckerTreeNode<T> > *node;
	flens::DenseVector<flens::Array<T> > evaluate;


	DimensionIndex vdim(1);
	vdim[0] = vardim;


	HTuckerTreeNode<T> *save;
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.end(); TIT >= tree.begin(); TIT--){
		node = TIT.getNode();
		save = node->getContent();
		//std::cout << "htucker.tcc Evaluate at node " << node->getContent()->getIndex() << std::endl;
		//std::cout << "htucker.tcc UorB " << node->getContent()->getUorB() << std::endl;
		int evaluatepos = save->getIndex()[0]-1;
		if(node->isRoot()){
			//std::cout << "root" << std::endl;
			DimensionIndex &lcindex = node->getlastChild()->getContent()->getIndex();
			DimensionIndex &fcindex = node->getfirstChild()->getContent()->getIndex();
			flens::DenseVector<flens::Array<T> > intermediateresult(save->getUorB().numRows());

			if(lcindex.intersect(vdim).length() > 0){
				flens::DenseVector<flens::Array<T> > intermediateresult(save->getUorB().numRows());
				intermediateresult = save->getUorB() * evalarrays[node->getfirstChild()->getContent()->getIndex()[0]-1];
				ret = specialdim * intermediateresult;
			} else {
				flens::DenseVector<flens::Array<T> > intermediateresult(save->getUorB().numCols());
				intermediateresult = flens::transpose(save->getUorB())*evalarrays[node->getlastChild()->getContent()->getIndex()[0]-1];
				ret = specialdim*intermediateresult;
			}
			delete [] evalarrays;
			//std::cout << "bevore return: " << ret << "  fi = " << ret.firstIndex() << "  li = " << ret.lastIndex() << std::endl;
			return ret;		
		} else if(node->isLeaf()){
			assert(save->getIndex().length() == 1);
			if(save->getIndex()[0] == vardim){
				ret = flens::DenseVector<flens::Array<T> >(save->getUorB().numRows());
				specialdim = save->getUorB();
			} else {
				int numcols = save->getUorB().numCols(); // in den Blättern geht das so, nicht aber in den inneren Knoten!
			evalarrays[evaluatepos] = flens::DenseVector<flens::Array<T> >(numcols,1);
			//std::cout << node->getContent()->getROA().A.minval << std::endl;
			//std::cout << (node->getContent()->getIndex()) << std::endl;
			//std::cout << (node->getContent()->getIndex())[0] - 1 << std::endl;
			int val = index[(save->getIndex())[0] - 1];
			// std::cout << "here" << std::endl;
			evalarrays[evaluatepos] = (save->getUorB())(val,_);
			}
			
		} else{
			//std::cout << "inner " << std::endl;
			//std::cout << node->getContent()->getIndex() << std::endl;
			DimensionIndex &lcindex = node->getlastChild()->getContent()->getIndex();
			DimensionIndex &fcindex = node->getfirstChild()->getContent()->getIndex();

			if(	lcindex.intersect(vdim).length() > 0){
				tmpmat = specialdim;
				int numblocks = save->getNumRows();
				specialdim = flens::GeMatrix<flens::FullStorage<T,cxxblas::RowMajor> >(tmpmat.numRows(),numblocks);
				int lcnumel = save->getLeftChildNumRows();
				int rcnumel = save->getRightChildNumRows();
				flens::DenseVector<flens::Array <T> > intermediateresult(rcnumel);
				flens::DenseVector<flens::Array <T> > &evalVectorfc = evalarrays[node->getfirstChild()->getContent()->getIndex()[0]-1];
				for(int i= 1; i<= numblocks; ++i){
					//std::cout << "     i = " << i << "  intermediateblock: " << rcnumel << " x " << node->getContent()->getUorB().numCols() << "  evaluate: " << (node->getfirstChild()->getContent()->getEvaluate()).length() << std::endl;
					intermediateresult = (save->getUorB())(_((i-1)*rcnumel + 1,i*rcnumel),_)  * evalVectorfc;
					//std::cout << "here: " << intermediateresult.firstIndex() << "  " << intermediateresult.lastIndex() << "  " << node->getlastChild()->getContent()->getEvaluate().firstIndex() << "  " << node->getlastChild()->getContent()->getEvaluate().lastIndex() << std::endl;
					specialdim(_,i) = tmpmat*intermediateresult;
				}
			} else if(fcindex.intersect(vdim).length() > 0){
				tmpmat = specialdim;
				int numblocks = save->getNumRows();
				specialdim = flens::GeMatrix<flens::FullStorage<T,cxxblas::RowMajor> >(tmpmat.numRows(),numblocks);
				int lcnumel = save->getLeftChildNumRows();
				int rcnumel = save->getRightChildNumRows();
				flens::DenseVector<flens::Array <T> > intermediateresult(rcnumel);
				flens::DenseVector<flens::Array <T> > &evalVectorlc = evalarrays[node->getlastChild()->getContent()->getIndex()[0]-1];
				for(int i= 1; i<= numblocks; ++i){
					//std::cout << "     i = " << i << "  intermediateblock: " << rcnumel << " x " << node->getContent()->getUorB().numCols() << "  evaluate: " << (node->getfirstChild()->getContent()->getEvaluate()).length() << std::endl;
					intermediateresult = flens::transpose((save->getUorB())(_((i-1)*rcnumel + 1,i*rcnumel),_))  * evalVectorlc;
					//std::cout << "here: " << intermediateresult.firstIndex() << "  " << intermediateresult.lastIndex() << "  " << node->getlastChild()->getContent()->getEvaluate().firstIndex() << "  " << node->getlastChild()->getContent()->getEvaluate().lastIndex() << std::endl;
					flens::DenseVector<flens::Array<T> > test;
					test = (tmpmat*intermediateresult);
					//std::cout << "evalVector length = " << evalVectorlc.length() << std::endl;
					//std::cout << "test.length = " << test.length() << std::endl;
					//std::cout << "tempmat: " << tmpmat.numRows() << " x " << tmpmat.numCols() << std::endl;
					//std::cout << "specialdim = " << specialdim.numRows() << " x " << specialdim.numCols() << std::endl;
					specialdim(_,i) = tmpmat*intermediateresult;
				}
			} else {
				int numblocks = save->getNumRows();
				evaluate = flens::DenseVector<flens::Array<T> >(numblocks,1);
				int lcnumel = save->getLeftChildNumRows();
				int rcnumel = save->getRightChildNumRows();

				flens::DenseVector<flens::Array <T> > intermediateresult(rcnumel);
				flens::DenseVector<flens::Array <T> > &evalVectorlc = evalarrays[node->getlastChild()->getContent()->getIndex()[0]-1];
				flens::DenseVector<flens::Array <T> > &evalVectorfc = evalarrays[node->getfirstChild()->getContent()->getIndex()[0]-1];
				for(int i= 1; i<= numblocks; ++i){
					//std::cout << "     i = " << i << "  intermediateblock: " << rcnumel << " x " << node->getContent()->getUorB().numCols() << "  evaluate: " << (node->getfirstChild()->getContent()->getEvaluate()).length() << std::endl;
					intermediateresult = (save->getUorB())(_((i-1)*rcnumel + 1,i*rcnumel),_)  * evalVectorfc;
					evaluate(i) = 0;
					//std::cout << "here: " << intermediateresult.firstIndex() << "  " << intermediateresult.lastIndex() << "  " << node->getlastChild()->getContent()->getEvaluate().firstIndex() << "  " << node->getlastChild()->getContent()->getEvaluate().lastIndex() << std::endl;
					int interli = intermediateresult.lastIndex();
					for(int j = intermediateresult.firstIndex(); j<= interli ; ++j){
						evaluate(i) += intermediateresult(j)*evalVectorlc(j);
					}
				}
				evalarrays[evaluatepos] = evaluate;
			
			}
		}
		//std::cout << "Evaluate: " << node->getContent()->getEvaluate() << std::endl;
	}
}



template <typename T>
void
HTuckerTree<T>::orthogonalize(){
    using flens::_;

	GeneralTreeNode<HTuckerTreeNode<T> > *node;
	HTuckerTreeNode<T> *save,*last,*first;
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.end(); TIT >= tree.begin(); TIT--){
		node = TIT.getNode();
		save = node->getContent();
		last = node->getlastChild()->getContent();
		first = node->getfirstChild()->getContent();
		//std::cout << "node = " << node->getContent()->getIndex() << std::endl;
		if(node->isLeaf()){
			//nothing happens
		} else {

			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > QR,temp1;
			flens::DenseVector<flens::Array<T> > tau;

			//std::cout << "orthogonalize 1. get UorB from right Child and transform into non-block-format" << std::endl;
			if(node->getlastChild()->isLeaf()){
				QR = last->getUorB();
			} else {
				int  numel = last->getNumRows();
				int lcnumel = last->getLeftChildNumRows();
				int rcnumel = last->getRightChildNumRows();
				QR = flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> >(lcnumel*rcnumel, numel);
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &lcUorB = node->getlastChild()->getContent()->getUorB();
				for(int i = 1; i <= lcnumel; ++i){
					for(int j = 1; j <= numel; ++j){
						QR(_((i-1)*rcnumel + 1, i*rcnumel),j) = lcUorB(_((j-1)*rcnumel + 1, j*rcnumel),i);
					}
				}
			}
			//std::cout << "orthogonalize 2. Comute QR decomposition of Right child" << std::endl;
			flens::qrf(QR,tau);
			//std::cout << "QR: " << QR.numRows() << " x " << QR.numCols() << std::endl;
			
			int basisnumel = save->getNumRows();
			int basisrcnumel = save->getRightChildNumRows();
			int basislcnumel = save->getLeftChildNumRows();

			int m = QR.numRows();
			int n = QR.numCols();

			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > intermediate(std::min(m,n)*basisnumel,basislcnumel);
			//std::cout << "orthogonalize 3. Compute R*B blockwise for all Blocks in B" << std::endl;
			for(int i = 1; i <= basisnumel; ++i){
				temp1 = (save->getUorB())(_((i-1)*basisrcnumel + 1, i*basisrcnumel),_);
				if(m == n){
					flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,QR.upper(),temp1);
					intermediate(_((i-1)*m + 1, i*m),_) = temp1;
				} else if(m > n){
					flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,QR(_(1,n),_).upper(),temp1);
					intermediate(_((i-1)*n+1,i*n),_) = temp1;
				} else {
					flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp(n,n);
					tmp(_(1,m),_) = QR;
					flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,tmp.upper(),temp1);
					intermediate(_((i-1)*m+1,i*m),_) = temp1(_(1,m),_);
				}
			}


			//R ist versorgt, nun berechnen wir Q und speichern das im Kindknoten!
			//std::cout << "orthogonalize 4. Compute Q of rigth child" << std::endl;
			if( n > m){
                temp1.resize(m, m);
				temp1 = QR(_(1,m),_(1,m));
                QR.resize(m, m);
				QR = temp1;
				flens::orgqr(QR,tau);
			} else {
				flens::orgqr(QR,tau);
			}
			// Jetzt muss QR in Blockformat gebracht werden
			//std::cout << "orthogonalize 5. Transform Q of right child into Block format" << std::endl;
			if(node->getlastChild()->isLeaf()){
				last->setUorB(QR);
				last->setNumRows(QR.numCols());
			} else {
				//std::cout << "     QR : " << QR.numRows() << " x " << QR.numCols() << std::endl;
				//std::cout << "     right_rcnumel = " << node->getlastChild()->getContent()->UorB_rcnumel << "   right_lcnumel = " <<  node->getlastChild()->getContent()->UorB_lcnumel << "  right_numel = " << node->getlastChild()->getContent()->UorB_numel <<  std::endl;
				int lcnumel = last->getLeftChildNumRows();
				int rcnumel = last->getRightChildNumRows();
                temp1.resize(std::min(n,m)*rcnumel, lcnumel);
				for(int i = 1; i <= std::min(n,m) ; ++i){
					for(int j = 1; j <= lcnumel; ++j){
						temp1(_((i-1)*rcnumel + 1, i*rcnumel),j) = QR(_((j-1)*rcnumel + 1,j*rcnumel),i);
					}
				}
				last->setUorB(temp1);
				last->setNumRows(std::min(n,m));
			}
			//hier sind wir mit dem rechten Kind fertig... weiter geht es mit dem linken...
			//std::cout << "orthogonalize 6. get UorB from left child and transform in non-block-format" << std::endl;
			if(node->getfirstChild()->isLeaf()){
                QR.resize(first->getUorB().numRows(),
                          first->getUorB().numCols());
				QR = first->getUorB();
			} else {
				int  numel = first->getNumRows();
				int lcnumel = first->getLeftChildNumRows();
				int rcnumel = first->getRightChildNumRows();
                QR.resize(lcnumel*rcnumel, numel);
				for(int i = 1; i <= lcnumel; ++i){
					for(int j = 1; j <= numel; ++j){
						QR(_((i-1)*rcnumel + 1, i*rcnumel),j) = (first->getUorB())(_((j-1)*rcnumel + 1, j*rcnumel),i);
					}
				}
			}
			flens::DenseVector<flens::Array<T> > tau1;
			//std::cout << "orthogonalize 7. Compute QR-Decomposition of left child " << std::endl;
			flens::qrf(QR,tau1);
			int ml = QR.numRows();
			int nl = QR.numCols();
			
			//std::cout << "orthogonalize 8. Compute B*R^T for all Blocks of B " << std::endl;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > final(std::min(m,n)*basisnumel, std::min(ml,nl));
			for(int i = 1; i <= basisnumel; ++i){
                temp1.resize(std::min(m,n), intermediate.numCols());
				temp1 = intermediate(_((i-1)*std::min(m,n) + 1, i*std::min(m,n)),_);
				if(nl == ml){
					flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,QR.upper(),temp1);
					final(_((i-1)*std::min(m,n) + 1,i*std::min(m,n)),_) = temp1;
				} else if (ml > nl){
					flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,QR(_(1,nl),_).upper(),temp1);
					final(_((i-1)*std::min(m,n)+1,i*std::min(m,n)),_) = temp1;
				} else {
					flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp(nl,nl);
					tmp(_(1,ml),_) = QR;
					flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,tmp.upper(),temp1);
					final(_((i-1) * std::min(m,n) + 1, i*std::min(m,n)),_) = temp1(_,_(1,ml));
				}
			}

			save->setUorB(final);
			save->setRightChildNumRows(std::min(m,n));
			save->setLeftChildNumRows(std::min(ml,nl));
			
			//std::cout << "orthogonalize 9. Compute Q of left child" << std::endl;
			//R ist versorgt, nun berechnen wir Q und speichern das im Kindknoten!
			if( nl > ml){
                temp1.resize(ml, ml);
				temp1 = QR(_(1,ml),_(1,ml));
                QR.resize(ml, ml);
				QR = temp1;
				orgqr(QR,tau1);
			} else {
				flens::orgqr(QR,tau1);
			}

			// Jetzt muss QR in Blockformat gebracht werden

			//std::cout << "orthogonalize 10. Transporm Q into Block format" << std::endl;
			if(node->getfirstChild()->isLeaf()){
				first->setUorB(QR);
				first->setNumRows(QR.numCols());
			} else {

				int lcnumel = first->getLeftChildNumRows();
				int rcnumel = first->getRightChildNumRows();
                temp1.resize(std::min(nl,ml)*rcnumel, lcnumel);
				for(int i = 1; i <= std::min(nl,ml) ; ++i){
					for(int j = 1; j <= lcnumel; ++j){
						temp1(_((i-1)*rcnumel + 1, i*rcnumel),j) = QR(_((j-1)*rcnumel + 1,j*rcnumel),i);
					}
				}
				first->setUorB(temp1);
				first->setNumRows(std::min(nl,ml));
			}
		}
	}
}


template <typename T>
void
HTuckerTree<T>::orthogonalize_svd(std::vector
                                  <flens::DenseVector
                                  <flens::Array<T> > >& sigmas,
                                  bool isorth)
{
    if (!sigmas.size()) sigmas.resize(dim());
    assert(sigmas.size()==(unsigned) dim());

    using flens::_;

    typedef typename flens::GeMatrix<
                     flens::FullStorage<T,cxxblas::ColMajor> >  Matrix;
    typedef typename flens::DenseVector<flens::Array<T> >       Vector;

    HTuckerTree<T> Stree;
    Stree.set_tree(*this);
    if (!isorth) {
       this->orthogonalize();
    }

    HTuckerTree<T> gram = gramians_orthogonal(*this);
    GeneralTreeNode<HTuckerTreeNode<T> > *node;

    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TITgram = gram.getGeneralTree().begin();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITs= Stree.getGeneralTree().begin();
    for(; TITgram <= gram.getGeneralTree().end(); TITgram ++,TITs++){
        Matrix& M = TITgram.getNode()->getContent()->getUorB();
        auto n    = M.numRows();
        auto m    = M.numCols();
        Matrix U(m, std::min(m, n)), Vt;
        Vector work, s(std::min(m, n));
        flens::lapack::svd(flens::lapack::SVD::Job::Save,
                           flens::lapack::SVD::Job::None,
                           M,
                           s,U,Vt,
                           work);

        if (TITgram.getNode()->isLeaf()) {
            sigmas[TITgram.getNode()
                   ->getContent()->getIndex()[0]-1].resize(s.length());
            sigmas[TITgram.getNode()->getContent()->getIndex()[0]-1] = s;
        }
		TITs.getNode()->getContent()->setUorB(U);
	}
    }

    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = this->getGeneralTree().begin();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITs = Stree.getGeneralTree().begin();
	for(; TIT <= this->getGeneralTree().end(); TIT++,TITs++){
		node = TIT.getNode();
		
		if(node->isLeaf()){
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp;
			flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,node->getContent()->getUorB(),TITs.getNode()->getContent()->getUorB(),0.0,tmp);
			node->getContent()->setUorB(tmp);
			node->getContent()->setNumRows(tmp.numCols());
		} else {
			int rcnumel = node->getContent()->getRightChildNumRows();
			int lcnumel = node->getContent()->getLeftChildNumRows();
			
			int numel = node->getContent()->getNumRows();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Sref = TITs.getNode()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Slref = TITs.getNode()->getfirstChild()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Srref = TITs.getNode()->getlastChild()->getContent()->getUorB();
			int newnumel = Sref.numCols();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Bref = node->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp(rcnumel*newnumel,lcnumel);
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Bt(rcnumel*lcnumel, numel);
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > BtSt(rcnumel*lcnumel, newnumel);

            /* Bt*St */
            for(int i = 1; i <= numel; ++i) {
                for(int j = 1; j <= lcnumel; ++j) {
                    Bt(_((j-1)*rcnumel+1, j*rcnumel), i) = Bref(_((i-1)*rcnumel+1, i*rcnumel), j);
                }
            }
            flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans, 1.0, Bt, Sref,0.0, BtSt);
            for(int i = 1; i <= newnumel; ++i) {
                for(int j = 1; j <= lcnumel; ++j) {
                    tmp(_((i-1)*rcnumel+1, i*rcnumel), j) = BtSt(_((j-1)*rcnumel+1, j*rcnumel), i);
                }
            }

			//hier haben wir nun als B_t*S_t berechnet. Fehlt noch (S_tl^T x S_tr^T)*B_t*S_t
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > newB(newnumel*Srref.numCols(),Slref.numCols());
			for(int i = 1; i <= newnumel; ++i){
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp1;
				flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,tmp(_((i-1)*rcnumel+1,i*rcnumel),_),Slref,0.0,tmp1);
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp2;
				flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,Srref,tmp1,0.0,tmp2);
				newB(_((i-1)*Srref.numCols()+1,i*Srref.numCols()),_) = tmp2;
			}
			node->getContent()->setUorB(newB);
			node->getContent()->setNumRows(newnumel);
			node->getContent()->setLeftChildNumRows(Slref.numCols());
			node->getContent()->setRightChildNumRows(Srref.numCols());
		
		}
		
		
	}
    }
}




template <typename T>
T
HTuckerTree<T>::L2norm() const{
	T erg = ScalarProduct(*this);
	return sqrt(fabs(erg));
}


template <typename T>
T
HTuckerTree<T>::L2normorthogonal() const
{
    return flens::blas::nrm2(tree.root->getContent()->getUorB().vectorView());
}


template <typename T>
T
HTuckerTree<T>::ScalarProduct(const HTuckerTree<T> & anothertree) const
{
    using flens::_;

    typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > Matrix;

    assert(this->getGeneralTree().root->getContent()->getIndex().length() ==
    anothertree.getGeneralTree().root->getContent()->getIndex().length());

    GeneralTree<HTuckerTreeNode<T> > tmp = anothertree.getGeneralTree();
    GeneralTreeNode<HTuckerTreeNode<T> > *nodethis, *nodeanother, *nodetmp;
    for (GeneralTreeIterator<HTuckerTreeNode<T> >
        TITthis = tree.end(),
        TITanother = anothertree.getGeneralTree().end(),
        TITtmp = tmp.end(); TITthis >= tree.begin();
        TITthis--,TITanother--,TITtmp--) {

        nodethis = TITthis.getNode();
        nodeanother = TITanother.getNode();
        nodetmp = TITtmp.getNode();

        if (nodethis->isLeaf()) {
            assert(nodethis->getContent()->getUorB().numRows() ==
                   nodeanother->getContent()->getUorB().numRows());
            Matrix tmpM;
            Matrix& Bleft  = nodethis->getContent()->getUorB();
            Matrix& Bright = nodeanother->getContent()->getUorB();
            flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                            1., Bleft, Bright,
                            0., tmpM);
            nodetmp->getContent()->setUorB(tmpM);

        } else {
            int rightnumel   = nodeanother->getContent()->getNumRows();
            int rightrcnumel = nodeanother->getContent()
                               ->getRightChildNumRows();

            int leftnumel    = nodethis->getContent()->getNumRows();
            int leftrcnumel  = nodethis->getContent()
                               ->getRightChildNumRows();
            int leftlcnumel  = nodethis->getContent()
                               ->getLeftChildNumRows();

            Matrix tensorblock(leftrcnumel*rightnumel, leftlcnumel);
            Matrix tensor(leftrcnumel*leftlcnumel, rightnumel);
            Matrix tmp(rightrcnumel, leftlcnumel);

            Matrix& Bleft  = nodethis->getContent()->getUorB();
            Matrix& Bright = nodeanother->getContent()->getUorB();

            Matrix Bleftnonblock(leftrcnumel*leftlcnumel, leftnumel);

            Matrix& Ml = nodetmp->getfirstChild()->getContent()->getUorB();
            Matrix& Mr = nodetmp->getlastChild()->getContent()->getUorB();

            Matrix result(leftnumel, rightnumel);

            for (int i=1; i<=rightnumel; ++i) {
                flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                1., Bright(_((i-1)*rightrcnumel+1,
                                              i*rightrcnumel), _),
                                Ml, 0., tmp);
                flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                1., Mr,
                                tmp, 0., tensorblock(_((i-1)*leftrcnumel+1,
                                                        i*leftrcnumel), _));
            }

            for (int i=1; i<=leftnumel; ++i) {
                for (int j=1; j<=leftrcnumel; ++j) {
                    Bleftnonblock(_((j-1)*leftlcnumel+1,
                                     j*leftlcnumel), i) =
                    Bleft((i-1)*leftrcnumel+j, _);
                }
            }

            for (int i=1; i<=rightnumel; ++i) {
                for (int j=1; j<=leftrcnumel; ++j) {
                    tensor(_((j-1)*leftlcnumel+1, j*leftlcnumel), i) =
                    tensorblock((i-1)*leftrcnumel+j, _);
                }
            }

            flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                            1., Bleftnonblock,
                            tensor, 0., result);

            nodetmp->getContent()->setUorB(result);
        }
    }

    return (tmp.root->getContent()->getUorB())(1,1);
}


template <typename T>
void 
HTuckerTree<T>::generateTofElementary(DenseVectorList<T> & list,const int k,const int d){ //flens::DenseVector list is non-const, because () operator in DVL is non-const
    typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> >      GEMatrix;
    flens::Underscore<typename GEMatrix::IndexType>        _;

	GeneralTreeNode<HTuckerTreeNode<T> > *node;
	HTuckerTreeNode<T> *content;
	this->d = d;
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.end(); TIT >= tree.begin(); TIT--){
		node = TIT.getNode();
		content = node->getContent();
		if(node->isLeaf()){
			int pos = (content->getIndex())[0];
			GEMatrix mat3((*list((pos-1)*k+1)).length(),k);
			for(int i = 1; i<= k; ++i){
				mat3(_,i) = *list((pos-1)*k+i);
			}
			content->setUorB(mat3);
			content->setNumRows(k);
		} else if(node->isRoot()) {
			GEMatrix mat(k,k);
			for(int i = 1; i<=k; ++i){
				mat(i,i) = 1.0;
			}
			content->setUorB(mat);
			content->setNumRows(1);
			content->setRightChildNumRows(k);
			content->setLeftChildNumRows(k);
		} else {
			GEMatrix mat2(k*k,k);
			for(int i = 1; i<= k; ++i){
				mat2((i-1)*k+i,i) = 1.0;
			}
			content->setUorB(mat2);
			content->setNumRows(k);
			content->setRightChildNumRows(k);
			content->setLeftChildNumRows(k);
		}
	}

}



template <typename T>
T 
HTuckerTree<T>::average_rank() const{
	int count = 0;
	T rank = 0.0;
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.end(); TIT >= tree.begin(); TIT--){
		if(TIT.getNode()->isLeaf()){
			rank += TIT.getNode()->getContent()->getUorB().numCols();
		} else  {
			rank += TIT.getNode()->getContent()->getNumRows();
		}
		count++;
	}

	return rank/(count +0.0);
}

template <typename T>
int
HTuckerTree<T>::max_rank() const{
	int rank = 0;
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.end(); TIT >= tree.begin(); TIT--){
		int newrank;
		if(TIT.getNode()->isLeaf()){
			 newrank = TIT.getNode()->getContent()->getUorB().numCols();
		} else  {
			 newrank = TIT.getNode()->getContent()->getNumRows();
		}
		if(newrank > rank){
			rank = newrank;
		}
	}
	return rank;
}


template <typename T>
int
HTuckerTree<T>::max_n() const{
    int rank = 0;
    for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.end(); TIT >= tree.begin(); TIT--){
        int newrank;
        if(TIT.getNode()->isLeaf()){
            newrank = TIT.getNode()->getContent()->getUorB().numRows();
        } else {
            break;
        }
        if(newrank > rank){
            rank = newrank;
        }
    }

    return rank;
}


template <typename T>
T 
HTuckerTree<T>::effective_rank() const{
	int storage = 0;
	int n = 0;
	int count = 0;
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.end(); TIT >= tree.begin(); TIT--){
		if(TIT.getNode()->isLeaf()){
			storage = storage + TIT.getNode()->getContent()->getUorB().numRows()*TIT.getNode()->getContent()->getUorB().numCols();
			n = n + TIT.getNode()->getContent()->getUorB().numRows();
			count ++;
		} else {
			storage = storage + TIT.getNode()->getContent()->getNumRows()*TIT.getNode()->getContent()->getLeftChildNumRows()*TIT.getNode()->getContent()->getRightChildNumRows();
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
}



template <typename T>
T
HTuckerTree<T>::Linfnorm(HTuckerTree<T> & anothertree,const DimensionIndex & minval,const DimensionIndex & maxval,const int n){
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
}

template <typename T>
T
HTuckerTree<T>::Linfnorm(const DimensionIndex & minval,const DimensionIndex & maxval,const int n){
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

}

template <typename T>
template <typename TensorFunction>
T
HTuckerTree<T>::Linfnorm(const tensor<TensorFunction> & tens,const DimensionIndex & minval, const DimensionIndex & maxval,const int n){
	T maxi = -1;
	int r;
	T res;
	DimensionIndex eval(minval.length());
	for(int i = 1; i<= n; ++i){
		eval.setRandom(minval,maxval);
		res = abs(evaluate(eval) - tens(eval));
		if(res > maxi){
			maxi = res;
		}
	}
	return maxi;
}

template <typename T>
DimensionIndex
HTuckerTree<T>::getmax() const{
	DimensionIndex DI(this->d);
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.end(); TIT >= tree.begin(); TIT--){
		if(TIT.getNode()->isLeaf()){
			DI[(TIT.getNode()->getContent()->getIndex())[0]-1] = TIT.getNode()->getContent()->getUorB().numRows();
		}
	}
	return DI;

}

template <typename T>
int
HTuckerTree<T>::dim() const{
	return d;
}

template <typename T>
void
HTuckerTree<T>::setdim(const int dim){
	d = dim;
}


template <typename T>
void
HTuckerTree<T>::scal(T alpha)
{
    typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Matrix;

    Matrix& B = const_cast<Matrix&>
                (getGeneralTree().root->getContent()->getUorB());
    B *= alpha;
}


template <typename T>
int
HTuckerTree<T>::depth() const
{
    return this->getGeneralTree().end().getNode()->level();
}


template <typename T>
HTuckerTree<T>&
HTuckerTree<T>::operator=(const HTuckerTree<T>& copy)
{
    d = copy.dim();
    tree = copy.getGeneralTree();
    return *this;
}


template <typename T>
template <typename TensorFunction>
void 
HTuckerTree<T>::approximate(const tensor<TensorFunction> & tf, const int rank, const int l){
	assert(tf.dim() == d);
	GeneralTree<ROA<T,TensorFunction> > roatree(tree);

	GeneralTreeNode<ROA<T,TensorFunction> > *node;
    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT2 = tree.begin();
	GeneralTreeIterator<ROA<T,TensorFunction> > TIT= roatree.begin();
    for(;TIT <= roatree.end(); TIT++,TIT2++){
					
		node = TIT.getNode();
		if(!(node->isRoot())){
			ROA<T,TensorFunction> roa(tf,TIT2.getNode()->getContent()->getIndex(),TIT2.getNode()->getnextSibling()->getContent()->getIndex(),d);
			node->setContent(roa);
			std::cout << std::endl <<  "HTuckerTree<T>::approximate: fatherpivs = " << node->getParent()->getContent()->pivots << std::endl <<std::endl;
			node->getContent()->approximate(rank,node->getParent()->getContent()->pivots, l);
					
			std::cout << "after approximation" << node->getContent()->pivots << std::endl;
		} else {
			ROA<T,TensorFunction> roa(tf,TIT2.getNode()->getContent()->getIndex(),DimensionIndex(),d);
			node->setContent(roa);
		}
	}
    }
	this->setUorB(roatree);


	//node->getContent()->approximate(rank,TIT2.getN->getParent()->getContent()->getROA().pivots);
}


template <typename T>
template <typename TensorFunction>
void 
HTuckerTree<T>::approximate(const tensor<TensorFunction> &tf, const double epsilon, const int l){
	
	

	assert(tf.dim() == d);
	GeneralTree<ROA<T,TensorFunction> > roatree(tree);

	GeneralTreeNode<ROA<T,TensorFunction> > *node;
    {
	GeneralTreeIterator<ROA<T,TensorFunction> > TIT;
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT2;
    for(TIT = roatree.begin(), TIT2 = tree.begin();
		TIT <= roatree.end(); TIT++,TIT2++){
					
		node = TIT.getNode();
		if(!(node->isRoot())){
			ROA<T,TensorFunction> roa(tf,TIT2.getNode()->getContent()->getIndex(),TIT2.getNode()->getnextSibling()->getContent()->getIndex(),d);
			node->setContent(roa);			
			node->getContent()->approximate(epsilon, node->getParent()->getContent()->pivots, l);
		} else {
			ROA<T,TensorFunction> roa(tf,TIT2.getNode()->getContent()->getIndex(),DimensionIndex(),d);
			node->setContent(roa);		
		}
	}
    }
	this->setUorB(roatree);
	//for(GeneralTreeIterator<HTuckerTreeNode<T,EVALTYPE> > TIT = tree.begin(); TIT <= tree.end(); TIT++){
	//	node = TIT.getNode();
	//	if(!(node->isRoot())){
	//		//std::cout << "HTUCKERTREE.tcc: approximate " << node->getContent()->getIndex() << std::endl;
	//		//std::cout << "HTUCKERTREE.tcc: roa " << std::endl;
	//		node->getContent()->getROA().print();
	//		node->getContent()->getROA().approximate(epsilon,node->getParent()->getContent()->getROA().pivots);
	//	}
	//}
}


template <typename T>
template <typename TensorFunction>
void
HTuckerTree<T>::setUorB(const GeneralTree<ROA<T,TensorFunction> > & gt){
    using flens::_;

	typedef double _T;
	double eps = gt.root->getContent()->SVD_epsilon;
	int d = this->dim();
	GeneralTreeNode<HTuckerTreeNode<double> > *node;
	GeneralTreeNode<ROA<_T,TensorFunction> > *roaNode;
    {
	GeneralTreeIterator<HTuckerTreeNode<double> > TIT;
    GeneralTreeIterator<ROA<_T,TensorFunction> > TIT2;
    for(TIT = this->getGeneralTree().end(),
		TIT2 = gt.end(); TIT >= this->getGeneralTree().begin(); TIT--,TIT2--){
		
		node = TIT.getNode();
		roaNode = TIT2.getNode();
		if(node->isLeaf()){
			int minval = roaNode->getContent()->A.getmin()[(node->getContent()->getIndex())[0]-1];
			int maxval = roaNode->getContent()->A.getmax()[(node->getContent()->getIndex())[0]-1];
			int pivlen = roaNode->getContent()->pivots.length();
			ROA<_T,TensorFunction> *roa = roaNode->getContent();
			const DimensionIndex &index = node->getContent()->getIndex();
			flens::GeMatrix<flens::FullStorage<_T,cxxblas::ColMajor> > U(maxval - minval + 1,pivlen);
			typedef typename flens::GeMatrix<flens::FullStorage<_T,cxxblas::ColMajor> >::IndexType Itype;
			
			for(int i = 0; i < pivlen; ++i){
				U(_,(Itype)(i+1)) = roa->evaluate(roa->pivots[i],index,(index)[0],minval,maxval);
			}
			flens::GeMatrix<flens::FullStorage<_T, cxxblas::ColMajor> > Utilde;
			Utilde = U*flens::transpose(roa->V);
			node->getContent()->setUorB(Utilde);
			node->getContent()->setNumRows(Utilde.numCols());
		} else {
			int lcnumel = roaNode->getfirstChild()->getContent()->pivots.length();
			int rcnumel = roaNode->getlastChild()->getContent()->pivots.length();
			int numel;
			
			node->getContent()->setLeftChildNumRows(lcnumel);
			node->getContent()->setRightChildNumRows(rcnumel);
			if(node->isRoot()){
				numel = 1;
			} else {
				numel = roaNode->getContent()->pivots.length();
			}
			node->getContent()->setNumRows(numel);
			
			flens::GeMatrix<flens::FullStorage<_T,cxxblas::ColMajor> > Tinner(numel * rcnumel,lcnumel);
			DimensionIndex idx(this->dim());
			flens::GeMatrix<flens::FullStorage<_T,cxxblas::ColMajor> > Amat(rcnumel,lcnumel);
			ROA<T,TensorFunction> *roa = roaNode->getContent();
			ROA<T,TensorFunction> *fcroa = roaNode->getfirstChild()->getContent();
			ROA<T,TensorFunction> *lcroa = roaNode->getlastChild()->getContent();
			const DimensionIndex &fcindex = node->getfirstChild()->getContent()->getIndex();
			const DimensionIndex &lcindex = node->getlastChild()->getContent()->getIndex();
			flens::DenseVector<flens::Array<_T > > &sl = fcroa->s;
			flens::DenseVector<flens::Array<_T > > &sr = lcroa->s;

			for(int i = 0; i < numel; ++i){
				
				if(!(node->isRoot())){
					idx = (roa->pivots)[i];
				} 
				for(int j = 0; j < lcnumel; ++j){
					for(int k = 0; k < rcnumel; ++k){
						idx.setValue(fcroa->pivots[j],fcindex);
						idx.setValue(lcroa->pivots[k],lcindex);
						Amat(k+1,j+1) = (roa->A)(idx);
					}
				}
				flens::GeMatrix<flens::FullStorage<_T,cxxblas::ColMajor> > AU, AVT;
				flens::DenseVector<flens::Array<_T> > As;
				flens::svd(Amat,As,AU,AVT);
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > O1,O2;
				O1 = flens::transpose(lcroa->U)*AU;
				O2 = AVT*(fcroa->U);
				flens::GeMatrix<flens::FullStorage<_T,cxxblas::ColMajor> > o1help(O1.numRows(),O2.numCols());
				int minrows = min(O1.numCols(),O2.numCols());
				for(int j = 1; j <= minrows; ++j){ //should be 02.numCols(), but, As has only as many nonzero-columns as O1!
					if(j <= As.length()){
						o1help(_,j) = O1(_,j)*As(j);
					} else {
						o1help(_,j) = O1(_,j)*0;
					}

				}
				O1 = o1help * O2;
				O2 = O1;

				
			
				for(int j = 1; j<= O2.numRows(); ++j){
					for(int k = 1; k<= O2.numCols(); ++k){
						if(abs(sl(k)) > eps && abs(sr(j)) > eps){
							O2(j,k) = O1(j,k) / sl(k)/sr(j);
						} else {
							O2(j,k) = 0;
						}
					}
				}
				Tinner(_(i*rcnumel+1, (i+1)*rcnumel),_) = O2;
				
			} //end for(i = 0; i < numel; ++i){
			
			if(!(node->isRoot())){
				flens::GeMatrix<flens::FullStorage<_T,cxxblas::ColMajor> > V;
				//V = flens::transpose(node->getContent()->getROA().V); hier wird die falsche copy-funktion aufgerufen => copy direkt
				flens::blas::copy(cxxblas::Trans,roaNode->getContent()->V,V);
				flens::GeMatrix<flens::FullStorage<_T,cxxblas::ColMajor> > Ttilde(V.numCols() * rcnumel,lcnumel);
				assert(V.numRows() == numel);
				for(int i1 = 1; i1 <= rcnumel; ++i1){
					for(int i2 = 1; i2<= V.numCols(); ++i2){
						for(int j = 1; j<= lcnumel; ++j){
							Ttilde((i2-1)*rcnumel + i1,j) = 0;
							for(int k = 1; k <= numel; ++k){
								Ttilde((i2-1)*rcnumel + i1,j) += Tinner( (k-1) * rcnumel + i1,j) * V(k, i2);
							} 
						}
					}
				}
				node->getContent()->setUorB(Ttilde);
				node->getContent()->setNumRows(V.numCols());
			} else {
				node->getContent()->setUorB(Tinner);
			}

			
		
		} //end else {
	} // end for
    }
}

template <typename T>
void
HTuckerTree<T>::spy(const char* filename, const char* terminal, const char* options, const double width, const double height, const double pointsize) const{
	DimensionIndex base(d);
	base.setValueAsc();
	spy(base,filename,terminal,options,width,height,pointsize);
}

template <typename T>
void 
HTuckerTree<T>::spy(const DimensionIndex &subknoten,const char* filename,const char* terminal, const char* options, const double width, const double height, const double pointsize) const{
	int low;
	T maximum,minimum;
	GeneralTreeNode<HTuckerTreeNode<T> > *node;
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.begin(); TIT <= tree.end(); TIT++){
		if(TIT.getNode()->getContent()->getIndex() == subknoten){
			std::cout << subknoten << std::endl;
			low = TIT.getNode()->level();
			std::cout << low << std::endl;
			node = TIT.getNode();
			maximum = TIT.getNode()->getContent()->getUorB()(1,1);
			minimum = maximum;
			break;
		}
	
	}
	int high = tree.end().getNode()->level();
	std::cout << "high = " << high << std::endl;
	
	T val;
	HTuckerTree<T> newtree = subtree(*this,subknoten[0],subknoten[subknoten.length()-1]);
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &UorB = newtree.getGeneralTree().root->getContent()->getUorB();
	
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = newtree.getGeneralTree().begin(); TIT <= newtree.getGeneralTree().end(); TIT++){
		UorB = TIT.getNode()->getContent()->getUorB();
		for(int i = 1; i <= UorB.numRows(); ++i){
			for(int j = 1; j<= UorB.numCols(); ++j){
				val = UorB(i,j);
				if(val > maximum){
					maximum = val;
				} else if(val < minimum){
					minimum = val;
				}
			}
		}
	}
	int d = newtree.dim();
	std::cout << minimum << " " << maximum << " " << d << std::endl;

	std::ofstream ofile((char*) (std::string(filename) + ".txt").c_str());


	ofile << "reset" << std::endl;
	ofile << "set terminal " << terminal << " " << options << std::endl;
	ofile << "set output '" << filename << "." << std::string(terminal) << "'" << std::endl;
	//ofile << "plot '-' using 1:2" << std::endl;
	//ofile << "e" << std::endl << std::endl;
	ofile << "set multiplot " << std::endl;
	ofile << "set origin 0,0" << std::endl;
	ofile << "set size " << width << ", " << height << std::endl;
	ofile << "unset border" << std::endl;
	
	GeneralTree<double> xpos(newtree.getGeneralTree());
	GeneralTree<double> ypos(newtree.getGeneralTree());
	double smallheight = height/(1.0*newtree.dim());
	int maxlevel = newtree.getGeneralTree().end().getNode()->level();
	double smallwidth = width/(1.25*maxlevel);
	ofile << "unset xtics" << std::endl;
	ofile << "unset ytics" << std::endl;
	ofile << std::endl << "set xrange [0:"<<width<<"]"<< std::endl;
	ofile << std::endl << "set yrange [0:"<< height << "]" << std::endl;
	ofile << "set lmargin 0 " << std::endl;
	ofile << "set bmargin 0 " << std::endl;
	ofile << "set rmargin 0 " << std::endl;
	ofile << "set tmargin 0 " << std::endl;
	ofile << std::endl << "plot '-' using 1:2 with lines lt 1 lw 2 lc rgb 'black' notitle, \\";
	ofile << std::endl << "'' using 1:2 with lines lt 1 lw 2 lc rgb 'black' notitle";
    {

    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = newtree.getGeneralTree().end();
	GeneralTreeIterator<double> TITx = xpos.end();
	GeneralTreeIterator<double> TITy = ypos.end();
    
    for(;TIT <= newtree.getGeneralTree().begin(); TIT--,TITx--,TITy--){
	
		if(TIT.getNode()->isLeaf()){
			TITx.getNode()->setContent(1.25*smallwidth*(maxlevel - TIT.getNode()->level()));
			TITy.getNode()->setContent(height - (TIT.getNode()->getContent()->getIndex())[0]*smallheight*1.0);
		} else {
			TITx.getNode()->setContent(1.25*smallwidth*(maxlevel - TIT.getNode()->level()));
			TITy.getNode()->setContent(0.5*((*(TITy.getNode()->getfirstChild()->getContent()))+(*(TITy.getNode()->getlastChild()->getContent()))));
			//ofile << (*(TITx.getNode()->getfirstChild()->getContent()) + smallwidth) << " " << (*(TITy.getNode()->getfirstChild()->getContent()) + 0.5*smallheight) << std::endl;
			//ofile << *(TITx.getNode()->getContent()) << " " << (*(TITy.getNode()->getContent()) + 0.5*smallheight) << "," << std::endl;
		}
		if(! TIT.getNode()->isRoot() && !TIT.getNode()->isLeaf()){
			ofile << ", \\" << std::endl << "'' using 1:2 with lines lt 1 lw 2 lc rgb 'black' notitle";
			ofile << ", \\" << std::endl << "'' using 1:2 with lines lt 1 lw 2 lc rgb 'black' notitle";
		}

	
	}
    }
	ofile << std::endl;

    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = newtree.getGeneralTree().end();
	GeneralTreeIterator<double> TITx = xpos.end();
	GeneralTreeIterator<double> TITy = ypos.end();
	for(;TIT <= newtree.getGeneralTree().begin(); TIT--,TITx--,TITy--){
	
		if(!TIT.getNode()->isLeaf()){
			ofile << (*(TITx.getNode()->getfirstChild()->getContent()) + smallwidth) << " " << (*(TITy.getNode()->getfirstChild()->getContent()) + 0.5*smallheight) << std::endl;
			ofile << *(TITx.getNode()->getContent()) << " " << (*(TITy.getNode()->getContent()) + 0.5*smallheight)  << std::endl;
			ofile << "e" << std::endl;
			ofile << (*(TITx.getNode()->getlastChild()->getContent()) + smallwidth) << " " << (*(TITy.getNode()->getlastChild()->getContent()) + 0.5*smallheight) << std::endl;
			ofile << *(TITx.getNode()->getContent()) << " " << (*(TITy.getNode()->getContent()) + 0.5*smallheight)  << std::endl;
			ofile << "e" << std::endl;
		}
	}
    }

	ofile << "set colorbox user origin " << width - 0.075 << ",0.1 size 0.05," << height - 0.2 << std::endl;
	ofile << "set cbrange [" << minimum << ":" << maximum << "]" << std::endl;
	ofile << "show colorbox" << std::endl;

    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = newtree.getGeneralTree().end();
	GeneralTreeIterator<double> TITx = xpos.end();
	GeneralTreeIterator<double> TITy = ypos.end();
	for(;TIT <= newtree.getGeneralTree().begin(); TIT--,TITx--,TITy--){
	
		if(!TIT.getNode()->isLeaf()){
			ofile << "set label '$\\{";
				for(int v = 0; v < TIT.getNode()->getContent()->getIndex().length(); ++v){
					if(v!= 0){
						ofile << "," << TIT.getNode()->getContent()->getIndex()[v];
					} else {
						ofile << TIT.getNode()->getContent()->getIndex()[v];
					}
				}
			ofile << "\\}$' at " << (*TITx.getNode()->getContent())+0.4*smallwidth << "," << (*TITy.getNode()->getContent()) + 1.25*smallheight<< " center" << std::endl;
		} else {
				ofile << "set label '$\\{";
				for(int v = 0; v < TIT.getNode()->getContent()->getIndex().length(); ++v){
					if(v!= 0){
						ofile << "," << TIT.getNode()->getContent()->getIndex()[v];
					} else {
						ofile << TIT.getNode()->getContent()->getIndex()[v];
					}
				}
			ofile << "\\}$' at " << (*TITx.getNode()->getContent()) + smallwidth*0.4 << "," << (*TITy.getNode()->getContent()) + 0.85*smallheight << " center" << std::endl;
		}
	}
    }
	ofile << "plot '-' using 1:2:3 with dots palette notitle" << std::endl;
	//ofile << *(xpos.root->getContent()) + smallwidth/2. << " " << *(ypos.root->getContent()) + smallheight/2. << " " << minimum << std::endl;
	//ofile << *(xpos.root->getContent()) + smallwidth/2. << " " << *(ypos.root->getContent()) + smallheight/2. << " " << maximum << std::endl;
	ofile << "e" << std::endl << std::endl;

	ofile << "unset label" << std::endl;	
	ofile << "unset key" << std::endl;
	
    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = newtree.getGeneralTree().end();
	GeneralTreeIterator<double> TITx = xpos.end();
	GeneralTreeIterator<double> TITy = ypos.end();
    for(;TIT <= newtree.getGeneralTree().begin(); TIT--,TITx--,TITy--){

			UorB = TIT.getNode()->getContent()->getUorB();
			int numel = TIT.getNode()->getContent()->getNumRows();
			int lcnumel = TIT.getNode()->getContent()->getLeftChildNumRows();
			int rcnumel = TIT.getNode()->getContent()->getRightChildNumRows();
			std::cout << std::endl << TIT.getNode()->getContent()->getIndex() << std::endl;

			if(TIT.getNode()->isLeaf()){
				ofile << std::endl << std::endl;
				ofile << std::endl << std::endl;
				
				
				ofile << "set origin " << *(TITx.getNode()->getContent()) << "," << *(TITy.getNode()->getContent())<< std::endl;
				ofile << "set size " << smallwidth*0.75 << "," << smallheight*0.75 << std::endl;
				ofile << "set pointsize " << pointsize << std::endl;
				
				//g1.cmd("unset border");
				ofile << "unset xtics" << std::endl;
				ofile << "unset ytics" << std::endl;
				ofile << "unset colorbox" << std::endl;
				ofile << "set border 15" << std::endl;
				ofile << "set xrange [0:" << UorB.numCols() + 1 << "]" << std::endl;
				ofile << "set yrange [0:" << UorB.numRows() + 1 << "]" << std::endl;
				ofile << "set lmargin 0.05 " << std::endl;
				ofile << "set bmargin 0.00 " << std::endl;
				ofile << "set rmargin 0.05 " << std::endl;
				ofile << "set tmargin 0.0 " << std::endl;
				if(abs(pointsize) < 0.00000000001){
					ofile << "plot '-' using 1:2:3 with dots palette notitle" << std::endl;
				} else {
					ofile << "plot '-' using 1:2:3 pt 9 palette notitle" << std::endl;
				}
				for(int i = 1; i <= UorB.numCols(); ++i){
					for(int j = 1; j<= UorB.numCols(); ++j){
						if(abs(UorB(i,j)) > 0.000000000000001){
							ofile <<  j << " " << i << " " << UorB(i,j) << std::endl;
						}
					}
				}
				ofile << "e" << std::endl;
			} else {
				ofile << std::endl << std::endl;
				

				ofile << "set origin " << *(TITx.getNode()->getContent()) << "," << *(TITy.getNode()->getContent()) - 0.5*smallheight  << std::endl;
				ofile << "set size " << smallwidth << ", " << 2*smallheight << std::endl;
				ofile << "set pointsize " << pointsize << std::endl;
				//ofile << "set border 31" << std::endl;
				
				ofile << "unset xtics" << std::endl;
				ofile << "unset ytics" << std::endl;
				ofile << "unset ztics" << std::endl;
				ofile << "unset colorbox" << std::endl;

				ofile << "set lmargin 0 " << std::endl;
				ofile << "set bmargin 0 " << std::endl;
				ofile << "set rmargin 0 " << std::endl;
				ofile << "set tmargin 0 " << std::endl;

				ofile << "set xrange [0:" << numel + 1 << "]" << std::endl;
				ofile << "set yrange [1:" << rcnumel << "]" << std::endl;
				ofile << "set zrange [0:" << lcnumel+1  << "]" << std::endl;
				if(abs(pointsize) < 0.00000000001){
					ofile << "splot '-' using 1:2:3:4 with dots palette notitle" << std::endl;
				} else {
					ofile << "splot '-' using 1:2:3:4 pt 9 palette notitle" << std::endl;
				}
				for(int i = 1; i <= numel; ++i){
					for(int j = 1; j<= rcnumel; ++j){
						for(int k = 1; k<=lcnumel; ++k){
							if(abs(UorB((i-1)*rcnumel + j,k) ) > 0.000000000000001){
								ofile << j << " " << i <<  " " << k << " " << UorB((i-1)*rcnumel + j,k) << std::endl;
							}
						}
					}
				}
				ofile << "e" << std::endl;
			}
	}
    }
	ofile << "unset multiplot" << std::endl;
	
	ofile.close();
	
	std::ifstream infile((char*) (std::string(filename) + ".txt").c_str());

	std::string STRING;
	Gnuplot::set_GNUPlotPath("C:\\Users\\Andi\\Documents\\Programmieren\\gnuplot\\gnuplot\\bin");
	Gnuplot g1 = Gnuplot("");
	while(!infile.eof()){
		getline(infile,STRING);
		if(STRING.length() > 0){
			g1.cmd(STRING);
		}
	}
	infile.close();

}

template <typename T>
void
HTuckerTree<T>::plot_sv(const char* filename, const char* terminal, const char* options, const double width, const double height, const double pointsize) const{
	DimensionIndex base(d);
	base.setValueAsc();
	plot_sv(base,filename,terminal,options,width,height,pointsize);
}

template <typename T>
void 
HTuckerTree<T>::plot_sv(const DimensionIndex &subknoten,const char* filename,const char* terminal, const char* options, const double width, const double height, const double pointsize) const{
	int low;
	T maximum,minimum;
	minimum = 1; 
	maximum = 0;
	GeneralTreeNode<HTuckerTreeNode<T> > *node;
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.begin(); TIT <= tree.end(); TIT++){
		if(TIT.getNode()->getContent()->getIndex() == subknoten){
			std::cout << subknoten << std::endl;
			low = TIT.getNode()->level();
			std::cout << low << std::endl;
			node = TIT.getNode();
			break;
		}
	
	}
	int high = tree.end().getNode()->level();
	int maxrank = max_rank();
	std::cout << "high = " << high << std::endl;
	
	T val;
	HTuckerTree<T> newtree = subtree(*this,subknoten[0],subknoten[subknoten.length()-1]);
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &UorB = newtree.getGeneralTree().root->getContent()->getUorB();
	

	int d = newtree.dim();
	

	std::ofstream ofile((char*) (std::string(filename) + ".txt").c_str());


	ofile << "reset" << std::endl;
	ofile << "set terminal " << terminal << " " << options << std::endl;
	ofile << "set output '" << filename << "." << std::string(terminal) << "'" << std::endl;
	//ofile << "plot '-' using 1:2" << std::endl;
	//ofile << "e" << std::endl << std::endl;
	ofile << "set multiplot " << std::endl;
	ofile << "set origin 0,0" << std::endl;
	ofile << "set size " << width << ", " << height << std::endl;
	ofile << "unset border" << std::endl;
	
	GeneralTree<double> xpos(newtree.getGeneralTree());
	GeneralTree<double> ypos(newtree.getGeneralTree());
	GeneralTree<flens::DenseVector<flens::Array<double> > > svtree(newtree.getGeneralTree());
	double smallheight = height/(1.0*newtree.dim());
	int maxlevel = newtree.getGeneralTree().end().getNode()->level();
	double smallwidth = width/(1.25*maxlevel);
	ofile << "unset xtics" << std::endl;
	ofile << "unset ytics" << std::endl;
	ofile << std::endl << "set xrange [0:"<<width<<"]"<< std::endl;
	ofile << std::endl << "set yrange [0:"<< height << "]" << std::endl;
	ofile << "set lmargin 0 " << std::endl;
	ofile << "set bmargin 0 " << std::endl;
	ofile << "set rmargin 0 " << std::endl;
	ofile << "set tmargin 0 " << std::endl;
	ofile << std::endl << "plot '-' using 1:2 with lines lt 1 lw 2 lc rgb 'black' notitle, \\";
	ofile << std::endl << "'' using 1:2 with lines lt 1 lw 2 lc rgb 'black' notitle";
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > U,V;
	flens::DenseVector<flens::Array<T> > s,rho;
	
    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = newtree.getGeneralTree().end();
	GeneralTreeIterator<double> TITx = xpos.end();
	GeneralTreeIterator<double> TITy = ypos.end();
	GeneralTreeIterator<flens::DenseVector<flens::Array<double > > > TITsv = svtree.end();
    for(;TIT <= newtree.getGeneralTree().begin(); TIT--,TITx--,TITy--,TITsv--){
	
		UorB = TIT.getNode()->getContent()->getUorB();
		flens::svd(UorB,s,U,V);
		flens::sort(s,rho);
		TITsv.getNode()->setContent(s);
		//std::cout << "store at node " << TIT.getNode()->getContent()->getIndex() << " values " << s << std::endl;
		if(minimum <= maximum){
		
			if(s(s.firstIndex()) > maximum){
				maximum = s(s.lastIndex());
			}
			if(s(s.lastIndex()) < minimum){
				minimum = s(s.firstIndex());
			}
		} else {
			maximum = s(s.lastIndex());
			minimum = s(s.firstIndex());
		}

		if(TIT.getNode()->isLeaf()){
			TITx.getNode()->setContent(1.25*smallwidth*(maxlevel - TIT.getNode()->level()));
			TITy.getNode()->setContent(height - (TIT.getNode()->getContent()->getIndex())[0]*smallheight*1.0);
		} else {
			TITx.getNode()->setContent(1.25*smallwidth*(maxlevel - TIT.getNode()->level()));
			TITy.getNode()->setContent(0.5*((*(TITy.getNode()->getfirstChild()->getContent()))+(*(TITy.getNode()->getlastChild()->getContent()))));
			//ofile << (*(TITx.getNode()->getfirstChild()->getContent()) + smallwidth) << " " << (*(TITy.getNode()->getfirstChild()->getContent()) + 0.5*smallheight) << std::endl;
			//ofile << *(TITx.getNode()->getContent()) << " " << (*(TITy.getNode()->getContent()) + 0.5*smallheight) << "," << std::endl;
		}
		if(! TIT.getNode()->isRoot() && !TIT.getNode()->isLeaf()){
			ofile << ", \\" << std::endl << "'' using 1:2 with lines lt 1 lw 2 lc rgb 'black' notitle";
			ofile << ", \\" << std::endl << "'' using 1:2 with lines lt 1 lw 2 lc rgb 'black' notitle";
		}

	
	}
    }
	ofile << std::endl;
    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = newtree.getGeneralTree().end();
	GeneralTreeIterator<double> TITx = xpos.end();
	GeneralTreeIterator<double> TITy = ypos.end();
	for(;TIT <= newtree.getGeneralTree().begin(); TIT--,TITx--,TITy--){
	
		if(!TIT.getNode()->isLeaf()){
			ofile << (*(TITx.getNode()->getfirstChild()->getContent()) + 0.9*smallwidth) << " " << (*(TITy.getNode()->getfirstChild()->getContent()) + 0.5*smallheight) << std::endl;
			ofile << *(TITx.getNode()->getContent())-0.05*smallwidth << " " << (*(TITy.getNode()->getContent()) + 0.5*smallheight)  << std::endl;
			ofile << "e" << std::endl;
			ofile << (*(TITx.getNode()->getlastChild()->getContent()) + 0.9*smallwidth) << " " << (*(TITy.getNode()->getlastChild()->getContent()) + 0.5*smallheight) << std::endl;
			ofile << *(TITx.getNode()->getContent())-0.05*smallwidth << " " << (*(TITy.getNode()->getContent()) + 0.5*smallheight)  << std::endl;
			ofile << "e" << std::endl;
		}
	}
    }

	ofile << "set colorbox user origin " << width - 0.075 << ",0.1 size 0.05," << height - 0.2 << std::endl;
	ofile << "set cbrange [" << minimum << ":" << maximum << "]" << std::endl;
	ofile << "show colorbox" << std::endl;
    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = newtree.getGeneralTree().end();
	GeneralTreeIterator<double> TITx = xpos.end();
	GeneralTreeIterator<double> TITy = ypos.end();
	for(;TIT <= newtree.getGeneralTree().begin(); TIT--,TITx--,TITy--){
	
			if(!TIT.getNode()->isRoot()){
				ofile << "set label '$\\{";
					for(int v = 0; v < TIT.getNode()->getContent()->getIndex().length(); ++v){
						if(v!= 0){
							ofile << "," << TIT.getNode()->getContent()->getIndex()[v];
						} else {
							ofile << TIT.getNode()->getContent()->getIndex()[v];
						}
					}
				ofile << "\\}$' at " << (*TITx.getNode()->getContent())+0.4*smallwidth << "," << (*TITy.getNode()->getContent()) + 0.95*smallheight<< " center" << std::endl;
			}
	}
    }
	ofile << "plot '-' using 1:2:3 with dots palette notitle" << std::endl;
	//ofile << *(xpos.root->getContent()) + smallwidth/2. << " " << *(ypos.root->getContent()) + smallheight/2. << " " << minimum << std::endl;
	//ofile << *(xpos.root->getContent()) + smallwidth/2. << " " << *(ypos.root->getContent()) + smallheight/2. << " " << maximum << std::endl;
	ofile << "e" << std::endl << std::endl;

	ofile << "unset label" << std::endl;	
	ofile << "unset key" << std::endl;
	
	flens::DenseVector<flens::Array<T> > &sref = *(svtree.root->getContent());
    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = newtree.getGeneralTree().end();
	GeneralTreeIterator<double> TITx = xpos.end();
	GeneralTreeIterator<double> TITy = ypos.end();
	GeneralTreeIterator<flens::DenseVector<flens::Array<double > > > TITsv = svtree.end();
	for(;TIT <= newtree.getGeneralTree().begin(); TIT--,TITx--,TITy--,TITsv--){

			if(!TIT.getNode()->isRoot()){
				

				sref = *(TITsv.getNode()->getContent());
				//std::cout << "at node " << TIT.getNode()->getContent()->getIndex() << "  recieve values " << sref << std::endl;
				int numel = TIT.getNode()->getContent()->getNumRows();
				int lcnumel = TIT.getNode()->getContent()->getLeftChildNumRows();
				int rcnumel = TIT.getNode()->getContent()->getRightChildNumRows();
				std::cout << std::endl << TIT.getNode()->getContent()->getIndex() << std::endl;

				ofile << std::endl << std::endl;
				ofile << std::endl << std::endl;
				
				
				ofile << "set origin " << *(TITx.getNode()->getContent()) << "," << *(TITy.getNode()->getContent())+0.125*smallheight<< std::endl;
				ofile << "set size " << smallwidth*0.75 << "," << smallheight*0.75 << std::endl;
				ofile << "set pointsize " << pointsize << std::endl;
				
				//g1.cmd("unset border");
				ofile << "unset xtics" << std::endl;
				ofile << "unset ytics" << std::endl;
				ofile << "unset colorbox" << std::endl;
				ofile << "set grid" << std::endl;
				ofile << "set border 15" << std::endl;
				ofile << "set xrange [1:" << maxrank << "]" << std::endl;
				ofile << "set yrange [0:" << sref(sref.lastIndex()) << "]" << std::endl;
				ofile << "set lmargin 0.0 " << std::endl;
				ofile << "set bmargin 0.00 " << std::endl;
				ofile << "set rmargin 0.0 " << std::endl;
				ofile << "set tmargin 0.0 " << std::endl;
				ofile << "plot '-' using 1:2 with lines palette notitle" << std::endl;
				
				for(int i = s.lastIndex(),j = 1; i >= s.firstIndex(); i--,j++){
					ofile << j << " " << sref(i) << std::endl;
				}
				ofile << "e" << std::endl;
			}
			
	}
    }
	ofile << "unset multiplot" << std::endl;
	
	ofile.close();
	
	std::ifstream infile((char*) (std::string(filename) + ".txt").c_str());

	std::string STRING;
	Gnuplot::set_GNUPlotPath("C:\\Users\\Andi\\Documents\\Programmieren\\gnuplot\\gnuplot\\bin");
	Gnuplot g1 = Gnuplot("");
	while(!infile.eof()){
		getline(infile,STRING);
		if(STRING.length() > 0){
			g1.cmd(STRING);
		}
	}
	infile.close();

}


template <typename T>
HTuckerTree<T>
reapproximate(HTuckerTree<T> & tree1, const double epsilon, const int l){
	DimensionIndex maxvals(tree1.dim());
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree1.getGeneralTree().end(); TIT >= tree1.getGeneralTree().begin(); TIT--){
		if(TIT.getNode()->isLeaf()){
			maxvals[(TIT.getNode()->getContent()->getIndex())[0]-1] = TIT.getNode()->getContent()->getUorB().numRows();
		}
	}
	DimensionIndex minvals(1,tree1.dim());
	tensor<HTuckerTree2Tensor<T> > httt(minvals,maxvals);
	httt.getTensorFunction().setTree(tree1);
	HTuckerTree<T> newtree(tree1); //changed 7.02.2013 from tree1.dim() to maintain the treestructure.
	
	newtree.approximate(httt,epsilon,l);
	return newtree;
}

template <typename T>
HTuckerTree<T>
reapproximate(HTuckerTree<T> & tree1, const int rank, const int l){
	DimensionIndex maxvals(tree1.dim());
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree1->getGeneralTree().end(); TIT >= tree1.getGeneralTree().begin(); TIT--){
		if(TIT.getNode()->isLeaf()){
			maxvals[(TIT.getNode()->getContent()->getIndex())[0]-1] = TIT.getNode()->getContent()->getUorB().numRows();
		}
	}
	DimensionIndex minvals(1,tree1.dim());
	tensor<HTuckerTree2Tensor<T> > httt(minvals,maxvals);
	httt.getTensorFunction().setTree(tree1);
	HTuckerTree<T> newtree(tree1); //changed 7.02.2013 from tree1.dim() to maintain the treestructure.
	newtree.approximate(httt,rank,l);
	return newtree;
}



template <typename T>
HTuckerTree<T>
add_truncate(const HTuckerTree<T> & tree1,
             const HTuckerTree<T> & tree2,
             const double eps)
{
    HTuckerTree<T> tree = tree1+tree2;

    using flens::_;

    typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Matrix;
    typedef flens::DenseVector<flens::Array<T> >                      Vec;
    typedef flens::DenseVector<flens::Array<int> >                    IntVec;

    GeneralTreeNode<HTuckerTreeNode<T> > *node, *nodem, *nodeleft, *noderight;
    GeneralTreeNode<HTuckerTreeNode<T> > *gramnode;
    HTuckerTreeNode<T> *parent;
    HTuckerTree<T> gram;
    gram.set_tree(tree);
    HTuckerTree<T> Mtree;
    Mtree.set_tree(tree);

    auto TIT      = tree.getGeneralTree().end();
    auto TITleft  = tree1.getGeneralTree().end();
    auto TITright = tree2.getGeneralTree().end();
    auto TITm     = Mtree.getGeneralTree().end();
    for (; TIT >= tree.getGeneralTree().begin(); TIT--, TITm--,
                                                 TITleft--, TITright--) {
        node      = TIT.getNode();
        nodem     = TITm.getNode();
        nodeleft  = TITleft.getNode();
        noderight = TITright.getNode();
        if(node->isLeaf()) {
            Matrix tmp;
            Matrix& U = node->getContent()->getUorB();
            flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                            1., U, U, 0., tmp);
            nodem->getContent()->setUorB(tmp);
        } else if (node->isInner()) {
            Matrix& Mr = nodem->getlastChild()->getContent()->getUorB();
            Matrix& Ml = nodem->getfirstChild()->getContent()->getUorB();

            Matrix& Bt     = node->getContent()->getUorB();
            Matrix& Bleft  = nodeleft->getContent()->getUorB();
            Matrix& Bright = noderight->getContent()->getUorB();

            int numel = node->getContent()->getNumRows();

            int rcnumelleft  = nodeleft->getContent()->getRightChildNumRows();
            int rcnumelright = noderight->getContent()->getRightChildNumRows();
            int lcnumelleft  = nodeleft->getContent()->getLeftChildNumRows();
            int lcnumelright = noderight->getContent()->getLeftChildNumRows();

            int numelleft  = nodeleft->getContent()->getNumRows();
            int numelright = noderight->getContent()->getNumRows();

            int nrleft  = Ml.numRows();
            int nrright = Mr.numRows();

            Matrix tensor(nrright*(numelleft+numelright), nrleft);
            Matrix tmp(std::max(rcnumelleft, rcnumelright), nrleft);
            /* Compute tensor product */
            for (int i=1; i<=numelleft; ++i) {
                flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                1.,
                                Bleft(_((i-1)*rcnumelleft+1, i*rcnumelleft),
                                      _),
                                Ml(_, _(1, lcnumelleft)),
                                0., tmp(_(1, rcnumelleft), _));
                flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                1.,
                                Mr(_, _(1, rcnumelleft)),
                                tmp(_(1, rcnumelleft), _),
                                0., tensor(_((i-1)*nrright+1, i*nrright), _));
            }

            for (int i=1; i<=numelright; ++i) {
                flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                1.,
                                Bright(_((i-1)*rcnumelright+1, i*rcnumelright),
                                       _),
                                Ml(_, _(lcnumelleft+1,
                                        lcnumelleft+lcnumelright)),
                                0., tmp(_(1, rcnumelright), _));
                flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                1.,
                                Mr(_, _(rcnumelleft+1,
                                        rcnumelleft+rcnumelright)),
                                tmp(_(1, rcnumelright), _),
                                0., tensor(_(numelleft*nrright+(i-1)*nrright+1,
                                             numelleft*nrright+i*nrright), _));
            }

            /* Compute B_t^H*tensor */
            Matrix B11((lcnumelleft+lcnumelright)*rcnumelleft, numelleft);
            Matrix B22((lcnumelleft+lcnumelright)*rcnumelright, numelright);

            Matrix Mt(numel, numel);
            Matrix tensornonblock(nrright*nrleft, numelleft+numelright);

            for (int i=1; i<=numelleft+numelright; ++i) {
                for (int j=1; j<=nrright; ++j) {
                    tensornonblock(_((j-1)*nrleft+1, j*nrleft), i) =
                    tensor((i-1)*nrright+j, _);
                }
            }

            for (int i=1; i<=numelleft; ++i) {
                for (int j=1; j<=rcnumelleft; ++j) {
                    B11(_((j-1)*(lcnumelleft+lcnumelright)+1,
                           j*(lcnumelleft+lcnumelright)), i) =
                    Bt((i-1)*(rcnumelleft+rcnumelright)+j, _);
                }
            }

            for (int i=1; i<=numelright; ++i) {
                for (int j=1; j<=rcnumelright; ++j) {
                    B22(_((j-1)*(lcnumelleft+lcnumelright)+1,
                           j*(lcnumelleft+lcnumelright)), i) =
                    Bt((rcnumelleft+rcnumelright)*numelleft+rcnumelleft+
                       (i-1)*(rcnumelleft+rcnumelright)+j, _);
                }
            }

            flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                            1., B11,
                            tensornonblock(_(1,
                            (lcnumelleft+lcnumelright)*rcnumelleft), _),
                            0., Mt(_(1, numelleft), _));
            flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                            1., B22,
                            tensornonblock(
                            _((lcnumelleft+lcnumelright)*rcnumelleft+1,
                            (lcnumelleft+lcnumelright)*
                            (rcnumelleft+rcnumelright)), _),
                            0., Mt(_(numelleft+1, numelleft+numelright), _));

            nodem->getContent()->setUorB(Mt);
        } else {
            Matrix& Mr = nodem->getlastChild()->getContent()->getUorB();
            Matrix& Ml = nodem->getfirstChild()->getContent()->getUorB();

            Matrix& Bt     = node->getContent()->getUorB();
            Matrix& Bleft  = nodeleft->getContent()->getUorB();
            Matrix& Bright = noderight->getContent()->getUorB();

            int numel = node->getContent()->getNumRows();

            int rcnumelleft  = nodeleft->getContent()->getRightChildNumRows();
            int rcnumelright = noderight->getContent()->getRightChildNumRows();
            int lcnumelleft  = nodeleft->getContent()->getLeftChildNumRows();
            int lcnumelright = noderight->getContent()->getLeftChildNumRows();

            int nrleft  = Ml.numRows();
            int nrright = Mr.numRows();

            Matrix tensor(nrright, nrleft);
            /* Compute tensor product */
            Matrix tmp(rcnumelleft+rcnumelright, nrleft);
            flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                            1.,
                            Bleft,
                            Ml(_, _(1, lcnumelleft)),
                            0., tmp(_(1, rcnumelleft), _));
            flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                            1.,
                            Bright,
                            Ml(_, _(lcnumelleft+1, lcnumelleft+lcnumelright)),
                            0., tmp(_(rcnumelleft+1,
                                      rcnumelleft+rcnumelright), _));

            flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                            1.,
                            Mr(_, _(1, rcnumelleft)),
                            tmp(_(1, rcnumelleft), _),
                            0., tensor);
            flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                            1.,
                            Mr(_, _(rcnumelleft+1, rcnumelleft+rcnumelright)),
                            tmp(_(rcnumelleft+1, rcnumelleft+rcnumelright), _),
                            1., tensor);

            /* Compute B_t^H*tensor */
            Matrix Bnonblock((rcnumelleft+rcnumelright)*
                             (lcnumelleft+lcnumelright), 1);
            Matrix Mt(numel, numel);
            Matrix tensornonblock(nrright*nrleft, 1);

            for (int j=1; j<=nrright; ++j) {
                tensornonblock(_((j-1)*nrleft+1, j*nrleft), 1) =
                tensor(j, _);
            }

            for (int i=1; i<=rcnumelleft+rcnumelright; ++i) {
                Bnonblock(_((i-1)*(lcnumelleft+lcnumelright)+1,
                             i*(lcnumelleft+lcnumelright)), 1) =
                Bt(i, _);
            }

            flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                            1., Bnonblock,
                            tensornonblock,
                            0., Mt);

            nodem->getContent()->setUorB(Mt);
        }
    }

    TIT          = tree.getGeneralTree().begin();
    auto GramTIT = gram.getGeneralTree().begin();
    TITm         = Mtree.getGeneralTree().begin();
    TITleft      = tree1.getGeneralTree().begin();
    TITright     = tree2.getGeneralTree().begin();
    for(; TIT <= tree.getGeneralTree().end(); TIT++, TITm++, GramTIT++,
                                              TITleft++, TITright++){
        node      = TIT.getNode();
        nodem     = TITm.getNode();
        gramnode  = GramTIT.getNode();
        nodeleft  = TITleft.getNode();
        noderight = TITright.getNode();
        if(node->isRoot()){
            Matrix M(1,1);
            M = 1;
            gramnode->getContent()->setUorB(M);
        } else if (node->getParent()->isRoot()) {

            parent = node->getParent()->getContent();

            int numel         = node->getContent()->getNumRows();
            int parentrcnumel = parent->getRightChildNumRows();
            int parentlcnumel = parent->getLeftChildNumRows();

            int rcnumelleft  = nodeleft->getParent()->
                               getContent()->getRightChildNumRows();
            int rcnumelright = noderight->getParent()->
                               getContent()->getRightChildNumRows();
            int lcnumelleft  = nodeleft->getParent()->
                               getContent()->getLeftChildNumRows();
            int lcnumelright = noderight->getParent()->
                               getContent()->getLeftChildNumRows();

            Matrix& Bleft  = nodeleft->getParent()->getContent()->getUorB();
            Matrix& Bright = noderight->getParent()->getContent()->getUorB();

            Matrix& Mr = nodem->getParent()->
                         getlastChild()->getContent()->getUorB();
            Matrix& Ml = nodem->getParent()->
                         getfirstChild()->getContent()->getUorB();

            Matrix G(numel,numel);
            if (node==node->getParent()->getlastChild()) {
                /* G_tr */

                /* Compute tensor product */
                Matrix tensor(parentlcnumel, parentrcnumel);
                Vec tensorblock(lcnumelleft+lcnumelright);

                for (int i=1; i<=rcnumelleft; ++i) {
                    flens::blas::mv(cxxblas::NoTrans,
                                    1., Ml(_, _(1, lcnumelleft)),
                                    Bleft(i, _), 0.,
                                    tensorblock);
                    tensor(_, i) = tensorblock;
                }

                for (int i=1; i<=rcnumelright; ++i) {
                    flens::blas::mv(cxxblas::NoTrans,
                                    1., Ml(_, _(lcnumelleft+1,
                                                lcnumelleft+lcnumelright)),
                                    Bright(i, _), 0.,
                                    tensorblock);
                    tensor(_, rcnumelleft+i) = tensorblock;
                }

                /* Compute G_tr */
                flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                1., Bleft,
                                tensor(_(1, lcnumelleft), _),
                                0.,
                                G(_(1, rcnumelleft), _));
                flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                1., Bright,
                                tensor(_(lcnumelleft+1,
                                         lcnumelleft+lcnumelright), _),
                                0.,
                                G(_(rcnumelleft+1, rcnumelleft+rcnumelright),
                                  _));
            } else {
                /* G_tl */

                /* Compute tensor product */
                Matrix tensor(parentrcnumel, parentlcnumel);
                Vec tensorblock(rcnumelleft+rcnumelright);

                for (int i=1; i<=lcnumelleft; ++i) {
                    flens::blas::mv(cxxblas::NoTrans,
                                    1., Mr(_, _(1, rcnumelleft)),
                                    Bleft(_, i), 0.,
                                    tensorblock);
                    tensor(_, i) = tensorblock;
                }

                for (int i=1; i<=lcnumelright; ++i) {
                    flens::blas::mv(cxxblas::NoTrans,
                                    1., Mr(_, _(rcnumelleft+1,
                                                rcnumelleft+rcnumelright)),
                                    Bright(_, i), 0.,
                                    tensorblock);
                    tensor(_, lcnumelleft+i) = tensorblock;
                }

                /* Compute G_tl */
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., Bleft,
                                tensor(_(1, rcnumelleft), _),
                                0.,
                                G(_(1, lcnumelleft),
                                  _));
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., Bright,
                                tensor(_(rcnumelleft+1,
                                         rcnumelleft+rcnumelright), _),
                                0.,
                                G(_(lcnumelleft+1, lcnumelleft+lcnumelright),
                                  _));
            }
            gramnode->getContent()->setUorB(G);

        } else {
            parent = node->getParent()->getContent();

            int numel         = node->getContent()->getNumRows();
            int parentnumel   = parent->getNumRows();
            int parentrcnumel = parent->getRightChildNumRows();
            int parentlcnumel = parent->getLeftChildNumRows();
            Matrix &parentG   = gramnode->getParent()->getContent()->getUorB();

            int rcnumelleft  = nodeleft->getParent()->
                               getContent()->getRightChildNumRows();
            int rcnumelright = noderight->getParent()->
                               getContent()->getRightChildNumRows();
            int lcnumelleft  = nodeleft->getParent()->
                               getContent()->getLeftChildNumRows();
            int lcnumelright = noderight->getParent()->
                               getContent()->getLeftChildNumRows();

            int numelleft  = nodeleft->getParent()->getContent()->getNumRows();
            int numelright = noderight->getParent()->
                             getContent()->getNumRows();

            Matrix& Bleft  = nodeleft->getParent()->getContent()->getUorB();
            Matrix& Bright = noderight->getParent()->getContent()->getUorB();

            Matrix& Mr = nodem->getParent()->
                         getlastChild()->getContent()->getUorB();
            Matrix& Ml = nodem->getParent()->
                         getfirstChild()->getContent()->getUorB();

            Matrix G(numel,numel);
            if (node==node->getParent()->getlastChild()) {
                /* G_tr */
                Matrix B11((lcnumelleft+lcnumelright)*numelleft,
                            rcnumelleft);
                Matrix B22((lcnumelleft+lcnumelright)*numelright,
                            rcnumelright);
                for (int i=1; i<=numelleft; ++i) {
                    for (int j=1; j<=rcnumelleft; ++j) {
                        B11(_((i-1)*(lcnumelleft+lcnumelright)+1,
                              (i-1)*(lcnumelleft+lcnumelright)+
                               lcnumelleft), j) =
                        Bleft((i-1)*rcnumelleft+j, _);
                    }
                }

                for (int i=1; i<=numelright; ++i) {
                    for (int j=1; j<=rcnumelright; ++j) {
                        B22(_((i-1)*(lcnumelleft+lcnumelright)+
                               +lcnumelleft+1,
                              i*(lcnumelleft+lcnumelright)), j) =
                        Bright((i-1)*rcnumelright+j, _);
                    }
                }

                /* Compute tensor product */
                Matrix B31blockleft(lcnumelleft, numelleft);
                Matrix B31blockright(lcnumelright, numelright);
                Matrix tensor(parentnumel*parentlcnumel, parentrcnumel);
                Matrix tensorblock(lcnumelleft+lcnumelright, parentnumel);

                Matrix tmp(std::max(lcnumelleft, lcnumelright), parentnumel);
                for (int i=1; i<=rcnumelleft; ++i) {
                    for (int j=1; j<=numelleft; ++j) {
                        B31blockleft(_, j) =
                        Bleft((j-1)*rcnumelleft+i, _);
                    }
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                    1., B31blockleft,
                                    parentG(_, _(1, numelleft)),
                                    0., tmp(_(1, lcnumelleft), _));
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                    1., Ml(_, _(1, lcnumelleft)),
                                    tmp(_(1, lcnumelleft), _), 0.,
                                    tensorblock);
                    for (int j=1; j<=parentnumel; ++j) {
                        tensor(_((j-1)*parentlcnumel+1, j*parentlcnumel), i) =
                        tensorblock(_, j);
                    }
                }

                for (int i=1; i<=rcnumelright; ++i) {
                    for (int j=1; j<=numelright; ++j) {
                        B31blockright(_, j) =
                        Bright((j-1)*rcnumelright+i, _);
                    }
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                    1., B31blockright,
                                    parentG(_,
                                            _(numelleft+1,
                                              numelleft+numelright)),
                                    0., tmp(_(1, lcnumelright), _));
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                    1., Ml(_, _(lcnumelleft+1,
                                                lcnumelleft+lcnumelright)),
                                    tmp(_(1, lcnumelright), _), 0.,
                                    tensorblock);
                    for (int j=1; j<=parentnumel; ++j) {
                        tensor(_((j-1)*parentlcnumel+1, j*parentlcnumel),
                                rcnumelleft+i) =
                        tensorblock(_, j);
                    }
                }

                /* Compute G_tr */
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., B11,
                                tensor(_(1,
                                        (lcnumelleft+lcnumelright)*numelleft),
                                       _),
                                0.,
                                G(_(1, rcnumelleft),
                                  _));
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., B22,
                                tensor(_((lcnumelleft+lcnumelright)*
                                         numelleft+1,
                                         (lcnumelleft+lcnumelright)*
                                         (numelleft+numelright)), _),
                                0.,
                                G(_(rcnumelleft+1, rcnumelleft+rcnumelright),
                                  _));
            } else {
                /* G_tl */
                Matrix B11((rcnumelleft+rcnumelright)*numelleft,
                            lcnumelleft);
                Matrix B22((rcnumelleft+rcnumelright)*numelright,
                            lcnumelright);
                for (int i=1; i<=numelleft; ++i) {
                    for (int j=1; j<=lcnumelleft; ++j) {
                        B11(_((i-1)*(rcnumelleft+rcnumelright)+1,
                              (i-1)*(rcnumelleft+rcnumelright)+
                               rcnumelleft), j) =
                        Bleft(_((i-1)*rcnumelleft+1, i*rcnumelleft), j);
                    }
                }

                for (int i=1; i<=numelright; ++i) {
                    for (int j=1; j<=lcnumelright; ++j) {
                        B22(_((i-1)*(rcnumelleft+rcnumelright)+
                               +rcnumelleft+1,
                              i*(rcnumelleft+rcnumelright)), j) =
                        Bright(_((i-1)*rcnumelright+1, i*rcnumelright), j);
                    }
                }

                /* Compute tensor product */
                Matrix B32blockleft(rcnumelleft, numelleft);
                Matrix B32blockright(rcnumelright, numelright);
                Matrix tensor(parentnumel*parentrcnumel, parentlcnumel);
                Matrix tensorblock(rcnumelleft+rcnumelright, parentnumel);

                Matrix tmp(std::max(rcnumelleft, rcnumelright), parentnumel);
                for (int i=1; i<=lcnumelleft; ++i) {
                    for (int j=1; j<=numelleft; ++j) {
                        B32blockleft(_, j) =
                        Bleft(_((j-1)*rcnumelleft+1, j*rcnumelleft),
                               i);
                    }
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                    1., B32blockleft,
                                    parentG(_, _(1, numelleft)),
                                    0., tmp(_(1, rcnumelleft), _));
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                    1., Mr(_, _(1, rcnumelleft)),
                                    tmp(_(1, rcnumelleft), _),
                                    0.,
                                    tensorblock);
                    for (int j=1; j<=parentnumel; ++j) {
                        tensor(_((j-1)*parentrcnumel+1, j*parentrcnumel), i) =
                        tensorblock(_, j);
                    }
                }

                for (int i=1; i<=lcnumelright; ++i) {
                    for (int j=1; j<=numelright; ++j) {
                        B32blockright(_, j) =
                        Bright(_((j-1)*rcnumelright+1, j*rcnumelright), i);
                    }
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                    1., B32blockright,
                                    parentG(_,
                                            _(numelleft+1,
                                              numelleft+numelright)),
                                    0., tmp(_(1, rcnumelright), _));
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                    1., Mr(_, _(rcnumelleft+1,
                                                rcnumelleft+rcnumelright)),
                                    tmp(_(1, rcnumelright), _),
                                    0.,
                                    tensorblock);
                    for (int j=1; j<=parentnumel; ++j) {
                        tensor(_((j-1)*parentrcnumel+1, j*parentrcnumel),
                                lcnumelleft+i) =
                        tensorblock(_, j);
                    }
                }

                /* Compute G_tl */
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., B11,
                                tensor(_(1,
                                        (rcnumelleft+rcnumelright)*numelleft), _),
                                0.,
                                G(_(1, lcnumelleft),
                                  _));
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., B22,
                                tensor(_((rcnumelleft+rcnumelright)*
                                         numelleft+1,
                                         (rcnumelleft+rcnumelright)*
                                         (numelleft+numelright)), _),
                                0.,
                                G(_(lcnumelleft+1, lcnumelleft+lcnumelright),
                                  _));
            }
            gramnode->getContent()->setUorB(G);
        }
    }

    HTuckerTree<T> Rtree;
    Rtree.set_tree(tree);

    Vec s;
    Matrix U, Vt;

    std::vector<Vec> svals;
    std::vector<int> ranks;

    std::vector<Matrix> Us;
    std::vector<Matrix> Urefs;
    std::vector<Matrix> Qs;
    std::vector<Matrix> QRmats;

    TIT          = tree.getGeneralTree().end();
    TITleft      = tree1.getGeneralTree().end();
    TITright     = tree2.getGeneralTree().end();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITR = Rtree.getGeneralTree().end();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITgram = gram.getGeneralTree().end();

    GeneralTreeIterator<HTuckerTreeNode<T> > save;
    GeneralTreeIterator<HTuckerTreeNode<T> > saver;
    GeneralTreeIterator<HTuckerTreeNode<T> > savegram;
    GeneralTreeIterator<HTuckerTreeNode<T> > saveleft;
    GeneralTreeIterator<HTuckerTreeNode<T> > saveright;

    T val = eps*eps/(T) (tree.depth()-1);
    Vec eps2vec(tree.depth()-1);
    eps2vec.fill(val);

    int  depth    = TIT.getNode()->level();
    int  level    = depth;
    bool sorted   = false;
    bool newlevel = true;
    int  counter  = 0;
    int  length   = 0;
    for(;TIT >= tree.getGeneralTree().begin();){
        node      = TIT.getNode();
        nodeleft  = TITleft.getNode();
        noderight = TITright.getNode();

        if (node->level()==level) {

            if (newlevel) {
                save      = TIT;
                saver     = TITR;
                savegram  = TITgram;
                saveleft  = TITleft;
                saveright = TITright;

                Us.clear();
                Urefs.clear();
                Qs.clear();
                QRmats.clear();
                svals.clear();
                ranks.clear();

                length   = 0;
                sorted   = false;
                newlevel = false;
            }

            if (node->isLeaf()) {
                //compute QR-Decomposition
                //hier ist die Kopie notwendig, weil qrf die Matrix verändert!
                Matrix Uref = node->getContent()->getUorB();
                int m = Uref.numRows();
                int n = Uref.numCols();
                Vec tau;
                qrf(Uref,tau);

                Matrix Q(m,std::min(m,n));
                Q = Uref(_,_(1,std::min(m,n)));
                orgqr(Q, tau);


                Matrix &Gref = TITgram.getNode()->getContent()->getUorB();
                Matrix tmp1,tmp2,tmp3,tmp4;

                tmp1 = Gref;

                //hier berechnen wir R_t*G_t*R_t^T   => R_t entsprciht U_ref, wenn m > n dann ist die Matrix höher als Breit und enthält nur 0en, welche rausgeschmissen werden können
                if (m < n){
                    //dieser Fall ist komplizierter, weil die Matrix R zerlegt werden muss..
                    tmp2 = Gref(_,_(1,m));
                    tmp3 = Gref(_,_(m+1,n));
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,Uref(_,_(1,m)).upper(),tmp2);
                    tmp1.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0,tmp3,Uref(_,_(m+1,n)),0.0,tmp1);
                    tmp3.resize(0, 0);
                    tmp3 = tmp1 + tmp2;
                    tmp1.resize(0, 0);
                    tmp2.resize(0, 0);
                    tmp1 = tmp3(_(1,m),_);
                    tmp2 = tmp3(_(m+1,n),_);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,Uref(_,_(1,m)).upper(),tmp1);
                    tmp3.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Uref(_,_(m+1,n)),tmp2,0.0,tmp3);
                    tmp2.resize(0, 0);
                    tmp1 = tmp1 + tmp3;
                } else {
                    //d.h. in diesem Fall ist R in n x n und Q in m x n
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,Uref(_(1,n),_).upper(),tmp1);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,Uref(_(1,n),_).upper(),tmp1);
                }

                flens::svd(tmp1,s,U,Vt);

                Us.push_back(U);
                Urefs.push_back(Uref);
                Qs.push_back(Q);
                svals.push_back(s);
                length += s.length();
            } else if (node->isRoot()) {
                sorted = true;
                --level;
                continue;
            } else {
                //compute (Stl^T Rtl x Str^T Rtr) B_t
                int numel        = node->getContent()->getNumRows();

                Matrix &leftMat  = TITR.getNode()->getfirstChild()->getContent()->getUorB();
                Matrix &rightMat = TITR.getNode()->getlastChild()->getContent()->getUorB();
                Matrix QRmat(leftMat.numCols()*rightMat.numCols(),numel);
                Matrix tmp1,tmp2,tmp3,tmp4;
                Matrix& Bleft  = nodeleft->getContent()->getUorB();
                Matrix& Bright = noderight->getContent()->getUorB();

                int rcnumelleft  = nodeleft->
                                   getContent()->getRightChildNumRows();
                int rcnumelright = noderight->
                                   getContent()->getRightChildNumRows();
                int lcnumelleft  = nodeleft->
                                   getContent()->getLeftChildNumRows();
                int lcnumelright = noderight->
                                   getContent()->getLeftChildNumRows();

                int numelleft  = nodeleft->getContent()->getNumRows();
                int numelright = noderight->getContent()->getNumRows();

                int lr = leftMat.numCols(); // wegen leftMat^T
                int rr = rightMat.numCols();

                for (int i = 1; i <= numelleft; ++i) {
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans, 1.,
                                    Bleft(_((i-1)*rcnumelleft+1,
                                             i*rcnumelleft), _),
                                    leftMat(_(1, lcnumelleft), _), 0., tmp1);
                    flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans, 1.,
                                    rightMat(_(1, rcnumelleft), _),
                                    tmp1, 0.,tmp2);
                    //jetzt noch in die neue Matrix einsortieren...
                    for(int j = 1; j <= lr; ++j){
                        QRmat(_((j-1)*rr+1,j*rr),i) = tmp2(_,j);
                    }
                }

                tmp1.resize(0, 0);

                for (int i = 1; i <= numelright; ++i) {
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans, 1.,
                                    Bright(_((i-1)*rcnumelright+1,
                                              i*rcnumelright), _),
                                    leftMat(_(lcnumelleft+1,
                                              lcnumelleft+lcnumelright),_),
                                         0., tmp1);
                    flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans, 1.,
                                    rightMat(_(rcnumelleft+1,
                                               rcnumelleft+rcnumelright), _),
                                    tmp1, 0.,tmp2);
                    //jetzt noch in die neue Matrix einsortieren...
                    for(int j = 1; j <= lr; ++j){
                        QRmat(_((j-1)*rr+1,j*rr),numelleft+i) = tmp2(_,j);
                    }
                }

                // jetzt die QR Zerlegung anwenden:
                int m = QRmat.numRows();
                int n = QRmat.numCols();
                Vec tau;
                qrf(QRmat,tau);

                Matrix Q(m,std::min(m,n));
                Q = QRmat(_,_(1,std::min(m,n)));
                orgqr(Q,tau);

                // Dann wie oben R_t G_t R_t^T berechnen
                Matrix &Gref = TITgram.getNode()->getContent()->getUorB();

                tmp1.resize(0, 0);
                tmp1 = Gref;
                if (m < n){
                    //dieser Fall ist komplizierter, weil die Matrix R zerlegt werden muss..
                    tmp2.resize(0, 0);
                    tmp3.resize(0, 0);
                    tmp2 = Gref(_,_(1,m));
                    tmp3 = Gref(_,_(m+1,n));
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,QRmat(_,_(1,m)).upper(),tmp2);
                    tmp1.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0,tmp3,QRmat(_,_(m+1,n)),0.0,tmp1);
                    tmp3.resize(0, 0);
                    tmp3 = tmp1 + tmp2;
                    tmp1.resize(0, 0);
                    tmp1 = tmp3(_(1,m),_);
                    tmp2.resize(0, 0);
                    tmp2 = tmp3(_(m+1,n),_);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,QRmat(_,_(1,m)).upper(),tmp1);
                    tmp3.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,QRmat(_,_(m+1,n)),tmp2,0.0,tmp3);
                    tmp2.resize(0, 0);
                    tmp1 = tmp1 + tmp3;
                } else {
                    //d.h. in diesem Fall ist R in n x n und Q in m x n 
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,QRmat(_(1,n),_).upper(),tmp1);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,QRmat(_(1,n),_).upper(),tmp1);
                }

                //Svd wie oben
                flens::svd(tmp1,s,U,Vt);

                Us.push_back(U);
                Qs.push_back(Q);
                QRmats.push_back(QRmat);
                svals.push_back(s);
                length += s.length();
            }

            TIT--;
            TITR--;
            TITgram--;
            TITleft--;
            TITright--;

        } else if (!sorted) {
            TIT      = save;
            TITR     = saver;
            TITgram  = savegram;
            TITleft  = saveleft;
            TITright = saveright;
            sorted   = true;
            counter  = 0;

            Vec s_all(length);
            int off = 1;
            for (std::size_t i=0; i<svals.size(); ++i) {
                s_all(_(off,off+svals[i].length()-1)) = svals[i];
                off += svals[i].length();
            }


            IntVec rho(s_all.length());
            flens::sort(s_all, rho);
            double delta = 0.;
            double eps2  = eps2vec(depth-level+1);
            off          = 0;
            for (; off<s_all.length(); ++off) {
                delta += s_all(off+1);
                if (delta>eps2) break;
            }

            IntVec dels(svals.size());
            for (int i=1; i<=off; ++i) {
                int right = 0;
                for (std::size_t j=0; j<svals.size(); ++j) {
                    right += svals[j].length();
                    if (rho(i)<=right) {
                        ++dels(j+1);
                        break;
                    }
                }
            }

            ranks.resize(svals.size());
            for (std::size_t i=0; i<svals.size(); ++i) {
                int rank = svals[i].length() - dels(i+1);
                if (rank==0) rank = 1;
                ranks[i] = rank;
            }

            --level;
        } else {
            if (node->isLeaf()) {

                Matrix tmp1, tmp2, tmp3, tmp4;
                Matrix& Uref = Urefs[counter];
                Matrix& U    = Us[counter];
                Matrix& Q    = Qs[counter];
                int rank     = ranks[counter];
                int m        = node->getContent()->getUorB().numRows();
                int n        = node->getContent()->getUorB().numCols();;

                if ( m < n){
                    tmp2 = U(_,_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank)));
                    tmp3 = tmp2;
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    Uref(_,_(1,m)).upper(),tmp2);
                    flens::blas::mm(cxxblas::Trans,
                                    cxxblas::NoTrans, 1.0,
                                    Uref(_,_(m+1,n)),tmp3,0.0,tmp4);
                    tmp1.resize(0, 0);
                    tmp1 = Matrix(tmp2.numRows()+tmp4.numRows(),tmp2.numCols());
                    tmp1(_(1,tmp2.numRows()),_) = tmp2;
                    tmp1(_(tmp2.numRows()+1,tmp1.numRows()),_) = tmp4;
                } else {
                    tmp1 = U
                           (_(1,n),_(1,std::min(U.numCols(),
                                                (FLENS_DEFAULT_INDEXTYPE)rank)));
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    Uref(_(1,n),_).upper(),tmp1);

                }

                TITR.getNode()->getContent()->setUorB(tmp1);

                tmp1.resize(0, 0);
                flens::blas::mm(cxxblas::NoTrans,
                                cxxblas::NoTrans, 1.0,
                                Q,
                                U
                                (_,_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank))),0.0,tmp1);
                node->getContent()->setUorB(tmp1);
                node->getContent()->setNumRows(tmp1.numCols());
            } else if (node->isRoot()) {
                Matrix &leftMat  = TITR.getNode()->getfirstChild()->getContent()->getUorB();
                Matrix &rightMat = TITR.getNode()->getlastChild()->getContent()->getUorB();
                Matrix &Bref     = node->getContent()->getUorB();

                int rcnumelleft  = nodeleft->
                                   getContent()->getRightChildNumRows();
                int rcnumelright = noderight->
                                   getContent()->getRightChildNumRows();
                int lcnumelleft  = nodeleft->
                                   getContent()->getLeftChildNumRows();
                int lcnumelright = noderight->
                                   getContent()->getLeftChildNumRows();

                Matrix tmp1(rcnumelleft+rcnumelright, leftMat.numCols()), tmp2;
                //numel sollte hier immer 1 sein. Ansonsten genau das selbe wie weiter unten...

                flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans, 1.,
                                Bref(_(1, rcnumelleft), _(1, lcnumelleft)),
                                leftMat(_(1, lcnumelleft), _),
                                0., tmp1(_(1, rcnumelleft), _));
                flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans, 1.,
                                Bref(_(rcnumelleft+1,
                                       rcnumelleft+rcnumelright),
                                     _(lcnumelleft+1,
                                       lcnumelleft+lcnumelright)),
                                leftMat(_(lcnumelleft+1,
                                          lcnumelleft+lcnumelright), _),
                                0., tmp1(_(rcnumelleft+1,
                                           rcnumelleft+rcnumelright), _));

                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans, 1.,
                                rightMat(_(1, rcnumelleft), _),
                                tmp1(_(1, rcnumelleft), _),
                                0.,tmp2);
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans, 1.,
                                rightMat(_(rcnumelleft+1,
                                           rcnumelleft+rcnumelright), _),
                                tmp1(_(rcnumelleft+1,
                                       rcnumelleft+rcnumelright), _),
                                1., tmp2);

                //einsortiert werden muss hier nicht, da wir nicht mit QR weiterverarbeiten ist dies schon die Blockstruktur!

                node->getContent()->setUorB(tmp2);
                node->getContent()->setNumRows(1);
                node->getContent()->setLeftChildNumRows(tmp2.numCols());
                node->getContent()->setRightChildNumRows(tmp2.numRows());

            } else {
                int numel        = node->getContent()->getNumRows();
                Matrix &leftMat  = TITR.getNode()->getfirstChild()
                                   ->getContent()->getUorB();
                Matrix &rightMat = TITR.getNode()->getlastChild()
                                   ->getContent()->getUorB();

                Matrix tmp1, tmp2, tmp3, tmp4;
                Matrix& QRmat = QRmats[counter];
                Matrix& U     = Us[counter];
                Matrix& Q     = Qs[counter];
                int rank      = ranks[counter];
                int m         = leftMat.numCols()*rightMat.numCols();
                int n         = numel;

                if ( m < n){
                    tmp2 = U(_,_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank)));
                    tmp3 = tmp2;
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    QRmat(_,_(1,m)).upper(),tmp2);
                    flens::blas::mm(cxxblas::Trans,
                                    cxxblas::NoTrans, 1.0,
                                    QRmat(_,_(m+1,n)),tmp3,0.0,tmp4);
                    tmp1 = Matrix(tmp2.numRows()+tmp4.numRows(),tmp2.numCols());
                    tmp1(_(1,tmp2.numRows()),_) = tmp2;
                    tmp1(_(tmp2.numRows()+1,tmp1.numRows()),_) = tmp4;
                } else {
                    tmp1 = U(_(1,n),_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank)));
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    QRmat(_(1,n),_).upper(),tmp1);

                }

                TITR.getNode()->getContent()->setUorB(tmp1);

                //Das hier ist Q_t*S_t, die muss nun zurück ins Blockformat!!!
                tmp1.resize(0, 0);
                flens::blas::mm(cxxblas::NoTrans,
                                cxxblas::NoTrans, 1.0,
                                Q,
                                U(_,_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank))),0.0,tmp1);
                int newnumel   = tmp1.numCols();
                int newlcnumel = node->getfirstChild()->getContent()->getNumRows();
                int newrcnumel = node->getlastChild()->getContent()->getNumRows();
                Matrix newB(newnumel*newrcnumel,newlcnumel);

                for(int i = 1; i <= newnumel; ++i){
                    for(int j = 1; j <= newlcnumel; ++j){
                        newB(_((i-1)*newrcnumel+1,i*newrcnumel),j) = tmp1(_((j-1)*newrcnumel + 1,j*newrcnumel) ,i);
                    }
                }

                node->getContent()->setUorB(newB);
                node->getContent()->setNumRows(newnumel);
                node->getContent()->setLeftChildNumRows(newlcnumel);
                node->getContent()->setRightChildNumRows(newrcnumel);
            }

            TIT--;
            TITR--;
            TITgram--;
            TITleft--;
            TITright--;
            ++counter;
            newlevel = true;
        }
    }

    return tree;
}


template <typename T>
HTuckerTree<T>
add_truncate_soft(const HTuckerTree<T> & tree1,
                  const HTuckerTree<T> & tree2,
                  const double eps)
{
    HTuckerTree<T> tree = tree1+tree2;

    using flens::_;

    typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Matrix;
    typedef flens::DenseVector<flens::Array<T> >                      Vec;

    GeneralTreeNode<HTuckerTreeNode<T> > *node, *nodem, *nodeleft, *noderight;
    GeneralTreeNode<HTuckerTreeNode<T> > *gramnode;
    HTuckerTreeNode<T> *parent;
    HTuckerTree<T> gram;
    gram.set_tree(tree);
    HTuckerTree<T> Mtree;
    Mtree.set_tree(tree);

    auto TIT      = tree.getGeneralTree().end();
    auto TITleft  = tree1.getGeneralTree().end();
    auto TITright = tree2.getGeneralTree().end();
    auto TITm     = Mtree.getGeneralTree().end();
    for (; TIT >= tree.getGeneralTree().begin(); TIT--, TITm--,
                                                 TITleft--, TITright--) {
        node      = TIT.getNode();
        nodem     = TITm.getNode();
        nodeleft  = TITleft.getNode();
        noderight = TITright.getNode();
        if(node->isLeaf()) {
            Matrix tmp;
            Matrix& U = node->getContent()->getUorB();
            flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                            1., U, U, 0., tmp);
            nodem->getContent()->setUorB(tmp);
        } else if (node->isInner()) {
            Matrix& Mr = nodem->getlastChild()->getContent()->getUorB();
            Matrix& Ml = nodem->getfirstChild()->getContent()->getUorB();

            Matrix& Bt     = node->getContent()->getUorB();
            Matrix& Bleft  = nodeleft->getContent()->getUorB();
            Matrix& Bright = noderight->getContent()->getUorB();

            int numel = node->getContent()->getNumRows();

            int rcnumelleft  = nodeleft->getContent()->getRightChildNumRows();
            int rcnumelright = noderight->getContent()->getRightChildNumRows();
            int lcnumelleft  = nodeleft->getContent()->getLeftChildNumRows();
            int lcnumelright = noderight->getContent()->getLeftChildNumRows();

            int numelleft  = nodeleft->getContent()->getNumRows();
            int numelright = noderight->getContent()->getNumRows();

            int nrleft  = Ml.numRows();
            int nrright = Mr.numRows();

            Matrix tensor(nrright*(numelleft+numelright), nrleft);
            Matrix tmp(std::max(rcnumelleft, rcnumelright), nrleft);
            /* Compute tensor product */
            for (int i=1; i<=numelleft; ++i) {
                flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                1.,
                                Bleft(_((i-1)*rcnumelleft+1, i*rcnumelleft),
                                      _),
                                Ml(_, _(1, lcnumelleft)),
                                0., tmp(_(1, rcnumelleft), _));
                flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                1.,
                                Mr(_, _(1, rcnumelleft)),
                                tmp(_(1, rcnumelleft), _),
                                0., tensor(_((i-1)*nrright+1, i*nrright), _));
            }

            for (int i=1; i<=numelright; ++i) {
                flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                1.,
                                Bright(_((i-1)*rcnumelright+1, i*rcnumelright),
                                       _),
                                Ml(_, _(lcnumelleft+1,
                                        lcnumelleft+lcnumelright)),
                                0., tmp(_(1, rcnumelright), _));
                flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                1.,
                                Mr(_, _(rcnumelleft+1,
                                        rcnumelleft+rcnumelright)),
                                tmp(_(1, rcnumelright), _),
                                0., tensor(_(numelleft*nrright+(i-1)*nrright+1,
                                             numelleft*nrright+i*nrright), _));
            }

            /* Compute B_t^H*tensor */
            Matrix B11((lcnumelleft+lcnumelright)*rcnumelleft, numelleft);
            Matrix B22((lcnumelleft+lcnumelright)*rcnumelright, numelright);

            Matrix Mt(numel, numel);
            Matrix tensornonblock(nrright*nrleft, numelleft+numelright);

            for (int i=1; i<=numelleft+numelright; ++i) {
                for (int j=1; j<=nrright; ++j) {
                    tensornonblock(_((j-1)*nrleft+1, j*nrleft), i) =
                    tensor((i-1)*nrright+j, _);
                }
            }

            for (int i=1; i<=numelleft; ++i) {
                for (int j=1; j<=rcnumelleft; ++j) {
                    B11(_((j-1)*(lcnumelleft+lcnumelright)+1,
                           j*(lcnumelleft+lcnumelright)), i) =
                    Bt((i-1)*(rcnumelleft+rcnumelright)+j, _);
                }
            }

            for (int i=1; i<=numelright; ++i) {
                for (int j=1; j<=rcnumelright; ++j) {
                    B22(_((j-1)*(lcnumelleft+lcnumelright)+1,
                           j*(lcnumelleft+lcnumelright)), i) =
                    Bt((rcnumelleft+rcnumelright)*numelleft+rcnumelleft+
                       (i-1)*(rcnumelleft+rcnumelright)+j, _);
                }
            }

            flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                            1., B11,
                            tensornonblock(_(1,
                            (lcnumelleft+lcnumelright)*rcnumelleft), _),
                            0., Mt(_(1, numelleft), _));
            flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                            1., B22,
                            tensornonblock(
                            _((lcnumelleft+lcnumelright)*rcnumelleft+1,
                            (lcnumelleft+lcnumelright)*
                            (rcnumelleft+rcnumelright)), _),
                            0., Mt(_(numelleft+1, numelleft+numelright), _));

            nodem->getContent()->setUorB(Mt);
        } else {
            Matrix& Mr = nodem->getlastChild()->getContent()->getUorB();
            Matrix& Ml = nodem->getfirstChild()->getContent()->getUorB();

            Matrix& Bt     = node->getContent()->getUorB();
            Matrix& Bleft  = nodeleft->getContent()->getUorB();
            Matrix& Bright = noderight->getContent()->getUorB();

            int numel = node->getContent()->getNumRows();

            int rcnumelleft  = nodeleft->getContent()->getRightChildNumRows();
            int rcnumelright = noderight->getContent()->getRightChildNumRows();
            int lcnumelleft  = nodeleft->getContent()->getLeftChildNumRows();
            int lcnumelright = noderight->getContent()->getLeftChildNumRows();

            int nrleft  = Ml.numRows();
            int nrright = Mr.numRows();

            Matrix tensor(nrright, nrleft);
            /* Compute tensor product */
            Matrix tmp(rcnumelleft+rcnumelright, nrleft);
            flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                            1.,
                            Bleft,
                            Ml(_, _(1, lcnumelleft)),
                            0., tmp(_(1, rcnumelleft), _));
            flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                            1.,
                            Bright,
                            Ml(_, _(lcnumelleft+1, lcnumelleft+lcnumelright)),
                            0., tmp(_(rcnumelleft+1,
                                      rcnumelleft+rcnumelright), _));

            flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                            1.,
                            Mr(_, _(1, rcnumelleft)),
                            tmp(_(1, rcnumelleft), _),
                            0., tensor);
            flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                            1.,
                            Mr(_, _(rcnumelleft+1, rcnumelleft+rcnumelright)),
                            tmp(_(rcnumelleft+1, rcnumelleft+rcnumelright), _),
                            1., tensor);

            /* Compute B_t^H*tensor */
            Matrix Bnonblock((rcnumelleft+rcnumelright)*
                             (lcnumelleft+lcnumelright), 1);
            Matrix Mt(numel, numel);
            Matrix tensornonblock(nrright*nrleft, 1);

            for (int j=1; j<=nrright; ++j) {
                tensornonblock(_((j-1)*nrleft+1, j*nrleft), 1) =
                tensor(j, _);
            }

            for (int i=1; i<=rcnumelleft+rcnumelright; ++i) {
                Bnonblock(_((i-1)*(lcnumelleft+lcnumelright)+1,
                             i*(lcnumelleft+lcnumelright)), 1) =
                Bt(i, _);
            }

            flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                            1., Bnonblock,
                            tensornonblock,
                            0., Mt);

            nodem->getContent()->setUorB(Mt);
        }
    }

    TIT          = tree.getGeneralTree().begin();
    auto GramTIT = gram.getGeneralTree().begin();
    TITm         = Mtree.getGeneralTree().begin();
    TITleft      = tree1.getGeneralTree().begin();
    TITright     = tree2.getGeneralTree().begin();
    for(; TIT <= tree.getGeneralTree().end(); TIT++, TITm++, GramTIT++,
                                              TITleft++, TITright++){
        node      = TIT.getNode();
        nodem     = TITm.getNode();
        gramnode  = GramTIT.getNode();
        nodeleft  = TITleft.getNode();
        noderight = TITright.getNode();
        if(node->isRoot()){
            Matrix M(1,1);
            M = 1;
            gramnode->getContent()->setUorB(M);
        } else if (node->getParent()->isRoot()) {

            parent = node->getParent()->getContent();

            int numel         = node->getContent()->getNumRows();
            int parentrcnumel = parent->getRightChildNumRows();
            int parentlcnumel = parent->getLeftChildNumRows();

            int rcnumelleft  = nodeleft->getParent()->
                               getContent()->getRightChildNumRows();
            int rcnumelright = noderight->getParent()->
                               getContent()->getRightChildNumRows();
            int lcnumelleft  = nodeleft->getParent()->
                               getContent()->getLeftChildNumRows();
            int lcnumelright = noderight->getParent()->
                               getContent()->getLeftChildNumRows();

            Matrix& Bleft  = nodeleft->getParent()->getContent()->getUorB();
            Matrix& Bright = noderight->getParent()->getContent()->getUorB();

            Matrix& Mr = nodem->getParent()->
                         getlastChild()->getContent()->getUorB();
            Matrix& Ml = nodem->getParent()->
                         getfirstChild()->getContent()->getUorB();

            Matrix G(numel,numel);
            if (node==node->getParent()->getlastChild()) {
                /* G_tr */

                /* Compute tensor product */
                Matrix tensor(parentlcnumel, parentrcnumel);
                Vec tensorblock(lcnumelleft+lcnumelright);

                for (int i=1; i<=rcnumelleft; ++i) {
                    flens::blas::mv(cxxblas::NoTrans,
                                    1., Ml(_, _(1, lcnumelleft)),
                                    Bleft(i, _), 0.,
                                    tensorblock);
                    tensor(_, i) = tensorblock;
                }

                for (int i=1; i<=rcnumelright; ++i) {
                    flens::blas::mv(cxxblas::NoTrans,
                                    1., Ml(_, _(lcnumelleft+1,
                                                lcnumelleft+lcnumelright)),
                                    Bright(i, _), 0.,
                                    tensorblock);
                    tensor(_, rcnumelleft+i) = tensorblock;
                }

                /* Compute G_tr */
                flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                1., Bleft,
                                tensor(_(1, lcnumelleft), _),
                                0.,
                                G(_(1, rcnumelleft), _));
                flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                1., Bright,
                                tensor(_(lcnumelleft+1,
                                         lcnumelleft+lcnumelright), _),
                                0.,
                                G(_(rcnumelleft+1, rcnumelleft+rcnumelright),
                                  _));
            } else {
                /* G_tl */

                /* Compute tensor product */
                Matrix tensor(parentrcnumel, parentlcnumel);
                Vec tensorblock(rcnumelleft+rcnumelright);

                for (int i=1; i<=lcnumelleft; ++i) {
                    flens::blas::mv(cxxblas::NoTrans,
                                    1., Mr(_, _(1, rcnumelleft)),
                                    Bleft(_, i), 0.,
                                    tensorblock);
                    tensor(_, i) = tensorblock;
                }

                for (int i=1; i<=lcnumelright; ++i) {
                    flens::blas::mv(cxxblas::NoTrans,
                                    1., Mr(_, _(rcnumelleft+1,
                                                rcnumelleft+rcnumelright)),
                                    Bright(_, i), 0.,
                                    tensorblock);
                    tensor(_, lcnumelleft+i) = tensorblock;
                }

                /* Compute G_tl */
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., Bleft,
                                tensor(_(1, rcnumelleft), _),
                                0.,
                                G(_(1, lcnumelleft),
                                  _));
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., Bright,
                                tensor(_(rcnumelleft+1,
                                         rcnumelleft+rcnumelright), _),
                                0.,
                                G(_(lcnumelleft+1, lcnumelleft+lcnumelright),
                                  _));
            }
            gramnode->getContent()->setUorB(G);

        } else {
            parent = node->getParent()->getContent();

            int numel         = node->getContent()->getNumRows();
            int parentnumel   = parent->getNumRows();
            int parentrcnumel = parent->getRightChildNumRows();
            int parentlcnumel = parent->getLeftChildNumRows();
            Matrix &parentG   = gramnode->getParent()->getContent()->getUorB();

            int rcnumelleft  = nodeleft->getParent()->
                               getContent()->getRightChildNumRows();
            int rcnumelright = noderight->getParent()->
                               getContent()->getRightChildNumRows();
            int lcnumelleft  = nodeleft->getParent()->
                               getContent()->getLeftChildNumRows();
            int lcnumelright = noderight->getParent()->
                               getContent()->getLeftChildNumRows();

            int numelleft  = nodeleft->getParent()->getContent()->getNumRows();
            int numelright = noderight->getParent()->
                             getContent()->getNumRows();

            Matrix& Bleft  = nodeleft->getParent()->getContent()->getUorB();
            Matrix& Bright = noderight->getParent()->getContent()->getUorB();

            Matrix& Mr = nodem->getParent()->
                         getlastChild()->getContent()->getUorB();
            Matrix& Ml = nodem->getParent()->
                         getfirstChild()->getContent()->getUorB();

            Matrix G(numel,numel);
            if (node==node->getParent()->getlastChild()) {
                /* G_tr */
                Matrix B11((lcnumelleft+lcnumelright)*numelleft,
                            rcnumelleft);
                Matrix B22((lcnumelleft+lcnumelright)*numelright,
                            rcnumelright);
                for (int i=1; i<=numelleft; ++i) {
                    for (int j=1; j<=rcnumelleft; ++j) {
                        B11(_((i-1)*(lcnumelleft+lcnumelright)+1,
                              (i-1)*(lcnumelleft+lcnumelright)+
                               lcnumelleft), j) =
                        Bleft((i-1)*rcnumelleft+j, _);
                    }
                }

                for (int i=1; i<=numelright; ++i) {
                    for (int j=1; j<=rcnumelright; ++j) {
                        B22(_((i-1)*(lcnumelleft+lcnumelright)+
                               +lcnumelleft+1,
                              i*(lcnumelleft+lcnumelright)), j) =
                        Bright((i-1)*rcnumelright+j, _);
                    }
                }

                /* Compute tensor product */
                Matrix B31blockleft(lcnumelleft, numelleft);
                Matrix B31blockright(lcnumelright, numelright);
                Matrix tensor(parentnumel*parentlcnumel, parentrcnumel);
                Matrix tensorblock(lcnumelleft+lcnumelright, parentnumel);

                Matrix tmp(std::max(lcnumelleft, lcnumelright), parentnumel);
                for (int i=1; i<=rcnumelleft; ++i) {
                    for (int j=1; j<=numelleft; ++j) {
                        B31blockleft(_, j) =
                        Bleft((j-1)*rcnumelleft+i, _);
                    }
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                    1., B31blockleft,
                                    parentG(_, _(1, numelleft)),
                                    0., tmp(_(1, lcnumelleft), _));
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                    1., Ml(_, _(1, lcnumelleft)),
                                    tmp(_(1, lcnumelleft), _), 0.,
                                    tensorblock);
                    for (int j=1; j<=parentnumel; ++j) {
                        tensor(_((j-1)*parentlcnumel+1, j*parentlcnumel), i) =
                        tensorblock(_, j);
                    }
                }

                for (int i=1; i<=rcnumelright; ++i) {
                    for (int j=1; j<=numelright; ++j) {
                        B31blockright(_, j) =
                        Bright((j-1)*rcnumelright+i, _);
                    }
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                    1., B31blockright,
                                    parentG(_,
                                            _(numelleft+1,
                                              numelleft+numelright)),
                                    0., tmp(_(1, lcnumelright), _));
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                    1., Ml(_, _(lcnumelleft+1,
                                                lcnumelleft+lcnumelright)),
                                    tmp(_(1, lcnumelright), _), 0.,
                                    tensorblock);
                    for (int j=1; j<=parentnumel; ++j) {
                        tensor(_((j-1)*parentlcnumel+1, j*parentlcnumel),
                                rcnumelleft+i) =
                        tensorblock(_, j);
                    }
                }

                /* Compute G_tr */
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., B11,
                                tensor(_(1,
                                        (lcnumelleft+lcnumelright)*numelleft),
                                       _),
                                0.,
                                G(_(1, rcnumelleft),
                                  _));
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., B22,
                                tensor(_((lcnumelleft+lcnumelright)*
                                         numelleft+1,
                                         (lcnumelleft+lcnumelright)*
                                         (numelleft+numelright)), _),
                                0.,
                                G(_(rcnumelleft+1, rcnumelleft+rcnumelright),
                                  _));
            } else {
                /* G_tl */
                Matrix B11((rcnumelleft+rcnumelright)*numelleft,
                            lcnumelleft);
                Matrix B22((rcnumelleft+rcnumelright)*numelright,
                            lcnumelright);
                for (int i=1; i<=numelleft; ++i) {
                    for (int j=1; j<=lcnumelleft; ++j) {
                        B11(_((i-1)*(rcnumelleft+rcnumelright)+1,
                              (i-1)*(rcnumelleft+rcnumelright)+
                               rcnumelleft), j) =
                        Bleft(_((i-1)*rcnumelleft+1, i*rcnumelleft), j);
                    }
                }

                for (int i=1; i<=numelright; ++i) {
                    for (int j=1; j<=lcnumelright; ++j) {
                        B22(_((i-1)*(rcnumelleft+rcnumelright)+
                               +rcnumelleft+1,
                              i*(rcnumelleft+rcnumelright)), j) =
                        Bright(_((i-1)*rcnumelright+1, i*rcnumelright), j);
                    }
                }

                /* Compute tensor product */
                Matrix B32blockleft(rcnumelleft, numelleft);
                Matrix B32blockright(rcnumelright, numelright);
                Matrix tensor(parentnumel*parentrcnumel, parentlcnumel);
                Matrix tensorblock(rcnumelleft+rcnumelright, parentnumel);

                Matrix tmp(std::max(rcnumelleft, rcnumelright), parentnumel);
                for (int i=1; i<=lcnumelleft; ++i) {
                    for (int j=1; j<=numelleft; ++j) {
                        B32blockleft(_, j) =
                        Bleft(_((j-1)*rcnumelleft+1, j*rcnumelleft),
                               i);
                    }
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                    1., B32blockleft,
                                    parentG(_, _(1, numelleft)),
                                    0., tmp(_(1, rcnumelleft), _));
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                    1., Mr(_, _(1, rcnumelleft)),
                                    tmp(_(1, rcnumelleft), _),
                                    0.,
                                    tensorblock);
                    for (int j=1; j<=parentnumel; ++j) {
                        tensor(_((j-1)*parentrcnumel+1, j*parentrcnumel), i) =
                        tensorblock(_, j);
                    }
                }

                for (int i=1; i<=lcnumelright; ++i) {
                    for (int j=1; j<=numelright; ++j) {
                        B32blockright(_, j) =
                        Bright(_((j-1)*rcnumelright+1, j*rcnumelright), i);
                    }
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                    1., B32blockright,
                                    parentG(_,
                                            _(numelleft+1,
                                              numelleft+numelright)),
                                    0., tmp(_(1, rcnumelright), _));
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                    1., Mr(_, _(rcnumelleft+1,
                                                rcnumelleft+rcnumelright)),
                                    tmp(_(1, rcnumelright), _),
                                    0.,
                                    tensorblock);
                    for (int j=1; j<=parentnumel; ++j) {
                        tensor(_((j-1)*parentrcnumel+1, j*parentrcnumel),
                                lcnumelleft+i) =
                        tensorblock(_, j);
                    }
                }

                /* Compute G_tl */
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., B11,
                                tensor(_(1,
                                        (rcnumelleft+rcnumelright)*numelleft), _),
                                0.,
                                G(_(1, lcnumelleft),
                                  _));
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., B22,
                                tensor(_((rcnumelleft+rcnumelright)*
                                         numelleft+1,
                                         (rcnumelleft+rcnumelright)*
                                         (numelleft+numelright)), _),
                                0.,
                                G(_(lcnumelleft+1, lcnumelleft+lcnumelright),
                                  _));
            }
            gramnode->getContent()->setUorB(G);
        }
    }

    HTuckerTree<T> Rtree;
    Rtree.set_tree(tree);

    Vec s;
    Matrix U, Vt;

    std::vector<Vec> svals;
    std::vector<Vec> scals;
    std::vector<int> ranks;

    std::vector<Matrix> Us;
    std::vector<Matrix> Urefs;
    std::vector<Matrix> Qs;
    std::vector<Matrix> QRmats;

    TIT          = tree.getGeneralTree().end();
    TITleft      = tree1.getGeneralTree().end();
    TITright     = tree2.getGeneralTree().end();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITR = Rtree.getGeneralTree().end();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITgram = gram.getGeneralTree().end();

    GeneralTreeIterator<HTuckerTreeNode<T> > save;
    GeneralTreeIterator<HTuckerTreeNode<T> > saver;
    GeneralTreeIterator<HTuckerTreeNode<T> > savegram;
    GeneralTreeIterator<HTuckerTreeNode<T> > saveleft;
    GeneralTreeIterator<HTuckerTreeNode<T> > saveright;

    T alpha2 = eps*eps;

    int  depth    = TIT.getNode()->level();
    int  level    = depth;
    bool sorted   = false;
    bool newlevel = true;
    int  counter  = 0;
    int  length   = 0;
    for(;TIT >= tree.getGeneralTree().begin();){
        node      = TIT.getNode();
        nodeleft  = TITleft.getNode();
        noderight = TITright.getNode();

        if (node->level()==level) {

            if (newlevel) {
                save      = TIT;
                saver     = TITR;
                savegram  = TITgram;
                saveleft  = TITleft;
                saveright = TITright;

                Us.clear();
                Urefs.clear();
                Qs.clear();
                QRmats.clear();
                svals.clear();
                scals.clear();
                ranks.clear();

                length   = 0;
                sorted   = false;
                newlevel = false;
            }

            if (node->isLeaf()) {
                //compute QR-Decomposition
                //hier ist die Kopie notwendig, weil qrf die Matrix verändert!
                Matrix Uref = node->getContent()->getUorB();
                int m = Uref.numRows();
                int n = Uref.numCols();
                Vec tau;
                qrf(Uref,tau);

                Matrix Q(m,std::min(m,n));
                Q = Uref(_,_(1,std::min(m,n)));
                orgqr(Q, tau);


                Matrix &Gref = TITgram.getNode()->getContent()->getUorB();
                Matrix tmp1,tmp2,tmp3,tmp4;

                tmp1 = Gref;

                //hier berechnen wir R_t*G_t*R_t^T   => R_t entsprciht U_ref, wenn m > n dann ist die Matrix höher als Breit und enthält nur 0en, welche rausgeschmissen werden können
                if (m < n){
                    //dieser Fall ist komplizierter, weil die Matrix R zerlegt werden muss..
                    tmp2 = Gref(_,_(1,m));
                    tmp3 = Gref(_,_(m+1,n));
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,Uref(_,_(1,m)).upper(),tmp2);
                    tmp1.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0,tmp3,Uref(_,_(m+1,n)),0.0,tmp1);
                    tmp3.resize(0, 0);
                    tmp3 = tmp1 + tmp2;
                    tmp1.resize(0, 0);
                    tmp2.resize(0, 0);
                    tmp1 = tmp3(_(1,m),_);
                    tmp2 = tmp3(_(m+1,n),_);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,Uref(_,_(1,m)).upper(),tmp1);
                    tmp3.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Uref(_,_(m+1,n)),tmp2,0.0,tmp3);
                    tmp2.resize(0, 0);
                    tmp1 = tmp1 + tmp3;
                } else {
                    //d.h. in diesem Fall ist R in n x n und Q in m x n
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,Uref(_(1,n),_).upper(),tmp1);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,Uref(_(1,n),_).upper(),tmp1);
                }

                flens::svd(tmp1,s,U,Vt);

                Us.push_back(U);
                Urefs.push_back(Uref);
                Qs.push_back(Q);
                svals.push_back(s);
                scals.push_back(s);
                length += s.length();
            } else if (node->isRoot()) {
                sorted = true;
                --level;
                continue;
            } else {
                //compute (Stl^T Rtl x Str^T Rtr) B_t
                int numel        = node->getContent()->getNumRows();

                Matrix &leftMat  = TITR.getNode()->getfirstChild()->getContent()->getUorB();
                Matrix &rightMat = TITR.getNode()->getlastChild()->getContent()->getUorB();
                Matrix QRmat(leftMat.numCols()*rightMat.numCols(),numel);
                Matrix tmp1,tmp2,tmp3,tmp4;
                Matrix& Bleft  = nodeleft->getContent()->getUorB();
                Matrix& Bright = noderight->getContent()->getUorB();

                int rcnumelleft  = nodeleft->
                                   getContent()->getRightChildNumRows();
                int rcnumelright = noderight->
                                   getContent()->getRightChildNumRows();
                int lcnumelleft  = nodeleft->
                                   getContent()->getLeftChildNumRows();
                int lcnumelright = noderight->
                                   getContent()->getLeftChildNumRows();

                int numelleft  = nodeleft->getContent()->getNumRows();
                int numelright = noderight->getContent()->getNumRows();

                int lr = leftMat.numCols(); // wegen leftMat^T
                int rr = rightMat.numCols();

                for (int i = 1; i <= numelleft; ++i) {
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans, 1.,
                                    Bleft(_((i-1)*rcnumelleft+1,
                                             i*rcnumelleft), _),
                                    leftMat(_(1, lcnumelleft), _), 0., tmp1);
                    flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans, 1.,
                                    rightMat(_(1, rcnumelleft), _),
                                    tmp1, 0.,tmp2);
                    //jetzt noch in die neue Matrix einsortieren...
                    for(int j = 1; j <= lr; ++j){
                        QRmat(_((j-1)*rr+1,j*rr),i) = tmp2(_,j);
                    }
                }

                tmp1.resize(0, 0);

                for (int i = 1; i <= numelright; ++i) {
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans, 1.,
                                    Bright(_((i-1)*rcnumelright+1,
                                              i*rcnumelright), _),
                                    leftMat(_(lcnumelleft+1,
                                              lcnumelleft+lcnumelright),_),
                                         0., tmp1);
                    flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans, 1.,
                                    rightMat(_(rcnumelleft+1,
                                               rcnumelleft+rcnumelright), _),
                                    tmp1, 0.,tmp2);
                    //jetzt noch in die neue Matrix einsortieren...
                    for(int j = 1; j <= lr; ++j){
                        QRmat(_((j-1)*rr+1,j*rr),numelleft+i) = tmp2(_,j);
                    }
                }

                // jetzt die QR Zerlegung anwenden:
                int m = QRmat.numRows();
                int n = QRmat.numCols();
                Vec tau;
                qrf(QRmat,tau);

                Matrix Q(m,std::min(m,n));
                Q = QRmat(_,_(1,std::min(m,n)));
                orgqr(Q,tau);

                // Dann wie oben R_t G_t R_t^T berechnen
                Matrix &Gref = TITgram.getNode()->getContent()->getUorB();

                tmp1.resize(0, 0);
                tmp1 = Gref;
                if (m < n){
                    //dieser Fall ist komplizierter, weil die Matrix R zerlegt werden muss..
                    tmp2.resize(0, 0);
                    tmp3.resize(0, 0);
                    tmp2 = Gref(_,_(1,m));
                    tmp3 = Gref(_,_(m+1,n));
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,QRmat(_,_(1,m)).upper(),tmp2);
                    tmp1.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0,tmp3,QRmat(_,_(m+1,n)),0.0,tmp1);
                    tmp3.resize(0, 0);
                    tmp3 = tmp1 + tmp2;
                    tmp1.resize(0, 0);
                    tmp1 = tmp3(_(1,m),_);
                    tmp2.resize(0, 0);
                    tmp2 = tmp3(_(m+1,n),_);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,QRmat(_,_(1,m)).upper(),tmp1);
                    tmp3.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,QRmat(_,_(m+1,n)),tmp2,0.0,tmp3);
                    tmp2.resize(0, 0);
                    tmp1 = tmp1 + tmp3;
                } else {
                    //d.h. in diesem Fall ist R in n x n und Q in m x n 
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,QRmat(_(1,n),_).upper(),tmp1);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,QRmat(_(1,n),_).upper(),tmp1);
                }

                //Svd wie oben
                flens::svd(tmp1,s,U,Vt);

                Us.push_back(U);
                Qs.push_back(Q);
                QRmats.push_back(QRmat);
                svals.push_back(s);
                scals.push_back(s);
                length += s.length();
            }

            TIT--;
            TITR--;
            TITgram--;
            TITleft--;
            TITright--;

        } else if (!sorted) {
            TIT      = save;
            TITR     = saver;
            TITgram  = savegram;
            TITleft  = saveleft;
            TITright = saveright;
            sorted   = true;
            counter  = 0;

            ranks.resize(svals.size());
            for (std::size_t i=0; i<svals.size(); ++i) {
                ranks[i] = 0;
                for (int j=svals[i].length(); j>=1; --j) {
                    if (svals[i](j)>alpha2) break;
                    else ++ranks[i];
                }
                ranks[i] = svals[i].length()-ranks[i];
                ranks[i] = std::max(ranks[i], 1);
            }

            for (std::size_t i=0; i<svals.size(); ++i) {
                for (int j=1; j<=ranks[i]; ++j) {
                    auto sigma = std::sqrt(svals[i](j));
                    scals[i](j) = (sigma-eps)/sigma;
                }
            }

            --level;
        } else {
            if (node->isLeaf()) {

                Matrix tmp1, tmp2, tmp3, tmp4;
                Matrix& Uref = Urefs[counter];
                Matrix& U    = Us[counter];
                Matrix& Q    = Qs[counter];
                int rank     = ranks[counter];
                int m        = node->getContent()->getUorB().numRows();
                int n        = node->getContent()->getUorB().numCols();;
                Vec ks       = scals[counter];

                for (int i=1; i<=rank; ++i) {
                    flens::blas::scal(ks(i), U(_, i));
                }

                if ( m < n){
                    tmp2 = U(_,_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank)));
                    tmp3 = tmp2;
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    Uref(_,_(1,m)).upper(),tmp2);
                    flens::blas::mm(cxxblas::Trans,
                                    cxxblas::NoTrans, 1.0,
                                    Uref(_,_(m+1,n)),tmp3,0.0,tmp4);
                    tmp1.resize(0, 0);
                    tmp1 = Matrix(tmp2.numRows()+tmp4.numRows(),tmp2.numCols());
                    tmp1(_(1,tmp2.numRows()),_) = tmp2;
                    tmp1(_(tmp2.numRows()+1,tmp1.numRows()),_) = tmp4;
                } else {
                    tmp1 = U
                           (_(1,n),_(1,std::min(U.numCols(),
                                                (FLENS_DEFAULT_INDEXTYPE)rank)));
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    Uref(_(1,n),_).upper(),tmp1);

                }

                TITR.getNode()->getContent()->setUorB(tmp1);

                tmp1.resize(0, 0);
                flens::blas::mm(cxxblas::NoTrans,
                                cxxblas::NoTrans, 1.0,
                                Q,
                                U
                                (_,_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank))),0.0,tmp1);
                node->getContent()->setUorB(tmp1);
                node->getContent()->setNumRows(tmp1.numCols());
            } else if (node->isRoot()) {
                Matrix &leftMat  = TITR.getNode()->getfirstChild()->getContent()->getUorB();
                Matrix &rightMat = TITR.getNode()->getlastChild()->getContent()->getUorB();
                Matrix &Bref     = node->getContent()->getUorB();

                int rcnumelleft  = nodeleft->
                                   getContent()->getRightChildNumRows();
                int rcnumelright = noderight->
                                   getContent()->getRightChildNumRows();
                int lcnumelleft  = nodeleft->
                                   getContent()->getLeftChildNumRows();
                int lcnumelright = noderight->
                                   getContent()->getLeftChildNumRows();

                Matrix tmp1(rcnumelleft+rcnumelright, leftMat.numCols()), tmp2;
                //numel sollte hier immer 1 sein. Ansonsten genau das selbe wie weiter unten...

                flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans, 1.,
                                Bref(_(1, rcnumelleft), _(1, lcnumelleft)),
                                leftMat(_(1, lcnumelleft), _),
                                0., tmp1(_(1, rcnumelleft), _));
                flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans, 1.,
                                Bref(_(rcnumelleft+1,
                                       rcnumelleft+rcnumelright),
                                     _(lcnumelleft+1,
                                       lcnumelleft+lcnumelright)),
                                leftMat(_(lcnumelleft+1,
                                          lcnumelleft+lcnumelright), _),
                                0., tmp1(_(rcnumelleft+1,
                                           rcnumelleft+rcnumelright), _));

                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans, 1.,
                                rightMat(_(1, rcnumelleft), _),
                                tmp1(_(1, rcnumelleft), _),
                                0.,tmp2);
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans, 1.,
                                rightMat(_(rcnumelleft+1,
                                           rcnumelleft+rcnumelright), _),
                                tmp1(_(rcnumelleft+1,
                                       rcnumelleft+rcnumelright), _),
                                1., tmp2);

                //einsortiert werden muss hier nicht, da wir nicht mit QR weiterverarbeiten ist dies schon die Blockstruktur!

                node->getContent()->setUorB(tmp2);
                node->getContent()->setNumRows(1);
                node->getContent()->setLeftChildNumRows(tmp2.numCols());
                node->getContent()->setRightChildNumRows(tmp2.numRows());

            } else {
                int numel        = node->getContent()->getNumRows();
                Matrix &leftMat  = TITR.getNode()->getfirstChild()
                                   ->getContent()->getUorB();
                Matrix &rightMat = TITR.getNode()->getlastChild()
                                   ->getContent()->getUorB();

                Matrix tmp1, tmp2, tmp3, tmp4;
                Matrix& QRmat = QRmats[counter];
                Matrix& U     = Us[counter];
                Matrix& Q     = Qs[counter];
                int rank      = ranks[counter];
                int m         = leftMat.numCols()*rightMat.numCols();
                int n         = numel;
                Vec ks        = scals[counter];

                for (int i=1; i<=rank; ++i) {
                    flens::blas::scal(ks(i), U(_, i));
                }


                if ( m < n){
                    tmp2 = U(_,_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank)));
                    tmp3 = tmp2;
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    QRmat(_,_(1,m)).upper(),tmp2);
                    flens::blas::mm(cxxblas::Trans,
                                    cxxblas::NoTrans, 1.0,
                                    QRmat(_,_(m+1,n)),tmp3,0.0,tmp4);
                    tmp1 = Matrix(tmp2.numRows()+tmp4.numRows(),tmp2.numCols());
                    tmp1(_(1,tmp2.numRows()),_) = tmp2;
                    tmp1(_(tmp2.numRows()+1,tmp1.numRows()),_) = tmp4;
                } else {
                    tmp1 = U(_(1,n),_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank)));
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    QRmat(_(1,n),_).upper(),tmp1);

                }

                TITR.getNode()->getContent()->setUorB(tmp1);

                //Das hier ist Q_t*S_t, die muss nun zurück ins Blockformat!!!
                tmp1.resize(0, 0);
                flens::blas::mm(cxxblas::NoTrans,
                                cxxblas::NoTrans, 1.0,
                                Q,
                                U(_,_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank))),0.0,tmp1);
                int newnumel   = tmp1.numCols();
                int newlcnumel = node->getfirstChild()->getContent()->getNumRows();
                int newrcnumel = node->getlastChild()->getContent()->getNumRows();
                Matrix newB(newnumel*newrcnumel,newlcnumel);

                for(int i = 1; i <= newnumel; ++i){
                    for(int j = 1; j <= newlcnumel; ++j){
                        newB(_((i-1)*newrcnumel+1,i*newrcnumel),j) = tmp1(_((j-1)*newrcnumel + 1,j*newrcnumel) ,i);
                    }
                }

                node->getContent()->setUorB(newB);
                node->getContent()->setNumRows(newnumel);
                node->getContent()->setLeftChildNumRows(newlcnumel);
                node->getContent()->setRightChildNumRows(newrcnumel);
            }

            TIT--;
            TITR--;
            TITgram--;
            TITleft--;
            TITright--;
            ++counter;
            newlevel = true;
        }
    }

    return tree;
}


template <typename T>
HTuckerTree<T>  operator+(const HTuckerTree<T>  & tree1, const HTuckerTree<T> & tree2){
    using flens::_;

	assert(tree1.getGeneralTree().root->getContent()->getIndex().length() == tree2.getGeneralTree().root->getContent()->getIndex().length());
	
	HTuckerTree<T>  tmp = tree1;
	GeneralTreeNode<HTuckerTreeNode<T> > *node1,*node2,*nodetmp;

    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT1 = tree1.getGeneralTree().end();
	GeneralTreeIterator<HTuckerTreeNode<T> > TIT2 = tree2.getGeneralTree().end();
	GeneralTreeIterator<HTuckerTreeNode<T> > TITtmp = tmp.getGeneralTree().end();
    for(;TIT1 >= tree1.getGeneralTree().begin(); TIT1--,TIT2--,TITtmp--){
		
		node1 = TIT1.getNode();
		node2 = TIT2.getNode();
		nodetmp = TITtmp.getNode();
		if(TIT1.getNode()->isLeaf()){
			flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > &node1UorB = node1->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > &node2UorB =  node2->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > tmpM(std::max(node1UorB.numRows(),node2UorB.numRows()),node1UorB.numCols()+node2UorB.numCols());
			tmpM(_(1,node1->getContent()->getUorB().numRows()),_(1,node1UorB.numCols())) = node1UorB;
			tmpM(_(1,node2->getContent()->getUorB().numRows()),_(node1UorB.numCols()+1,node1UorB.numCols()+ node2UorB.numCols())) = node2UorB;
			nodetmp->getContent()->setUorB(tmpM);
			nodetmp->getContent()->setNumRows(node1UorB.numCols()+node2UorB.numCols());
		} else if(TIT1.getNode()->isRoot()){
			//tree1.print_w_UorB();
			//tree2.print_w_UorB();
			int rcnumel1 = node1->getContent()->getRightChildNumRows();
			int rcnumel2 = node2->getContent()->getRightChildNumRows();
			int lcnumel1 = node1->getContent()->getLeftChildNumRows();
			int lcnumel2 = node2->getContent()->getLeftChildNumRows();
			// the variables numel1 and numel2 should be = 1;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > UorBsum(rcnumel1+rcnumel2,lcnumel1+lcnumel2);
			UorBsum(_(1,rcnumel1),_(1,lcnumel1)) = node1->getContent()->getUorB();
			UorBsum(_(rcnumel1+1, rcnumel1+rcnumel2),_(lcnumel1+1,lcnumel1+lcnumel2)) = node2->getContent()->getUorB();
			nodetmp->getContent()->setUorB(UorBsum);
			nodetmp->getContent()->setNumRows(1);
			nodetmp->getContent()->setLeftChildNumRows(lcnumel1 + lcnumel2);
			nodetmp->getContent()->setRightChildNumRows(rcnumel1 + rcnumel2);
			
		} else {
			//std::cout << "else-Teil" << std::endl;
			int rcnumel1 = node1->getContent()->getRightChildNumRows();
			int rcnumel2 = node2->getContent()->getRightChildNumRows();
			int lcnumel1 = node1->getContent()->getLeftChildNumRows();
			int lcnumel2 = node2->getContent()->getLeftChildNumRows();
			int numel1 = node1->getContent()->getNumRows();
			int numel2 = node2->getContent()->getNumRows();
			
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > UorBsum((rcnumel1+rcnumel2)*(numel1+numel2),lcnumel1 + lcnumel2);
			for(int i = 1; i <= numel1; ++i){
				UorBsum(_((i-1)*(rcnumel1 + rcnumel2) + 1,(i-1)*(rcnumel1 + rcnumel2) + rcnumel1),_(1,lcnumel1)) = (node1->getContent()->getUorB())(_((i-1)*rcnumel1+1,i*rcnumel1),_);
			}
			for(int i = 1; i<= numel2; ++i){
				UorBsum(_(numel1*(rcnumel1 + rcnumel2) + (i-1) * (rcnumel1 + rcnumel2) + rcnumel1 + 1,numel1*(rcnumel1 + rcnumel2) + i * (rcnumel1 + rcnumel2)),_(lcnumel1+1,lcnumel1 + lcnumel2)) = (node2->getContent()->getUorB())(_((i-1)*rcnumel2+1,i*rcnumel2),_);
			}
			nodetmp->getContent()->setUorB(UorBsum);
			nodetmp->getContent()->setNumRows(numel1 + numel2);
			nodetmp->getContent()->setLeftChildNumRows(lcnumel1 + lcnumel2);
			nodetmp->getContent()->setRightChildNumRows(rcnumel1 + rcnumel2);
		}

	}
    }

	return tmp;
}


template <typename T>
HTuckerTree<T> operator*(const double number, const HTuckerTree<T>  & tree){
	HTuckerTree<T>  tmp = tree;
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmpmat;
	
	tmpmat = tmp.getGeneralTree().root->getContent()->getUorB();
	//blas::scal(number,tmpmat);
	tmpmat *= number;
	//tmpmat = number * tmp.getGeneralTree().root->getContent()->getUorB();
	tmp.getGeneralTree().root->getContent()->setUorB(tmpmat);
	return tmp;
}

template <typename T>
HTuckerTree<T> operator*(const HTuckerTree<T> & tree, const double number){
	return number * tree;
}

template <typename T>
HTuckerTree<T> operator-(const HTuckerTree<T>  & tree1, const HTuckerTree<T>  & tree2){
	return tree1 + (-1)*tree2;
}

template <typename T>
HTuckerTree<T> operator*(const HTuckerTree<T> & tree1, const HTuckerTree<T> & tree2){
    using flens::_;

	assert(tree1.dim() == tree2.dim());

	HTuckerTree<T>  tmp = tree1;
	GeneralTreeNode<HTuckerTreeNode<T> > *node1,*node2,*nodetmp;
	
    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT1 = tree1.getGeneralTree().end();
	GeneralTreeIterator<HTuckerTreeNode<T> > TIT2 = tree2.getGeneralTree().end();
	GeneralTreeIterator<HTuckerTreeNode<T> > TITtmp = tmp.getGeneralTree().end();
    for(;TIT1 >= tree1.getGeneralTree().begin(); TIT1--,TIT2--,TITtmp--){
		node1 = TIT1.getNode();
		node2 = TIT2.getNode();
		nodetmp = TITtmp.getNode();
		if(TIT1.getNode()->isLeaf()){
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > & node1UorB = node1->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > & node2UorB = node2->getContent()->getUorB();
			int len = node1UorB.numRows() / node2UorB.numRows();
			assert(node1UorB.numRows() == node2UorB.numRows() * len);
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Rneu(len,node1UorB.numCols() * node2UorB.numCols());
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > matversion(len,node2UorB.numRows());
			int node1nc = node1->getContent()->getUorB().numCols();
			int node2nr = node2->getContent()->getUorB().numRows();
			for(int i = 1; i <= node1nc; ++i){
				for(int j = 1; j <= node2nr; ++j){
					matversion(_,j) = node1UorB(_((j-1)*len + 1, j * len),i);
				}
				Rneu(_,_((i-1)*node2UorB.numCols()+1,i * node2UorB.numCols())) = matversion * node2UorB;
			}
			nodetmp->getContent()->setUorB(Rneu);
			nodetmp->getContent()->setNumRows(node1UorB.numCols() * node2UorB.numCols());
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
}

template <typename T>
HTuckerTree<T> ones(const DimensionIndex &dims){
	DenseVectorList<T> list;
	for(int i = 0; i < dims.length(); ++i){
		flens::DenseVector<flens::Array<T> > dv(dims[i]);
		for(int j = dv.firstIndex(); j <= dv.lastIndex(); ++j){
			dv(j) = 1;
		}
		list.add(dv);
	}

	DimensionIndex minis(1,dims.length());
	HTuckerTree<T> tree(dims.length());
	tree.generateTofElementary(list,1,dims.length());
	return tree;
}

template <typename T>
HTuckerTree<T> identity(const DimensionIndex &dims){
	DenseVectorList<T> list;
	DimensionIndex maxi(dims.length());
	for(int i = 0; i < dims.length(); ++i){
		maxi[i] = dims[i]*dims[i];
		flens::DenseVector<flens::Array<T> > dv(dims[i]*dims[i]);	
		for(int j = 1; j <= dims[i]; ++j){
			dv((j-1)*dims[i] + j-1 + dv.firstIndex()) = 1;
		}
		list.add(dv);
	}

	DimensionIndex minis(1,dims.length());
	HTuckerTree<T> tree(dims.length());
	tree.generateTofElementary(list,1,dims.length());
	return tree;
}

template <typename T> 
HTuckerTree<T>
transpose(const HTuckerTree<T> & tree){
    using flens::_;

	HTuckerTree<T> copy = tree;
	//tree.print_w_UorB();
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT1 = copy.getGeneralTree().end(); TIT1 >= copy.getGeneralTree().begin(); TIT1--){
		if(TIT1.getNode()->isLeaf()){
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > & Uref = TIT1.getNode()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > U(Uref.numRows(),Uref.numCols());
			int sqrtN = sqrt(Uref.numRows()+0.0); //Assuming the matrix is a square Matrix!
			for(int i = 1; i <= sqrtN; ++i){
				for(int j = 1; j <= sqrtN; ++j){
					U((i-1)*sqrtN + j,_) = Uref((j-1)*Uref + i,_);
				}
			}
			TIT1.getNode()->getContent()->setUorB(U);
		}
	
	}
	return copy;

}

template <typename T>
HTuckerTree<T> vec2diag(const HTuckerTree<T> &tree){
    using flens::_;

	HTuckerTree<T> copy = tree;
	//tree.print_w_UorB();
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT1 = copy.getGeneralTree().end(); TIT1 >= copy.getGeneralTree().begin(); TIT1--){
		if(TIT1.getNode()->isLeaf()){
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > & Uref = TIT1.getNode()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > U(Uref.numRows()*Uref.numRows(),Uref.numCols());
			for(int i = 1; i <= Uref.numRows(); ++i){
				U((i-1)*Uref.numRows() + i,_) = Uref(i,_);
			}
			TIT1.getNode()->getContent()->setUorB(U);
		}
	
	}
	return copy;
}



template <typename T>
HTuckerTree<T> concat(const HTuckerTree<T> &tree1, const HTuckerTree<T> &tree2){
	//std::cout << tree1.getGeneralTree().root->getContent()->getIndex() << std::endl;

	int maximum = 0;
	DimensionIndex & indref1 = tree1.getGeneralTree().root->getContent()->getIndex();
	DimensionIndex & indref2 = tree2.getGeneralTree().root->getContent()->getIndex();
	for(int i = 0; i < indref1.length(); ++i){
		if(indref1[i] > maximum){
			maximum = indref1[i];
		}
	}
	int dneu = indref1.length() + indref2.length();
	
	//std::cout << "tree1.minval = " <<  tree1.getGeneralTree().root->getContent()->getROA().A.minval << "  tree2.minval = " <<  tree2.getGeneralTree().root->getContent()->getROA().A.minval << std::endl;
	//std::cout << "tree1.maxval = " <<  tree1.getGeneralTree().root->getContent()->getROA().A.maxval << "  tree2.maxval = " <<  tree2.getGeneralTree().root->getContent()->getROA().A.maxval << std::endl;

		//std::cout << "maximum = " << maximum << std::endl;
	//std::cout << "dneu = " << dneu << std::endl;
	//std::cout << "minval = " << tree1.getGeneralTree().root->getContent()->getROA().A.minval << std::endl;
	//std::cout << "maxval = " << tree1.getGeneralTree().root->getContent()->getROA().A.maxval << std::endl;
	
	DimensionIndex dn(dneu);
	dn.setValueAsc();

	//std::cout << dn << std::endl;
	
	HTuckerTree<T> rettree;
	
	rettree.setdim(dneu);
	HTuckerTreeNode<T> rootnode(dn);

	rettree.getGeneralTree().root = new GeneralTreeNode<HTuckerTreeNode<T> >(rootnode);
	GeneralTreeNode<HTuckerTreeNode<T> > * node = rettree.getGeneralTree().root;
	GeneralTreeNode<HTuckerTreeNode<T> > * t1node = tree1.getGeneralTree().root;
	GeneralTreeNode<HTuckerTreeNode<T> > * t2node = tree2.getGeneralTree().root;
	
	node->appendChild(*(t1node->getContent() ));
	node->appendChild(*(t2node->getContent() ));
	DimensionIndex idx =  node->getfirstChild()->getnextSibling()->getContent()->getIndex();
	for(int i = 0; i < idx.length(); ++i){
		idx[i] = idx[i] + indref1.length();
	}
	node->getfirstChild()->getnextSibling()->getContent()->setIndex(idx);

	GeneralTreeNode<HTuckerTreeNode<T> > *levelup = NULL;
	GeneralTreeNode<HTuckerTreeNode<T> > *retlevelup = NULL;
	GeneralTreeNode<HTuckerTreeNode<T> > * tnode = node->getfirstChild();
	GeneralTreeNode<HTuckerTreeNode<T> > * child;
	
	while(true){
		child = t1node->getfirstChild();
		while(child != NULL){
			//std::cout << "wird angehängt: " << child->getContent()->getIndex() << std::endl;
			tnode->appendChild(*(child->getContent()));
			//std::cout << "wird überschrieben: " << tnode->getfirstChild()->getpreviousSibling()->getContent()->getIndex() << std::endl;
			child = child->getnextSibling();
			if(child == child->getParent()->getfirstChild()){
				break;
			}
		}
		
	
		if(levelup == NULL && t1node->getfirstChild() != NULL){
			levelup = t1node->getfirstChild();
			retlevelup = tnode->getfirstChild();
		}

		if(t1node->getlevelright() != NULL){
			t1node = t1node->getlevelright();
			tnode = tnode->getlevelright();
		} else if(levelup != NULL){
			t1node = levelup;
			levelup = NULL;
			tnode = retlevelup;
		} else {
			break;
		}
	}



	tnode = node->getfirstChild()->getnextSibling();
	
	//tnode->getContent()->getROA().A.minval = empty.minval;
	//tnode->getContent()->getROA().A.maxval = empty.maxval;
	t1node = t2node;
	HTuckerTreeNode<T> save;

	while(true){
		child = t1node->getfirstChild();
		while(child != NULL){
			tnode->appendChild(*(child->getContent()));
			idx = tnode->getfirstChild()->getpreviousSibling()->getContent()->getIndex();
			for(int i = 0; i < idx.length(); ++i){
				idx[i] = idx[i] + indref1.length();
			}
			tnode->getfirstChild()->getpreviousSibling()->getContent()->setIndex(idx);
		
			child = child->getnextSibling();
			if(child == child->getParent()->getfirstChild()){
				break;
			}
		}
		
	
		if(levelup == NULL && t1node->getfirstChild() != NULL){
			levelup = t1node->getfirstChild();
			retlevelup = tnode->getfirstChild();
		}

		if(t1node->getlevelright() != NULL){
			t1node = t1node->getlevelright();
			tnode = tnode->getlevelright();
		} else if(levelup != NULL){
			t1node = levelup;
			levelup = NULL;
			tnode = retlevelup;
		} else {
			break;
		}
	}


	node->getContent()->setRightChildNumRows(node->getlastChild()->getContent()->getNumRows());
	node->getContent()->setLeftChildNumRows(node->getfirstChild()->getContent()->getNumRows());
	node->getContent()->setNumRows(1);

	flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > g(1,1);
	g = 1;
	node->getContent()->setUorB(g);

	return rettree;
}




template <typename T>
HTuckerTree<T> subtree(const HTuckerTree<T> &tree, const int mindim, const int maxdim){
	HTuckerTree<T> copy = tree;
	GeneralTreeNode<HTuckerTreeNode<T> > * node;
	DimensionIndex index;
	bool found = false;
	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT1 = copy.getGeneralTree().end(); TIT1 >= copy.getGeneralTree().begin(); TIT1--){
		index = TIT1.getNode()->getContent()->getIndex();
		if(index[0] == mindim && index[index.length() - 1] == maxdim){
			node = TIT1.getNode();
			found = true;
			break;
		}
	}

	if(!found || node->isRoot()){
		std::cout << "return copy" << std::endl;
		return copy;
	}

	copy.setdim(maxdim-mindim+1);
	
	GeneralTreeNode<HTuckerTreeNode<T> > * htroot = copy.getGeneralTree().root;
	GeneralTreeNode<HTuckerTreeNode<T> > * hnode = node;
	
	
	while(hnode != htroot){
		if(hnode->getParent()->getfirstChild() == hnode){
		    hnode->getParent()->removeChild(2);
		} else {
			hnode->getParent()->removeChild(1);
		}
		hnode = hnode->getParent();
	}

	/*std::cout << " copy ----------------------" << std::endl;

	copy.print_w_UorB();

	std::cout << "copy ende" << std::endl << std::endl;*/
	
	node->getParent()->setfirstChild(0);
	node->getParent()->setlastChild(0);
	node->setParent(0);

	if(!(htroot->getfirstChild() == 0)){
		htroot->removeChild(1);
	}
	delete htroot;
	copy.getGeneralTree().root = node;


	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT1 = copy.getGeneralTree().end(); TIT1 >= copy.getGeneralTree().begin(); TIT1--){
		index = TIT1.getNode()->getContent()->getIndex();
		for(int i = 0; i < index.length(); ++i){
			index[i] -= mindim - 1;
		}
		TIT1.getNode()->getContent()->setIndex(index);
	}
	return copy;
		
}


template <typename T>
HTuckerTree<T> gramians_orthogonal3(const HTuckerTree<T>& tree)
{
    using flens::_;

    typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Matrix;

    GeneralTreeNode<HTuckerTreeNode<T> > *node;
    GeneralTreeNode<HTuckerTreeNode<T> > *gramnode;
    HTuckerTreeNode<T> *parent;
    HTuckerTree<T> gram;
    gram.set_tree(tree);

    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.getGeneralTree().begin();
    GeneralTreeIterator<HTuckerTreeNode<T> > GramTIT = gram.getGeneralTree().begin();

    std::chrono::duration<T> diff;
    auto tmm   = diff.count();
    tmm        = 0.;
    auto talloc = tmm;
    auto tcopy = tmm;
    auto tdata = tmm;

    for(; TIT <= tree.getGeneralTree().end(); TIT++,GramTIT++){
        node = TIT.getNode();
        gramnode = GramTIT.getNode();
        if(node->isRoot()){
            auto t0 = std::chrono::high_resolution_clock::now();

            Matrix M(1,1);
            M = 1;
            gramnode->getContent()->setUorB(M);

            auto t1 = std::chrono::high_resolution_clock::now();
            diff    = t1 - t0;
            talloc  += diff.count();
        } else {
            auto t0 = std::chrono::high_resolution_clock::now();

            parent = node->getParent()->getContent();

            int numel         = node->getContent()->getNumRows();
            int parentnumel   = parent->getNumRows();
            int parentrcnumel = parent->getRightChildNumRows();
            int parentlcnumel = parent->getLeftChildNumRows();
            Matrix &parentG   = gramnode->getParent()->getContent()->getUorB();
            Matrix &parentB   = parent->getUorB();
            Matrix &G         = gramnode->getContent()->getUorB();
            G.resize(numel, numel);

            auto t1 = std::chrono::high_resolution_clock::now();
            diff    = t1 - t0;
            talloc  += diff.count();
            if(node == node->getParent()->getlastChild()){
                t0 = std::chrono::high_resolution_clock::now();
                /* G_tr */
                Matrix B31(parentnumel*parentlcnumel, parentrcnumel);
                t1      = std::chrono::high_resolution_clock::now();
                diff    = t1 - t0;
                talloc  += diff.count();

                t0 = std::chrono::high_resolution_clock::now();
                for (int i=1; i<=parentnumel; ++i) {
                    for (int j=1; j<=parentlcnumel; ++j) {
                        B31((i-1)*parentlcnumel+j, _) =
                        parentB(_((i-1)*parentrcnumel+1, i*parentrcnumel), j);
                    }
                }
                t1      = std::chrono::high_resolution_clock::now();
                diff    = t1 - t0;
                tcopy  += diff.count();

                /* Compute tensor product */
                t0 = std::chrono::high_resolution_clock::now();
                Matrix B31block(parentlcnumel, parentnumel);
                Matrix tensor(parentnumel*parentlcnumel, parentrcnumel);
                Matrix tmp(parentlcnumel, parentnumel);
                t1      = std::chrono::high_resolution_clock::now();
                diff    = t1 - t0;
                talloc  += diff.count();

                for (int i=1; i<=parentrcnumel; ++i) {
                    t0 = std::chrono::high_resolution_clock::now();
                    for (int j=1; j<=parentnumel; ++j) {
                        B31block(_, j) =
                        B31(_((j-1)*parentlcnumel+1, j*parentlcnumel), i);
                    }
                    t1      = std::chrono::high_resolution_clock::now();
                    diff    = t1 - t0;
                    tcopy  += diff.count();

                    t0 = std::chrono::high_resolution_clock::now();
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                    1., B31block, parentG, 0., tmp);
                    t1      = std::chrono::high_resolution_clock::now();
                    diff    = t1 - t0;
                    tmm  += diff.count();

                    t0 = std::chrono::high_resolution_clock::now();
                    for (int j=1; j<=parentnumel; ++j) {
                        tensor(_((j-1)*parentlcnumel+1, j*parentlcnumel), i) =
                        tmp(_, j);
                    }
                    t1      = std::chrono::high_resolution_clock::now();
                    diff    = t1 - t0;
                    tcopy  += diff.count();
                }

                /* Compute G_tr */
                t0 = std::chrono::high_resolution_clock::now();
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., B31, tensor, 0., G);
                t1 = std::chrono::high_resolution_clock::now();
                diff    = t1 - t0;
                tmm  += diff.count();
            } else {
                /* G_tl */
                /* Compute tensor product */
                t0 = std::chrono::high_resolution_clock::now();
                Matrix B32block(parentrcnumel, parentnumel);
                Matrix tensor(parentnumel*parentrcnumel, parentlcnumel);
                Matrix tmp(parentrcnumel, parentnumel);
                t1 = std::chrono::high_resolution_clock::now();
                diff    = t1 - t0;
                talloc  += diff.count();
                for (int i=1; i<=parentlcnumel; ++i) {
                    t0 = std::chrono::high_resolution_clock::now();
                    for (int j=1; j<=parentnumel; ++j) {
                        B32block(_, j) =
                        parentB(_((j-1)*parentrcnumel+1, j*parentrcnumel), i);
                    }
                    t1 = std::chrono::high_resolution_clock::now();
                    diff    = t1 - t0;
                    tcopy  += diff.count();

                    t0 = std::chrono::high_resolution_clock::now();
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                    1., B32block, parentG, 0., tmp);
                    t1    = std::chrono::high_resolution_clock::now();
                    diff  = t1 - t0;
                    tmm  += diff.count();

                    t0 = std::chrono::high_resolution_clock::now();
                    for (int j=1; j<=parentnumel; ++j) {
                        tensor(_((j-1)*parentrcnumel+1, j*parentrcnumel), i) =
                        tmp(_, j);
                    }
                    t1    = std::chrono::high_resolution_clock::now();
                    diff  = t1 - t0;
                    tcopy += diff.count();
                }

                /* Compute G_tl */
                t0 = std::chrono::high_resolution_clock::now();
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., parentB, tensor, 0., G);
                t1    = std::chrono::high_resolution_clock::now();
                diff  = t1 - t0;
                tmm += diff.count();
            }
        }
    }

    tdata = tcopy+talloc;
    std::cout << tdata  << " secs for data\n";
    std::cout << talloc << " secs to allocate\n";
    std::cout << tcopy  << " secs to copy\n";
    std::cout << tmm    << " secs to perform MM\n";
    return gram;
}


template <typename T>
HTuckerTree<T> gramians_orthogonal2(const HTuckerTree<T>& tree)
{
    using flens::_;

    typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Matrix;

    GeneralTreeNode<HTuckerTreeNode<T> > *node;
    GeneralTreeNode<HTuckerTreeNode<T> > *gramnode;
    HTuckerTreeNode<T> *parent;
    HTuckerTree<T> gram;
    gram.set_tree(tree);

    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.getGeneralTree().begin();
    GeneralTreeIterator<HTuckerTreeNode<T> > GramTIT = gram.getGeneralTree().begin();

    std::chrono::duration<T> diff;
    auto tmm   = diff.count();
    tmm        = 0.;
    auto tdata = tmm;
    auto talloc = tmm;
    auto tcopy = tmm;

    for(; TIT <= tree.getGeneralTree().end(); TIT++,GramTIT++){
        node = TIT.getNode();
        gramnode = GramTIT.getNode();
        if(node->isRoot()){
            auto t0 = std::chrono::high_resolution_clock::now();

            Matrix M(1,1);
            M = 1;
            gramnode->getContent()->setUorB(M);

            auto t1 = std::chrono::high_resolution_clock::now();
            diff    = t1 - t0;
            talloc  += diff.count();
        } else {
            auto t0 = std::chrono::high_resolution_clock::now();

            parent = node->getParent()->getContent();

            int numel         = node->getContent()->getNumRows();
            int parentnumel   = parent->getNumRows();
            int parentrcnumel = parent->getRightChildNumRows();
            int parentlcnumel = parent->getLeftChildNumRows();
            Matrix &parentG   = gramnode->getParent()->getContent()->getUorB();
            Matrix &parentB   = parent->getUorB();

            Matrix G(numel,numel);
            Matrix B3(parentnumel, parentlcnumel*parentrcnumel);
            Matrix Bmod(parentnumel, parentlcnumel*parentrcnumel);

            auto t1 = std::chrono::high_resolution_clock::now();
            diff    = t1 - t0;
            talloc  += diff.count();

            t0 = std::chrono::high_resolution_clock::now();
            for (int it = 1; it <= parentnumel; ++it) {
                for (int itr = 1; itr <= parentrcnumel; ++itr) {
                    B3(it, _((itr-1)*parentlcnumel + 1, itr*parentlcnumel))
                    = parentB((it-1)*parentrcnumel + itr, _);
                }
            }
            t1 = std::chrono::high_resolution_clock::now();
            diff    = t1 - t0;
            tcopy  += diff.count();

            t0 = std::chrono::high_resolution_clock::now();
            flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                            1., parentG, B3, 0., Bmod);
            t1      = std::chrono::high_resolution_clock::now();
            diff    = t1 - t0;
            tmm  += diff.count();

            if(node == node->getParent()->getlastChild()){
                /* G_tr */
                t0 = std::chrono::high_resolution_clock::now();
                Matrix Bt13(parentnumel*parentlcnumel, parentrcnumel);
                t1      = std::chrono::high_resolution_clock::now();
                diff    = t1 - t0;
                talloc  += diff.count();

                t0 = std::chrono::high_resolution_clock::now();
                for (int i = 1; i <= parentnumel; ++i) {
                    for (int j = 1; j <= parentlcnumel; ++j) {
                        Bt13((i-1)*parentlcnumel + j, _) =
                        parentB(_((i-1)*parentrcnumel + 1, i*parentrcnumel), j);
                    }
                }
                t1      = std::chrono::high_resolution_clock::now();
                diff    = t1 - t0;
                tcopy  += diff.count();

                t0 = std::chrono::high_resolution_clock::now();
                Matrix Bmod13(parentnumel*parentlcnumel, parentrcnumel);
                t1      = std::chrono::high_resolution_clock::now();
                diff    = t1 - t0;
                talloc  += diff.count();

                t0 = std::chrono::high_resolution_clock::now();
                for (int i = 1; i <= parentnumel; ++i) {
                    for (int j = 1; j <= parentrcnumel; ++j) {
                        Bmod13(_((i-1)*parentlcnumel + 1, i*parentlcnumel), j)
                        = Bmod(i, _((j-1)*parentlcnumel + 1, j*parentlcnumel));
                    }
                }
                t1      = std::chrono::high_resolution_clock::now();
                diff    = t1 - t0;
                tcopy  += diff.count();

                /* Compute G_tr */
                t0 = std::chrono::high_resolution_clock::now();
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., Bt13, Bmod13, 0., G);
                t1      = std::chrono::high_resolution_clock::now();
                diff    = t1 - t0;
                tmm  += diff.count();
            } else {
                /* G_tl */
                /* Compute tensor product */
                t0 = std::chrono::high_resolution_clock::now();
                Matrix Bmod23(parentnumel*parentrcnumel, parentlcnumel);
                t1      = std::chrono::high_resolution_clock::now();
                diff    = t1 - t0;
                talloc  += diff.count();

                t0 = std::chrono::high_resolution_clock::now();
                for (int i = 1; i <= parentnumel; ++i) {
                    for (int j = 1; j <= parentrcnumel; ++j) {
                        Bmod23((i-1)*parentrcnumel+j, _)
                        = Bmod(i, _((j-1)*parentlcnumel + 1, j*parentlcnumel));
                    }
                }
                t1      = std::chrono::high_resolution_clock::now();
                diff    = t1 - t0;
                tcopy  += diff.count();

                /* Compute G_tl */
                t0 = std::chrono::high_resolution_clock::now();
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., parentB, Bmod23, 0., G);
                t1      = std::chrono::high_resolution_clock::now();
                diff    = t1 - t0;
                tmm  += diff.count();
            }
            t0 = std::chrono::high_resolution_clock::now();
            gramnode->getContent()->setUorB(G);
            t1 = std::chrono::high_resolution_clock::now();
            diff    = t1 - t0;
            tcopy  += diff.count();
            talloc  += diff.count();
        }
    }

    tdata = tcopy+talloc;
    std::cout << tdata  << " secs for data\n";
    std::cout << talloc << " secs to allocate\n";
    std::cout << tcopy  << " secs to copy\n";
    std::cout << tmm   << " secs to perform MM\n";
    return gram;
}


template <typename T>
HTuckerTree<T> gramians_orthogonal(const HTuckerTree<T>& tree)
{
    using flens::_;

    typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Matrix;
    typedef typename Matrix::View                                     MView;

    GeneralTreeNode<HTuckerTreeNode<T> > *node;
    GeneralTreeNode<HTuckerTreeNode<T> > *gramnode;
    HTuckerTreeNode<T> *parent;
    HTuckerTree<T> gram;
    gram.set_tree(tree);

    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.getGeneralTree().begin();
    GeneralTreeIterator<HTuckerTreeNode<T> > GramTIT = gram.getGeneralTree().begin();

    std::chrono::duration<T> diff;
    auto tmm   = diff.count();
    tmm        = 0.;
    auto talloc = tmm;
    auto tcopy = tmm;
    auto tdata = tmm;

    auto t0 = std::chrono::high_resolution_clock::now();
    auto r = tree.max_rank();
    Matrix Buffer(r*r, 4*r);
    auto t1 = std::chrono::high_resolution_clock::now();
    diff    = t1 - t0;
    talloc  += diff.count();

    for(; TIT <= tree.getGeneralTree().end(); TIT++,GramTIT++){
        node = TIT.getNode();
        gramnode = GramTIT.getNode();
        if(node->isRoot()){
            t0 = std::chrono::high_resolution_clock::now();

            Matrix M(1,1);
            M = 1;
            gramnode->getContent()->setUorB(M);

            t1 = std::chrono::high_resolution_clock::now();
            diff    = t1 - t0;
            talloc  += diff.count();
        } else {
            t0 = std::chrono::high_resolution_clock::now();

            parent = node->getParent()->getContent();

            int numel         = node->getContent()->getNumRows();
            int parentnumel   = parent->getNumRows();
            int parentrcnumel = parent->getRightChildNumRows();
            int parentlcnumel = parent->getLeftChildNumRows();
            Matrix &parentG   = gramnode->getParent()->getContent()->getUorB();
            Matrix &parentB   = parent->getUorB();
            Matrix &G         = gramnode->getContent()->getUorB();
            G.resize(numel, numel);

            t1 = std::chrono::high_resolution_clock::now();
            diff    = t1 - t0;
            talloc  += diff.count();
            if(node == node->getParent()->getlastChild()){
                /* G_tr */
                MView B31 = Buffer(_(1, parentnumel*parentlcnumel),
                                   _(1, parentrcnumel));

                t0 = std::chrono::high_resolution_clock::now();
                for (int i=1; i<=parentnumel; ++i) {
                    for (int j=1; j<=parentlcnumel; ++j) {
                        B31((i-1)*parentlcnumel+j, _) =
                        parentB(_((i-1)*parentrcnumel+1, i*parentrcnumel), j);
                    }
                }
                t1      = std::chrono::high_resolution_clock::now();
                diff    = t1 - t0;
                tcopy  += diff.count();

                /* Compute tensor product */
                MView B31block = Buffer(_(1, parentlcnumel),
                                        _(r+1, r+parentnumel));
                MView tensor   = Buffer(_(1, parentnumel*parentlcnumel),
                                        _(2*r+1, 2*r+parentrcnumel));
                MView tmp      = Buffer(_(1, parentlcnumel), 
                                        _(3*r+1, 3*r+parentnumel));

                for (int i=1; i<=parentrcnumel; ++i) {
                    t0 = std::chrono::high_resolution_clock::now();
                    for (int j=1; j<=parentnumel; ++j) {
                        B31block(_, j) =
                        B31(_((j-1)*parentlcnumel+1, j*parentlcnumel), i);
                    }
                    t1      = std::chrono::high_resolution_clock::now();
                    diff    = t1 - t0;
                    tcopy  += diff.count();

                    t0 = std::chrono::high_resolution_clock::now();
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                    1., B31block, parentG, 0., tmp);
                    t1      = std::chrono::high_resolution_clock::now();
                    diff    = t1 - t0;
                    tmm  += diff.count();

                    t0 = std::chrono::high_resolution_clock::now();
                    for (int j=1; j<=parentnumel; ++j) {
                        tensor(_((j-1)*parentlcnumel+1, j*parentlcnumel), i) =
                        tmp(_, j);
                    }
                    t1      = std::chrono::high_resolution_clock::now();
                    diff    = t1 - t0;
                    tcopy  += diff.count();
                }

                /* Compute G_tr */
                t0 = std::chrono::high_resolution_clock::now();
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., B31, tensor, 0., G);
                t1 = std::chrono::high_resolution_clock::now();
                diff    = t1 - t0;
                tmm  += diff.count();
            } else {
                /* G_tl */
                /* Compute tensor product */
                MView B32block = Buffer(_(1, parentrcnumel),
                                        _(1, parentnumel));
                MView tensor   = Buffer(_(1, parentnumel*parentrcnumel),
                                        _(r+1, r+parentlcnumel));
                MView tmp      = Buffer(_(1, parentrcnumel),
                                        _(2*r+1, 2*r+parentnumel));
                for (int i=1; i<=parentlcnumel; ++i) {
                    t0 = std::chrono::high_resolution_clock::now();
                    for (int j=1; j<=parentnumel; ++j) {
                        B32block(_, j) =
                        parentB(_((j-1)*parentrcnumel+1, j*parentrcnumel), i);
                    }
                    t1 = std::chrono::high_resolution_clock::now();
                    diff    = t1 - t0;
                    tcopy  += diff.count();

                    t0 = std::chrono::high_resolution_clock::now();
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                    1., B32block, parentG, 0., tmp);
                    t1    = std::chrono::high_resolution_clock::now();
                    diff  = t1 - t0;
                    tmm  += diff.count();

                    t0 = std::chrono::high_resolution_clock::now();
                    for (int j=1; j<=parentnumel; ++j) {
                        tensor(_((j-1)*parentrcnumel+1, j*parentrcnumel), i) =
                        tmp(_, j);
                    }
                    t1    = std::chrono::high_resolution_clock::now();
                    diff  = t1 - t0;
                    tcopy += diff.count();
                }

                /* Compute G_tl */
                t0 = std::chrono::high_resolution_clock::now();
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., parentB, tensor, 0., G);
                t1    = std::chrono::high_resolution_clock::now();
                diff  = t1 - t0;
                tmm += diff.count();
            }
        }
    }

//    tdata = tcopy+talloc;
//    std::cout << tdata  << " secs for data\n";
//    std::cout << talloc << " secs to allocate\n";
//    std::cout << tcopy  << " secs to copy\n";
//    std::cout << tmm    << " secs to perform MM\n";
    return gram;
}


template <typename T>
HTuckerTree<T> gramians_nonorthogonal(const HTuckerTree<T> & tree)
{
    using flens::_;

    typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Matrix;

    GeneralTreeNode<HTuckerTreeNode<T> > *node, *nodem;
    GeneralTreeNode<HTuckerTreeNode<T> > *gramnode;
    HTuckerTreeNode<T> *parent;
    HTuckerTree<T> gram;
    gram.set_tree(tree);
    HTuckerTree<T> Mtree;
    Mtree.set_tree(tree);

    /* Compute M's */
    auto TIT      = tree.getGeneralTree().end();
    auto TITm     = Mtree.getGeneralTree().end();
    for (; TIT >= tree.getGeneralTree().begin(); TIT--, TITm--) {
        node      = TIT.getNode();
        nodem     = TITm.getNode();
        if(node->isLeaf()) {
            Matrix tmp;
            Matrix& U = node->getContent()->getUorB();
            flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                            1., U, U, 0., tmp);
            nodem->getContent()->setUorB(tmp);
        } else if (node->isInner()) {
            Matrix& Mr = nodem->getlastChild()->getContent()->getUorB();
            Matrix& Ml = nodem->getfirstChild()->getContent()->getUorB();

            Matrix& Bt     = node->getContent()->getUorB();

            int numel = node->getContent()->getNumRows();
            int rcnumel = TIT.getNode()->getContent()->getRightChildNumRows();
            int lcnumel = TIT.getNode()->getContent()->getLeftChildNumRows();

            Matrix tensor(rcnumel*numel, lcnumel);
            Matrix tmp(rcnumel, lcnumel);
            /* Compute tensor product */
            for (int i=1; i<=numel; ++i) {
                flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                1.,
                                Bt(_((i-1)*rcnumel+1, i*rcnumel), _),
                                Ml,
                                0., tmp);
                flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                1.,
                                Mr,
                                tmp,
                                0., tensor(_((i-1)*rcnumel+1, i*rcnumel), _));
            }

            /* Compute B_t^H*tensor */
            Matrix Mt(numel, numel);
            Matrix tensornonblock(rcnumel*lcnumel, numel);
            Matrix Bnonblock(rcnumel*lcnumel, numel);

            for (int i=1; i<=numel; ++i) {
                for (int j=1; j<=rcnumel; ++j) {
                    tensornonblock(_((j-1)*lcnumel+1, j*lcnumel), i) =
                    tensor((i-1)*rcnumel+j, _);

                    Bnonblock(_((j-1)*lcnumel+1, j*lcnumel), i) =
                    Bt((i-1)*rcnumel+j, _);
                }
            }

            flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                            1., Bnonblock,
                            tensornonblock,
                            0., Mt);

            nodem->getContent()->setUorB(Mt);
        }
    }

    /* Compute gram */
    TIT          = tree.getGeneralTree().begin();
    auto GramTIT = gram.getGeneralTree().begin();
    TITm         = Mtree.getGeneralTree().begin();
    for (; TIT <= tree.getGeneralTree().end(); TIT++, GramTIT++, TITm++) {
        node     = TIT.getNode();
        nodem    = TITm.getNode();
        gramnode = GramTIT.getNode();
        if(node->isRoot()){
            Matrix M(1,1);
            M = 1;
            gramnode->getContent()->setUorB(M);
        } else {
            parent = node->getParent()->getContent();

            auto numel         = node->getContent()->getNumRows();
            auto parentnumel   = parent->getNumRows();
            auto parentrcnumel = parent->getRightChildNumRows();
            auto parentlcnumel = parent->getLeftChildNumRows();
            Matrix &parentG   = gramnode->getParent()->getContent()->getUorB();
            Matrix &parentB   = parent->getUorB();

            Matrix& Mr = nodem->getParent()->
                         getlastChild()->getContent()->getUorB();
            Matrix& Ml = nodem->getParent()->
                         getfirstChild()->getContent()->getUorB();

            Matrix G(numel,numel);
            if(node == node->getParent()->getlastChild()){
                /* G_tr */
                Matrix B31(parentnumel*parentlcnumel, parentrcnumel);
                for (int i=1; i<=parentnumel; ++i) {
                    for (int j=1; j<=parentlcnumel; ++j) {
                        B31((i-1)*parentlcnumel+j, _) =
                        parentB(_((i-1)*parentrcnumel+1, i*parentrcnumel), j);
                    }
                }

                /* Compute tensor product */
                Matrix B31block(parentlcnumel, parentnumel);
                Matrix tensor(parentnumel*parentlcnumel, parentrcnumel);
                Matrix tmp2(parentlcnumel, parentnumel);
                Matrix tmp(parentlcnumel, parentnumel);
                for (int i=1; i<=parentrcnumel; ++i) {
                    for (int j=1; j<=parentnumel; ++j) {
                        B31block(_, j) =
                        B31(_((j-1)*parentlcnumel+1, j*parentlcnumel), i);
                    }
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                    1., B31block, parentG, 0., tmp2);
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                    1., Ml, tmp2, 0., tmp);
                    for (int j=1; j<=parentnumel; ++j) {
                        tensor(_((j-1)*parentlcnumel+1, j*parentlcnumel), i) =
                        tmp(_, j);
                    }
                }

                /* Compute G_tr */
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., B31, tensor, 0., G);
            } else {
                /* G_tl */
                /* Compute tensor product */
                Matrix B32block(parentrcnumel, parentnumel);
                Matrix tensor(parentnumel*parentrcnumel, parentlcnumel);
                Matrix tmp2(parentrcnumel, parentnumel);
                Matrix tmp(parentrcnumel, parentnumel);
                for (int i=1; i<=parentlcnumel; ++i) {
                    for (int j=1; j<=parentnumel; ++j) {
                        B32block(_, j) =
                        parentB(_((j-1)*parentrcnumel+1, j*parentrcnumel), i);
                    }
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                    1., B32block, parentG, 0., tmp2);
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                    1., Mr, tmp2, 0., tmp);
                    for (int j=1; j<=parentnumel; ++j) {
                        tensor(_((j-1)*parentrcnumel+1, j*parentrcnumel), i) =
                        tmp(_, j);
                    }
                }

                /* Compute G_tl */
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., parentB, tensor, 0., G);
            }
            gramnode->getContent()->setUorB(G);
        }
    }

    return gram;
}


template <typename T>
HTuckerTree<T>
contract(const HTuckerTree<T> &tree1,
         const HTuckerTree<T> &tree2)
{
    using flens::_;

    typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Matrix;

    GeneralTreeNode<HTuckerTreeNode<T> > *node1, *node2, *nodem;
    GeneralTreeNode<HTuckerTreeNode<T> > *gramnode;
    HTuckerTreeNode<T> *parent1, *parent2;
    HTuckerTree<T> gram;
    gram.set_tree(tree1);
    HTuckerTree<T> Mtree;
    Mtree.set_tree(tree1);

    /* Compute M's */
    auto TIT1     = tree1.getGeneralTree().end();
    auto TIT2     = tree2.getGeneralTree().end();
    auto TITm     = Mtree.getGeneralTree().end();
    for (; TIT1 >= tree1.getGeneralTree().begin(); TIT1--, TIT2--, TITm--) {
        node1      = TIT1.getNode();
        node2      = TIT2.getNode();
        nodem     = TITm.getNode();
        if(node1->isLeaf()) {
            Matrix tmp;
            Matrix& U1 = node1->getContent()->getUorB();
            Matrix& U2 = node2->getContent()->getUorB();
            flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                            1., U2, U1, 0., tmp);
            nodem->getContent()->setUorB(tmp);
        } else if (node1->isInner()) {
            Matrix& Mr = nodem->getlastChild()->getContent()->getUorB();
            Matrix& Ml = nodem->getfirstChild()->getContent()->getUorB();

            Matrix& Bt1 = node1->getContent()->getUorB();
            Matrix& Bt2 = node2->getContent()->getUorB();

            int numel1   = node1->getContent()->getNumRows();
            int numel2   = node2->getContent()->getNumRows();
            int rcnumel1 = TIT1.getNode()->getContent()->getRightChildNumRows();
            int rcnumel2 = TIT2.getNode()->getContent()->getRightChildNumRows();
            int lcnumel2 = TIT2.getNode()->getContent()->getLeftChildNumRows();

            Matrix tensor(rcnumel2*numel1, lcnumel2);
            Matrix tmp(rcnumel1, lcnumel2);
            /* Compute tensor product */
            for (int i=1; i<=numel1; ++i) {
                flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                1.,
                                Bt1(_((i-1)*rcnumel1+1, i*rcnumel1), _),
                                Ml,
                                0., tmp);
                flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                1.,
                                Mr,
                                tmp,
                                0., tensor(_((i-1)*rcnumel2+1, i*rcnumel2), _));
            }

            /* Compute B_t^H*tensor */
            Matrix Mt(numel2, numel1);
            Matrix tensornonblock(rcnumel2*lcnumel2, numel1);
            Matrix Bnonblock(rcnumel2*lcnumel2, numel2);

            for (int i=1; i<=numel1; ++i) {
                for (int j=1; j<=rcnumel2; ++j) {
                    tensornonblock(_((j-1)*lcnumel2+1, j*lcnumel2), i) =
                    tensor((i-1)*rcnumel2+j, _);
                }
            }

            for (int i=1; i<=numel2; ++i) {
                for (int j=1; j<=rcnumel2; ++j) {
                    Bnonblock(_((j-1)*lcnumel2+1, j*lcnumel2), i) =
                    Bt2((i-1)*rcnumel2+j, _);
                }
            }

            flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                            1., Bnonblock,
                            tensornonblock,
                            0., Mt);

            nodem->getContent()->setUorB(Mt);
        }
    }

    /* Compute gram */
    TIT1          = tree1.getGeneralTree().begin();
    TIT2          = tree2.getGeneralTree().begin();
    auto GramTIT  = gram.getGeneralTree().begin();
    TITm          = Mtree.getGeneralTree().begin();
    for (; TIT1 <= tree1.getGeneralTree().end(); TIT1++, TIT2++, GramTIT++, TITm++) {
        node1     = TIT1.getNode();
        node2     = TIT2.getNode();
        nodem     = TITm.getNode();
        gramnode  = GramTIT.getNode();
        if(node1->isRoot()){
            Matrix M(1,1);
            M = 1;
            gramnode->getContent()->setUorB(M);
        } else {
            parent1 = node1->getParent()->getContent();
            parent2 = node2->getParent()->getContent();

            auto numel1         = node1->getContent()->getNumRows();
            auto numel2         = node2->getContent()->getNumRows();
            auto parentnumel1   = parent1->getNumRows();
            auto parentnumel2   = parent2->getNumRows();
            auto parentrcnumel1 = parent1->getRightChildNumRows();
            auto parentrcnumel2 = parent2->getRightChildNumRows();
            auto parentlcnumel1 = parent1->getLeftChildNumRows();
            auto parentlcnumel2 = parent2->getLeftChildNumRows();
            Matrix &parentG     = gramnode->getParent()->getContent()->getUorB();
            Matrix &parentB1    = parent1->getUorB();
            Matrix &parentB2    = parent2->getUorB();

            Matrix& Mr = nodem->getParent()->
                         getlastChild()->getContent()->getUorB();
            Matrix& Ml = nodem->getParent()->
                         getfirstChild()->getContent()->getUorB();

            Matrix G(numel2, numel1);
            if(node1 == node1->getParent()->getlastChild()){
                /* G_tr */
                Matrix B31(parentnumel1*parentlcnumel1, parentrcnumel1);
                Matrix B31_2(parentnumel2*parentlcnumel2, parentrcnumel2);
                for (int i=1; i<=parentnumel1; ++i) {
                    for (int j=1; j<=parentlcnumel1; ++j) {
                        B31((i-1)*parentlcnumel1+j, _) =
                        parentB1(_((i-1)*parentrcnumel1+1, i*parentrcnumel1), j);
                    }
                }

                for (int i=1; i<=parentnumel2; ++i) {
                    for (int j=1; j<=parentlcnumel2; ++j) {
                        B31_2((i-1)*parentlcnumel2+j, _) =
                        parentB2(_((i-1)*parentrcnumel2+1, i*parentrcnumel2), j);
                    }
                }

                /* Compute tensor product */
                Matrix B31block(parentlcnumel1, parentnumel1);
                Matrix tensor(parentnumel2*parentlcnumel2, parentrcnumel1);
                Matrix tmp2(parentlcnumel1, parentnumel2);
                Matrix tmp(parentlcnumel2, parentnumel2);
                for (int i=1; i<=parentrcnumel1; ++i) {
                    for (int j=1; j<=parentnumel1; ++j) {
                        B31block(_, j) =
                        B31(_((j-1)*parentlcnumel1+1, j*parentlcnumel1), i);
                    }
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                    1., B31block, parentG, 0., tmp2);
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                    1., Ml, tmp2, 0., tmp);
                    for (int j=1; j<=parentnumel2; ++j) {
                        tensor(_((j-1)*parentlcnumel2+1, j*parentlcnumel2), i) =
                        tmp(_, j);
                    }
                }

                /* Compute G_tr */
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., B31_2, tensor, 0., G);
            } else {
                /* G_tl */
                /* Compute tensor product */
                Matrix B32block(parentrcnumel1, parentnumel1);
                Matrix tensor(parentnumel2*parentrcnumel2, parentlcnumel1);
                Matrix tmp2(parentrcnumel1, parentnumel2);
                Matrix tmp(parentrcnumel2, parentnumel2);
                for (int i=1; i<=parentlcnumel1; ++i) {
                    for (int j=1; j<=parentnumel1; ++j) {
                        B32block(_, j) =
                        parentB1(_((j-1)*parentrcnumel1+1, j*parentrcnumel1), i);
                    }
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans,
                                    1., B32block, parentG, 0., tmp2);
                    flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans,
                                    1., Mr, tmp2, 0., tmp);
                    for (int j=1; j<=parentnumel2; ++j) {
                        tensor(_((j-1)*parentrcnumel2+1, j*parentrcnumel2), i) =
                        tmp(_, j);
                    }
                }

                /* Compute G_tl */
                flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,
                                1., parentB2, tensor, 0., G);
            }
            gramnode->getContent()->setUorB(G);
        }
    }

    return gram;
}


template <typename T>
flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> >
contract_leaf(const HTuckerTree<T> &tree,
              const HTuckerTree<T> &gram,
              const int j)
{
    assert(j>0 && j<=tree.dim());
    htucker::DimensionIndex idx(1);
    idx[0] = j;

    typedef typename flens::GeMatrix
                    <flens::FullStorage
                    <T, flens::ColMajor> >      Matrix;
    Matrix ret;

    for (auto itt=tree.getGeneralTree().end(),
              itg=gram.getGeneralTree().end();
              itt>=tree.getGeneralTree().begin();
              itt--, itg--)
    {
        if (itt.getNode()->getContent()->getIndex()==idx) {
            const Matrix &U  = itt.getNode()->getContent()->getUorB();
            const Matrix &Gt = itg.getNode()->getContent()->getUorB();
            flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans, 1.,
                            U, Gt, 0., ret);
            return ret;
        }
    }

    std::cerr << "Error: idx not found\n";
    return ret;
}


template <typename T>
flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> >
contract(const HTuckerTree<T> &tree1,
         const HTuckerTree<T> &tree2,
         const int j)
{
    assert(j>0 && j<=tree1.dim());
    htucker::DimensionIndex idx(1);
    idx[0] = j;

    typedef typename flens::GeMatrix
                    <flens::FullStorage
                    <T, flens::ColMajor> >      Matrix;
    Matrix ret;
    HTuckerTree<T> gram = contract(tree1, tree2);

    for (auto itt=tree1.getGeneralTree().end(),
              itg=gram.getGeneralTree().end();
              itt>=tree1.getGeneralTree().begin();
              itt--, itg--)
    {
        if (itt.getNode()->getContent()->getIndex()==idx) {
            const Matrix &U  = itt.getNode()->getContent()->getUorB();
            const Matrix &Gt = itg.getNode()->getContent()->getUorB();
            flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1.,
                            U, Gt, 0., ret);
            return ret;
        }
    }

    std::cerr << "Error: idx not found\n";
    return ret;
}


template <typename T>
flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> >
projection(const HTuckerTree<T> &tree1,
           const HTuckerTree<T> &tree2,
           const int j)
{
    assert(j>0 && j<=tree1.dim());
    htucker::DimensionIndex idx(1);
    idx[0] = j;

    typedef typename flens::GeMatrix
                    <flens::FullStorage
                    <T, flens::ColMajor> >      Matrix;
    Matrix ret;
    HTuckerTree<T> gram = contract(tree1, tree2);

    for (auto itg=gram.getGeneralTree().end();
              itg>=gram.getGeneralTree().begin();
              itg--)
    {
        if (itg.getNode()->getContent()->getIndex()==idx) {
            const Matrix &Gt = itg.getNode()->getContent()->getUorB();
            return Gt;
        }
    }

    std::cerr << "Error: idx not found\n";
    return ret;
}


template <typename T>
const flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> >&
extract(const HTuckerTree<T> &tree,
        const DimensionIndex idx)
{
    for (auto itt=tree.getGeneralTree().end();
              itt>=tree.getGeneralTree().begin();
              itt--)
    {
        if (itt.getNode()->getContent()->getIndex()==idx) {
            return itt.getNode()->getContent()->getUorB();
        }
    }

    std::cerr << "Error: idx not found\n";
    exit(1); // High end exception handling.
}


template <typename T, typename I>
void
extract_core(const HTuckerTree<T>& tree,
                   std::vector<flens::GeMatrix<
                               flens::FullStorage<
                               T, flens::ColMajor> > >& C,
                   flens::DenseVector<
                   flens::Array<I> >&                   ranks)
{
    using flens::_;

    typedef typename flens::GeMatrix<
                     flens::FullStorage<
                     T, flens::ColMajor> >                      Matrix;
    typedef typename std::vector<Matrix>::size_type             size_type;
    typedef typename flens::DenseVector<
                     flens::Array<I> >                          Vector;

    size_type count=0;
    for (auto itt=tree.getGeneralTree().end();
              itt>=tree.getGeneralTree().begin();
              itt--)
    {
        auto node = itt.getNode();
        if (!node->isLeaf()) {
            auto &U   = node->getContent()->getUorB();
            auto rank = node->getContent()->getNumRows();

            if (C.size()>=count+1) C[count] = U;
            else                   C.push_back(U);

            if ((unsigned) ranks.length()>=count+1) {
                ranks(count+1) = rank;
            } else {
                Vector tmp = ranks;
                ranks.resize(count+1);
                ranks(_(1, tmp.length())) = tmp;
            }

            ++count;
        }
    }
}


template <typename T, typename I>
void
insert_core(      HTuckerTree<T>& tree,
            const std::vector<flens::GeMatrix<
                               flens::FullStorage<
                               T, flens::ColMajor> > >& C,
            const flens::DenseVector<
                   flens::Array<I> >&                   ranks)
{
    using flens::_;

    typedef typename flens::GeMatrix<
                     flens::FullStorage<
                     T, flens::ColMajor> >                      Matrix;
    typedef typename std::vector<Matrix>::size_type             size_type;

    size_type count=0;
    for (auto itt=tree.getGeneralTree().end();
              itt>=tree.getGeneralTree().begin();
              itt--)
    {
        auto node = itt.getNode();
        if (!node->isLeaf()) {
            auto     numel   = ranks(count+1);
            auto     lcnumel = C[count].numCols();
            unsigned rcnumel = C[count].numRows()/numel;
            node->getContent()->setUorB(C[count]);
            node->getContent()->setNumRows(numel);
            node->getContent()->setLeftChildNumRows(lcnumel);
            node->getContent()->setRightChildNumRows(rcnumel);

            ++count;
        }
    }
}


template <typename T>
void
insert(HTuckerTree<T> &tree,
       const flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> >& B,
       const DimensionIndex idx)
{
    for (auto itt=tree.getGeneralTree().end();
              itt>=tree.getGeneralTree().begin();
              itt--)
    {
        if (itt.getNode()->getContent()->getIndex()==idx) {
            itt.getNode()->getContent()->getUorB() = B;
            return;
        }
    }

    std::cerr << "Error: idx not found\n";
}


template <typename T>
HTuckerTree<T> gramians_nonorthogonal2(const HTuckerTree<T> & tree){
    using flens::_;

	HTuckerTree<T> Mtree(tree);

    {
	GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.getGeneralTree().end();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITm = Mtree.getGeneralTree().end();
    for(; TIT >= tree.getGeneralTree().begin(); TIT--,TITm--){
		if(TIT.getNode()->isLeaf()){
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Uref = TIT.getNode()->getContent()->getUorB();
			flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,Uref,Uref,0.0,tmp);
			TITm.getNode()->getContent()->setUorB(tmp);
		} else {
			// berechne M_tr*B_t*M_tl^T
			int rcnumel = TIT.getNode()->getContent()->getRightChildNumRows();
			int lcnumel = TIT.getNode()->getContent()->getLeftChildNumRows();
			int numel = TIT.getNode()->getContent()->getNumRows();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Mr = TITm.getNode()->getlastChild()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Ml = TITm.getNode()->getfirstChild()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Bt = TIT.getNode()->getContent()->getUorB();

			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp(numel*rcnumel,lcnumel);
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp2;
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp3;

			for(int i = 1; i<= numel; ++i){
				flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0,Bt(_((i-1)*rcnumel+1,i*rcnumel),_),Ml,0.0,tmp2);
				flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Mr,tmp2,0.0,tmp3);
				tmp(_((i-1)*rcnumel+1,i*rcnumel),_) = tmp3;
			}

			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Mnew(numel,numel);
			for(int i = 1; i <= numel; ++i){
				for(int j = 1; j <= numel; ++j){
					T sum = 0.0;
					for(int k = 1; k<= rcnumel; ++k){
						for(int l = 1; l<= lcnumel; ++l){
							sum += Bt((i-1)*rcnumel + k,l)*tmp((j-1)*rcnumel + k,l);
						}
					}
					Mnew(i,j) = sum;
				}
			}
			TITm.getNode()->getContent()->setUorB(Mnew);
		
		}
	}
    }

	GeneralTreeNode<HTuckerTreeNode<T> > *node;
	GeneralTreeNode<HTuckerTreeNode<T> > *gramnode;
	HTuckerTreeNode<T> *save,*last,*first,*parent;
	HTuckerTree<T> gram(tree);
	
    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = tree.getGeneralTree().begin();
    GeneralTreeIterator<HTuckerTreeNode<T> > GramTIT = gram.getGeneralTree().begin();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITm = Mtree.getGeneralTree().begin();
    for(; TIT <= tree.getGeneralTree().end(); TIT++,GramTIT++,TITm++){
		node = TIT.getNode();
		gramnode = GramTIT.getNode();
		save = node->getContent();
		last = node->getlastChild()->getContent();
		first = node->getfirstChild()->getContent();
		if(node->isRoot()){
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > M(1,1);
			M = 1;
			gramnode->getContent()->setUorB(M);
		} else { 
			parent = node->getParent()->getContent();
			
			int numel = node->getContent()->getNumRows();
			int parentnumel = parent->getNumRows();
			int lcnumel = save->getLeftChildNumRows();
			int rcnumel = save->getRightChildNumRows();
			int parentrcnumel = node->getParent()->getContent()->getRightChildNumRows();
			int parentlcnumel = node->getParent()->getContent()->getLeftChildNumRows();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &parentG = gramnode->getParent()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &parentB = parent->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > B(numel,numel);
			if(node == node->getParent()->getlastChild()){
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Ml= TITm.getNode()->getParent()->getfirstChild()->getContent()->getUorB();

				for(int i = 1; i <= numel; ++i){
					for(int j = 1; j <= numel; ++j){
						T sum = 0;
						for(int k = 1; k <= parentnumel; ++k){
							for(int l = 1; l <= parentnumel; ++l){
								for(int m1 = 1; m1 <= parentlcnumel; ++m1){
									for(int m2 = 1; m2 <= parentlcnumel; ++m2){
										sum += parentG(k,l)*parentB((k-1)*parentrcnumel + i, m1)*parentB((l-1)*parentrcnumel + j, m2)*Ml(m1,m2);
									}
								}
							}
						}
	
						B(i,j) = sum;
					}
				}
				
				////B31 im Blockformat abspeichern
				//flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > B31(parentnumel*parentrcnumel,parentlcnumel);
				//for(int i = 1; i <= parentnumel; ++i){
				//	for(int j = 1; j <= parentlcnumel; ++j){
				//		for(int k = 1; k <= parentrcnumel; ++k){
				//			B31((k-1)*parentnumel + i, j) = parentB((i-1)*parentrcnumel + k,j);
				//		}
				//	}
				//}

				////for each Block
				//flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > temp4(1,1),temp5(1,1);
				//for(int i = 1; i <=  parentrcnumel; ++i){
				//	flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0,B31(_((i-1)*parentnumel+1,i*parentnumel),_),Ml,0.0,temp4);
				//	flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,parentG,temp4,0.0,temp5);
				//	for(int j = 1; j <= parentrcnumel; ++j){
				//		T sum = 0;
				//		for(int k = 1; k <= parentnumel; ++k){
				//			for(int l = 1; l <= parentlcnumel; ++l){
				//				sum+= B31((j-1)*parentnumel + k, l)*temp5(k,l);
				//			}
				//		}
				//		B(i,j) = sum;
				//	}
				//}
				
				

			} else {
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Mr = TITm.getNode()->getParent()->getlastChild()->getContent()->getUorB();
				for(int i = 1; i <= numel; ++i){
					for(int j = 1; j <= numel; ++j){
						T sum = 0;
						for(int k = 1; k <= parentnumel; ++k){
							for(int l = 1; l <= parentnumel; ++l){
								for(int m1 = 1; m1 <= parentrcnumel; ++m1){
									for(int m2 = 1; m2 <= parentrcnumel; ++m2){
										sum += parentG(k,l)*parentB((k-1)*parentrcnumel + m2, i)*parentB((l-1)*parentrcnumel + m1, j)*Mr(m1,m2);
									}
								}
							}
						}
						B(i,j) = sum;
					}
				}
				////B32 im Blockformat abspeichern
				//flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > B32(parentnumel*parentrcnumel,parentlcnumel);
				//for(int i = 1; i <= parentnumel; ++i){
				//	for(int j = 1; j <= parentlcnumel; ++j){
				//		for(int k = 1; k <= parentrcnumel; ++k){
				//			B32((k-1)*parentnumel + i, j) = parentB((i-1)*parentrcnumel + k,j);
				//		}
				//	}
				//}

				//flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > temp4(1,1),temp5(1,1);
				//for(int i = 1; i <=  parentlcnumel; ++i){
				//	flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0,B32(_((i-1)*parentnumel+1,i*parentnumel),_),Mr,0.0,temp4);
				//	flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,parentG,temp4,0.0,temp5);
				//	for(int j = 1; j <= parentlcnumel; ++j){
				//		T sum = 0;
				//		for(int k = 1; k <= parentnumel; ++k){
				//			for(int l = 1; l <= parentrcnumel; ++l){
				//				sum+= B32((j-1)*parentnumel + k, l)*temp5(k,l);
				//			}
				//		}
				//		B(i,j) = sum;
				//	}
				//}
			
			}

			
			gramnode->getContent()->setUorB(B);
		}
	}
    }

	return gram;

}


template <typename T>
HTuckerTree<T> gramians_elementary(const HTuckerTree<T>& tree)
{
    /* Compute reduced gramians */
    typedef GeneralTreeIterator<htucker::HTuckerTreeNode<T>>          TreeIt;
    typedef GeneralTreeNode<htucker::HTuckerTreeNode<T>>              TreeNode;
    typedef HTuckerTreeNode<T>                                        HTNode;
    typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor>>  Matrix;

    HTuckerTree<T> Mtree;
    Mtree.set_tree(tree);


    TreeIt TIT  = tree.getGeneralTree().end();
    TreeIt TITm = Mtree.getGeneralTree().end();
    for (; TIT >= tree.getGeneralTree().begin(); TIT--,TITm--) {
        TreeNode *node  = TIT.getNode();
        TreeNode *nodem = TITm.getNode();
        if (node->isLeaf()) {
            Matrix  tmp;
            Matrix& U = TIT.getNode()->getContent()->getUorB();
            flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans,1.0,U, U,0.0,tmp);
            nodem->getContent()->setUorB(tmp);
        } else if (node->isInner()) {
            int numel   = node->getContent()->getNumRows();
            Matrix&  Mr = nodem->getlastChild()->getContent()->getUorB();
            Matrix&  Ml = nodem->getfirstChild()->getContent()->getUorB();
            Matrix   Mt(numel, numel);

            for (int i = 1; i <= numel; ++i) {
                int mapi = numel*(i-1)+i;
                for (int j = 1; j <= numel; ++j) {
                    int row = (mapi-1)%numel+1;
                    int col = (mapi-1)/numel+1;
                    Mt(i, j) = Ml(col, j)*Mr(row, j);
                }
            }
            nodem->getContent()->setUorB(Mt);
        }
    }

    htucker::HTuckerTree<T> gram;
    gram.set_tree(tree);

    TIT            = tree.getGeneralTree().begin();
    TreeIt GramTIT = gram.getGeneralTree().begin();
    TITm           = Mtree.getGeneralTree().begin();
    for(; TIT <= tree.getGeneralTree().end(); TIT++,GramTIT++,TITm++) {
        TreeNode *node, *gramnode, *nodem;
        HTNode *parent;

        node     = TIT.getNode();
        nodem    = TITm.getNode();
        gramnode = GramTIT.getNode();
        if (node->isRoot()) {
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > M(1,1);
            M = 1;
            gramnode->getContent()->setUorB(M);
        } else if (node->getParent()->isRoot()) {
            if (node == node->getParent()->getlastChild()) {
                Matrix& Ml = nodem->getParent()->
                             getfirstChild()->getContent()->getUorB();
                gramnode->getContent()->setUorB(Ml);
            } else {
                Matrix& Mr = nodem->getParent()->
                             getlastChild()->getContent()->getUorB();
                gramnode->getContent()->setUorB(Mr);
            }
        } else {
            parent          = node->getParent()->getContent();
            int numel       = node->getContent()->getNumRows();
            int parentnumel = parent->getNumRows();
            Matrix& parentG = gramnode->getParent()->
                                getContent()->getUorB();
            Matrix G(numel, numel);

            if (node == node->getParent()->getlastChild()) {
                Matrix& Ml = nodem->getParent()->
                             getfirstChild()->getContent()->getUorB();

                for (int i = 1; i <= numel; ++i) {
                    int mapi   = numel*(i-1)+i;
                    int blocki = (mapi-1)/parentnumel+1;
                    int rowG   = (mapi-1)%parentnumel+1;
                    for (int j = 1; j <= numel; ++j) {
                        int mapj   = numel*(j-1)+j;
                        int blockj = (mapj-1)/parentnumel+1;
                        int colG   = (mapj-1)%parentnumel+1;
                        G(i, j)    = Ml(blocki, blockj)*parentG(rowG, colG);
                    }
                }
            } else {
                Matrix& Mr = nodem->getParent()->
                             getlastChild()->getContent()->getUorB();
                for (int i = 1; i <= numel; ++i) {
                    int mapi   = numel*(i-1)+i;
                    int blocki = (mapi-1)/parentnumel+1;
                    int rowG   = (mapi-1)%parentnumel+1;
                    for (int j = 1; j <= numel; ++j) {
                        int mapj   = numel*(j-1)+j;
                        int blockj = (mapj-1)/parentnumel+1;
                        int colG   = (mapj-1)%parentnumel+1;
                        G(i, j)    = Mr(blocki, blockj)*parentG(rowG, colG);
                    }
                }
            }
            gramnode->getContent()->setUorB(G);
        }
    }

    return gram;
}


template <typename T>
void
HTuckerTree<T>::truncate(const int rank, bool isorth)
{
    using flens::_;

    HTuckerTree<T> Stree;
    Stree.set_tree(*this);
    if (!isorth) {
       this->orthogonalize();
    }

    HTuckerTree<T> gram = gramians_orthogonal(*this);
    flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > U(1,1);
    flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > Vt(1,1);
    flens::DenseVector<flens::Array<T> > s;
    GeneralTreeNode<HTuckerTreeNode<T> > *node;

    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TITgram = gram.getGeneralTree().begin();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITs= Stree.getGeneralTree().begin();
    for(; TITgram <= gram.getGeneralTree().end(); TITgram ++,TITs++){
		flens::svd(TITgram.getNode()->getContent()->getUorB(),s,U,Vt);
		TITs.getNode()->getContent()->setUorB(U(_,_(1,std::min(U.numRows(),
                                                    (FLENS_DEFAULT_INDEXTYPE)rank))));
	}
    }

    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = this->getGeneralTree().begin();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITs = Stree.getGeneralTree().begin();
	for(; TIT <= this->getGeneralTree().end(); TIT++,TITs++){
		node = TIT.getNode();
		
		if(node->isLeaf()){
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp;
			flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,node->getContent()->getUorB(),TITs.getNode()->getContent()->getUorB(),0.0,tmp);
			node->getContent()->setUorB(tmp);
			node->getContent()->setNumRows(tmp.numCols());
		} else {
			int rcnumel = node->getContent()->getRightChildNumRows();
			int lcnumel = node->getContent()->getLeftChildNumRows();
			
			int numel = node->getContent()->getNumRows();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Sref = TITs.getNode()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Slref = TITs.getNode()->getfirstChild()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Srref = TITs.getNode()->getlastChild()->getContent()->getUorB();
			int newnumel = Sref.numCols();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Bref = node->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp(rcnumel*newnumel,lcnumel);
			
			for(int l = 1; l <= newnumel; ++l){
				for(int j = 1; j <= rcnumel; ++j){
					for(int k = 1; k <= lcnumel; ++k){
						T sum = 0.0;
						for(int i = 1; i <= numel; ++i){
							sum += Bref((i-1)*rcnumel + j,k)*Sref(i,l);
						}
						tmp((l-1)*rcnumel + j,k) = sum;
					}
				}
			}
			//hier haben wir nun als B_t*S_t berechnet. Fehlt noch (S_tl^T x S_tr^T)*B_t*S_t
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > newB(newnumel*Srref.numCols(),Slref.numCols());
			for(int i = 1; i <= newnumel; ++i){
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp1;
				flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,tmp(_((i-1)*rcnumel+1,i*rcnumel),_),Slref,0.0,tmp1);
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp2;
				flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,Srref,tmp1,0.0,tmp2);
				newB(_((i-1)*Srref.numCols()+1,i*Srref.numCols()),_) = tmp2;
			}
			node->getContent()->setUorB(newB);
			node->getContent()->setNumRows(newnumel);
			node->getContent()->setLeftChildNumRows(Slref.numCols());
			node->getContent()->setRightChildNumRows(Srref.numCols());
		
		}
		
		
	}
    }
}


template <typename T>
void
HTuckerTree<T>::truncate(double eps, bool isorth)
{
    using flens::_;
    eps *= eps;

    HTuckerTree<T> Stree;
    Stree.set_tree(*this);
    if (!isorth) {
       this->orthogonalize();
    }

    HTuckerTree<T> gram = gramians_orthogonal(*this);

    flens::lapack::SVD::Job jobU  = flens::lapack::SVD::Save;
    flens::lapack::SVD::Job jobVt = flens::lapack::SVD::None;

    GeneralTreeNode<HTuckerTreeNode<T> > *node;

    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TITgram = gram.getGeneralTree().begin();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITs= Stree.getGeneralTree().begin();
    std::vector<flens::DenseVector<flens::Array<T> > >  svals;
    std::vector<flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> >>
                                                        Us;
    int count = 0;
    for(; TITgram <= gram.getGeneralTree().end(); TITgram ++){
        flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > U;
        flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > Vt;
        flens::DenseVector<flens::Array<T> > s;
        flens::DenseVector<flens::Array<T> > work;

        auto& Gt = TITgram.getNode()->getContent()->getUorB();
        auto  m  = Gt.numRows();
        auto  r  = std::min(m, Gt.numCols());
        U.resize(m, r);
        s.resize(r);

        flens::lapack::svd(jobU, jobVt,
                           Gt, s, U, Vt,
                           work);

        Us.push_back(U);
        svals.push_back(s);
        count += s.length();
    }

    flens::DenseVector<flens::Array<T> > s_all(count);
    count = 1;
    for (std::size_t i=0; i<svals.size(); ++i) {
        s_all(_(count,count+svals[i].length()-1)) = svals[i];
        count += svals[i].length();
    }


    flens::DenseVector<flens::Array<int> > rho(s_all.length());
    flens::sort(s_all, rho);
    double delta = 0.;
    count = 0;
    for (; count<s_all.length(); ++count) {
        delta += s_all(count+1);
        if (delta>eps) break;
    }

    flens::DenseVector<flens::Array<std::size_t> > dels(svals.size());
    for (int i=1; i<=count; ++i) {
        int right = 0;
        for (std::size_t j=0; j<svals.size(); ++j) {
            right += svals[j].length();
            if (rho(i)<=right) {
                ++dels(j+1);
                break;
            }
        }
    }

    for(std::size_t i=1; TITs <= Stree.getGeneralTree().end();
                         TITs++, ++i){
        int rank = svals[i-1].length() - dels(i);
        if (rank==0) rank = 1;
        TITs.getNode()->getContent()->
            setUorB(Us[i-1](_,_(1,std::min(Us[i-1].numRows(),
                               (FLENS_DEFAULT_INDEXTYPE)rank))));
    }
    }

    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = this->getGeneralTree().begin();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITs = Stree.getGeneralTree().begin();
    for(; TIT <= this->getGeneralTree().end(); TIT++,TITs++){
        node = TIT.getNode();

        if(node->isLeaf()){
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp;
            flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,node->getContent()->getUorB(),TITs.getNode()->getContent()->getUorB(),0.0,tmp);
            node->getContent()->setUorB(tmp);
            node->getContent()->setNumRows(tmp.numCols());
        } else {
            int rcnumel = node->getContent()->getRightChildNumRows();
            int lcnumel = node->getContent()->getLeftChildNumRows();

            int numel = node->getContent()->getNumRows();
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Sref = TITs.getNode()->getContent()->getUorB();
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Slref = TITs.getNode()->getfirstChild()->getContent()->getUorB();
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Srref = TITs.getNode()->getlastChild()->getContent()->getUorB();
            int newnumel = Sref.numCols();
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Bref = node->getContent()->getUorB();
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp(rcnumel*newnumel,lcnumel);
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Bt(rcnumel*lcnumel, numel);
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > BtSt(rcnumel*lcnumel, newnumel);

            /* Bt*St */
            for(int i = 1; i <= numel; ++i) {
                for(int j = 1; j <= lcnumel; ++j) {
                    Bt(_((j-1)*rcnumel+1, j*rcnumel), i) = Bref(_((i-1)*rcnumel+1, i*rcnumel), j);
                }
            }
            flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans, 1.0, Bt, Sref,0.0, BtSt);
            for(int i = 1; i <= newnumel; ++i) {
                for(int j = 1; j <= lcnumel; ++j) {
                    tmp(_((i-1)*rcnumel+1, i*rcnumel), j) = BtSt(_((j-1)*rcnumel+1, j*rcnumel), i);
                }
            }

            //hier haben wir nun als B_t*S_t berechnet. Fehlt noch (S_tl^T x S_tr^T)*B_t*S_t
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > newB(newnumel*Srref.numCols(),Slref.numCols());
            for(int i = 1; i <= newnumel; ++i){
                flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp1;
                flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,tmp(_((i-1)*rcnumel+1,i*rcnumel),_),Slref,0.0,tmp1);
                flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,Srref,tmp1,0.0,newB(_((i-1)*Srref.numCols()+1,i*Srref.numCols()),_));
            }
            node->getContent()->setUorB(newB);
            node->getContent()->setNumRows(newnumel);
            node->getContent()->setLeftChildNumRows(Slref.numCols());
            node->getContent()->setRightChildNumRows(Srref.numCols());
        }
    }
    }
}


template <typename T>
void
HTuckerTree<T>::truncate2(double eps, bool isorth)
{
    using flens::_;
    eps *= eps;
    eps /= 2.*d-3.;

    if (!isorth) {
       this->orthogonalize();
    }

    HTuckerTree<T> gram = gramians_orthogonal(*this);


    flens::lapack::SVD::Job jobU  = flens::lapack::SVD::Save;
    flens::lapack::SVD::Job jobVt = flens::lapack::SVD::None;

    GeneralTreeNode<HTuckerTreeNode<T> > *node;

    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TITgram = gram.getGeneralTree().begin();

    for(; TITgram <= gram.getGeneralTree().end(); TITgram++){
        flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > U;
        flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > Vt;
        flens::DenseVector<flens::Array<T> > s;
        flens::DenseVector<flens::Array<T> > work;

        auto& Gt = TITgram.getNode()->getContent()->getUorB();
        auto  m  = Gt.numRows();
        auto  r  = std::min(m, Gt.numCols());
        U.resize(m, r);
        s.resize(r);

        flens::lapack::svd(jobU, jobVt,
                           Gt, s, U, Vt,
                           work);

        T error                      = 0.;
        FLENS_DEFAULT_INDEXTYPE rank = s.length();
        for (; rank>=1; --rank) {
            error += s(rank);
            if (error > eps || rank == 1) {
                break;
            }
        }

        TITgram.getNode()->getContent()->
            setUorB(U(_,_(1,std::min(U.numRows(), rank))));
    }
    }

    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT     = this->getGeneralTree().begin();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITgram = gram.getGeneralTree().begin();
    for(; TIT <= this->getGeneralTree().end(); TIT++,TITgram++){
        node = TIT.getNode();

        if(node->isLeaf()){
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp;
            flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,node->getContent()->getUorB(),TITgram.getNode()->getContent()->getUorB(),0.0,tmp);
            node->getContent()->setUorB(tmp);
            node->getContent()->setNumRows(tmp.numCols());
        } else {
            int rcnumel = node->getContent()->getRightChildNumRows();
            int lcnumel = node->getContent()->getLeftChildNumRows();

            int numel = node->getContent()->getNumRows();
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Sref = TITgram.getNode()->getContent()->getUorB();
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Slref = TITgram.getNode()->getfirstChild()->getContent()->getUorB();
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Srref = TITgram.getNode()->getlastChild()->getContent()->getUorB();
            int newnumel = Sref.numCols();
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Bref = node->getContent()->getUorB();
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp(rcnumel*newnumel,lcnumel);
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Bt(rcnumel*lcnumel, numel);
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > BtSt(rcnumel*lcnumel, newnumel);

            /* Bt*St */
            for(int i = 1; i <= numel; ++i) {
                for(int j = 1; j <= lcnumel; ++j) {
                    Bt(_((j-1)*rcnumel+1, j*rcnumel), i) = Bref(_((i-1)*rcnumel+1, i*rcnumel), j);
                }
            }

            flens::blas::mm(cxxblas::NoTrans, cxxblas::NoTrans, 1.0, Bt, Sref,0.0, BtSt);
            for(int i = 1; i <= newnumel; ++i) {
                for(int j = 1; j <= lcnumel; ++j) {
                    tmp(_((i-1)*rcnumel+1, i*rcnumel), j) = BtSt(_((j-1)*rcnumel+1, j*rcnumel), i);
                }
            }

            //hier haben wir nun als B_t*S_t berechnet. Fehlt noch (S_tl^T x S_tr^T)*B_t*S_t
            flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > newB(newnumel*Srref.numCols(),Slref.numCols());
            for(int i = 1; i <= newnumel; ++i){
                flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp1;
                flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,tmp(_((i-1)*rcnumel+1,i*rcnumel),_),Slref,0.0,tmp1);
                flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,Srref,tmp1,0.0,newB(_((i-1)*Srref.numCols()+1,i*Srref.numCols()),_));
            }
            node->getContent()->setUorB(newB);
            node->getContent()->setNumRows(newnumel);
            node->getContent()->setLeftChildNumRows(Slref.numCols());
            node->getContent()->setRightChildNumRows(Srref.numCols());
        }
    }
    }
}


template <typename T>
void
HTuckerTree<T>::truncate_hsvd(const T eps)
{
    auto gram = gramians_nonorthogonal(*this);
    int  size = this->depth()-1;
    flens::DenseVector<flens::Array<T> > eps2vec(size);
    T val = eps*eps/(T) size;
    eps2vec.fill(val);
    this->truncate(gram, eps2vec);
}


template <typename T>
void
HTuckerTree<T>::truncate_hsvd_soft(const T eps)
{
    auto gram = gramians_nonorthogonal(*this);
    this->truncate_soft(gram, eps);
}


template <typename T>
void
HTuckerTree<T>::truncate(const HTuckerTree<T>& gram, const int rank){
    using flens::_;

	HTuckerTree<T> Rtree;
    Rtree.set_tree(*this);


	flens::DenseVector<flens::Array<T> > s;
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > U, Vt;
	
	GeneralTreeNode<HTuckerTreeNode<T> > *node;
	
    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = this->getGeneralTree().end();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITR = Rtree.getGeneralTree().end();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITgram = gram.getGeneralTree().end();
    for(;TIT >= this->getGeneralTree().begin(); TIT--,TITR--,TITgram--){
		node = TIT.getNode();
		if(node->isLeaf()){
			//compute QR-Decomposition
			//hier ist die Kopie notwendig, weil qrf die Matrix verändert!
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Uref = node->getContent()->getUorB();
			int m = Uref.numRows();
			int n = Uref.numCols();
			flens::DenseVector<flens::Array<T> > tau;
			qrf(Uref,tau);
			
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Q(m,std::min(m,n));
			Q(_(1,m),_(1,std::min(n,m))) = Uref(_,_(1,std::min(m,n)));
			orgqr(Q,tau);


			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Gref = TITgram.getNode()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp1,tmp2,tmp3,tmp4;

            tmp1 = Gref;

			//hier berechnen wir R_t*G_t*R_t^T   => R_t entsprciht U_ref, wenn m > n dann ist die Matrix höher als Breit und enthält nur 0en, welche rausgeschmissen werden können
			if (m < n){
				//dieser Fall ist komplizierter, weil die Matrix R zerlegt werden muss..
				tmp2 = Gref(_,_(1,m));
				tmp3 = Gref(_,_(m+1,n));
				flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,Uref(_,_(1,m)).upper(),tmp2);
                tmp1.resize(0, 0);
				flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0,tmp3,Uref(_,_(m+1,n)),0.0,tmp1);
                tmp3.resize(0, 0);
				tmp3 = tmp1 + tmp2;
                tmp1.resize(0, 0);
                tmp2.resize(0, 0);
				tmp1 = tmp3(_(1,m),_);
				tmp2 = tmp3(_(m+1,n),_);
				flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,Uref(_,_(1,m)).upper(),tmp1);
                tmp3.resize(0, 0);
				flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Uref(_,_(m+1,n)),tmp2,0.0,tmp3);
                tmp2.resize(0, 0);
				tmp2 = tmp1 + tmp3;
                tmp1.resize(0, 0);
				tmp1 = tmp2;
			} else {
				//d.h. in diesem Fall ist R in n x n und Q in m x n 
				flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,Uref(_(1,n),_).upper(),tmp1);
				flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,Uref(_(1,n),_).upper(),tmp1);
			}

			//flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,Uref(_(1,std::min(m,n)),_).upper(),tmp1);
			//flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,Uref(_(1,std::min(m,n)),_).upper(),tmp1);

			//flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0, Gref,Uref(_(1,std::min(m,n)),_).upper(),0.0,tmp1);
			//flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Uref(_(1,std::min(m,n)),_).upper(),tmp1,0.0,tmp2);

			flens::svd(tmp1,s,U,Vt);
			

			//Die S_t^T*  R_t werden hier transponiert gespeicher! => also R_t^T*S_t
			if ( m < n){
                tmp2.resize(0, 0);
                tmp3.resize(0, 0);
				tmp2 = U(_,_(1,std::min(U.numCols(),rank)));
				tmp3 = tmp2;
				flens::blas::mm(cxxblas::Left,cxxblas::Trans,1.0,Uref(_,_(1,m)).upper(),tmp2);
				flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,Uref(_,_(m+1,n)),tmp3,0.0,tmp4);
                tmp1.resize(0, 0);
				tmp1 = flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> >(tmp2.numRows()+tmp4.numRows(),tmp2.numCols());
				tmp1(_(1,tmp2.numRows()),_) = tmp2;
				tmp1(_(tmp2.numRows()+1,tmp1.numRows()),_) = tmp4;
			} else {
                tmp1.resize(0, 0);
				tmp1 = U(_(1,n),_(1,std::min(U.numCols(),rank)));
				flens::blas::mm(cxxblas::Left,cxxblas::Trans,1.0,Uref(_(1,n),_).upper(),tmp1);

			}
			
			
			//flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,U(_,_(1,std::min(U.numCols(),rank))),Uref(_(1,std::min(m,n)),_).upper(),0.0,tmp1);
			TITR.getNode()->getContent()->setUorB(tmp1);

            tmp1.resize(0, 0);
			flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Q,U(_,_(1,std::min(U.numCols(),rank))),0.0,tmp1);
			node->getContent()->setUorB(tmp1);
			node->getContent()->setNumRows(tmp1.numCols());
		} else if(node->isRoot()){
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &leftMat = TITR.getNode()->getfirstChild()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &rightMat = TITR.getNode()->getlastChild()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Bref = node->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp1, tmp2;

			//numel sollte hier immer 1 sein. Ansonsten genau das selbe wie weiter unten...
			

			flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Bref,leftMat,0.0,tmp1);
			flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,rightMat,tmp1,0.0,tmp2);
			
			//einsortiert werden muss hier nicht, da wir nicht mit QR weiterverarbeiten ist dies schon die Blockstruktur!

			node->getContent()->setUorB(tmp2);
			node->getContent()->setNumRows(1);
			node->getContent()->setLeftChildNumRows(tmp2.numCols());
			node->getContent()->setRightChildNumRows(tmp2.numRows());
			
		} else {
			//compute (Stl^T Rtl x Str^T Rtr) B_t
			int numel = node->getContent()->getNumRows();
			int rcnumel = node->getContent()->getRightChildNumRows();
			
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &leftMat = TITR.getNode()->getfirstChild()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &rightMat = TITR.getNode()->getlastChild()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Bref = node->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > QRmat(leftMat.numCols()*rightMat.numCols(),numel);
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp1,tmp2,tmp3,tmp4;

			int lr = leftMat.numCols(); // wegen leftMat^T
			int rr = rightMat.numCols();
		

			for(int i = 1; i <= numel; ++i){
				flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Bref(_((i-1)*rcnumel + 1,i*rcnumel),_),leftMat,0.0,tmp1);
				flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,rightMat,tmp1,0.0,tmp2);
				//jetzt noch in die neue Matrix einsortieren...
				for(int j = 1; j <= lr; ++j){
					QRmat(_((j-1)*rr+1,j*rr),i) = tmp2(_,j);
				}
			}
			// jetzt die QR Zerlegung anwenden:
			int m = QRmat.numRows();
			int n = QRmat.numCols();
			flens::DenseVector<flens::Array<T> > tau;
			qrf(QRmat,tau);
			
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Q(m,std::min(m,n));
			Q(_(1,m),_) = QRmat(_,_(1,std::min(m,n)));
			orgqr(Q,tau);

			// Dann wie oben R_t G_t R_t^T berechnen
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Gref = TITgram.getNode()->getContent()->getUorB();

            tmp1.resize(0, 0);
            tmp1 = Gref;
			if (m < n){
				//dieser Fall ist komplizierter, weil die Matrix R zerlegt werden muss..
                tmp2.resize(0, 0);
                tmp3.resize(0, 0);
				tmp2 = Gref(_,_(1,m));
				tmp3 = Gref(_,_(m+1,n));
				flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,QRmat(_,_(1,m)).upper(),tmp2);
                tmp1.resize(0, 0);
				flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0,tmp3,QRmat(_,_(m+1,n)),0.0,tmp1);
                tmp3.resize(0, 0);
				tmp3 = tmp1 + tmp2;
                tmp1.resize(0, 0);
				tmp1 = tmp3(_(1,m),_);
                tmp2.resize(0, 0);
				tmp2 = tmp3(_(m+1,n),_);
				flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,QRmat(_,_(1,m)).upper(),tmp1);
                tmp3.resize(0, 0);
				flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,QRmat(_,_(m+1,n)),tmp2,0.0,tmp3);
                tmp2.resize(0, 0);
				tmp2 = tmp1 + tmp3;
                tmp1.resize(0, 0);
				tmp1 = tmp2;
			} else {
				//d.h. in diesem Fall ist R in n x n und Q in m x n 
				flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,QRmat(_(1,n),_).upper(),tmp1);
				flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,QRmat(_(1,n),_).upper(),tmp1);
			}
			//flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0, Gref,QRmat(_(1,std::min(m,n)),_).upper(),0.0,tmp1);
			//flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,QRmat(_(1,std::min(m,n)),_).upper(),tmp1,0.0,tmp2);

			//Svd wie oben

			flens::svd(tmp1,s,U,Vt);
			//tmp1 = U(_,_(1,std::min(U.numCols(),rank)));
			if ( m < n){
                tmp2.resize(0, 0);
				tmp2 = U(_,_(1,std::min(U.numCols(),rank)));
                tmp3.resize(0, 0);
				tmp3 = tmp2;
				flens::blas::mm(cxxblas::Left,cxxblas::Trans,1.0,QRmat(_,_(1,m)).upper(),tmp2);
				flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,QRmat(_,_(m+1,n)),tmp3,0.0,tmp4);
                tmp1.resize(0, 0);
				tmp1 = flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> >(tmp2.numRows()+tmp4.numRows(),tmp2.numCols());
				tmp1(_(1,tmp2.numRows()),_) = tmp2;
				tmp1(_(tmp2.numRows()+1,tmp1.numRows()),_) = tmp4;
			} else {
                tmp1.resize(0, 0);
				tmp1 = U(_(1,n),_(1,std::min(U.numCols(),rank)));
				flens::blas::mm(cxxblas::Left,cxxblas::Trans,1.0,QRmat(_(1,n),_).upper(),tmp1);

			}
			//flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,U(_,_(1,std::min(U.numCols(),rank))),Uref(_(1,std::min(m,n)),_).upper(),0.0,tmp1);
			TITR.getNode()->getContent()->setUorB(tmp1);

			//Das hier ist Q_t*S_t, die muss nun zurück ins Blockformat!!!
            tmp1.resize(0, 0);
			flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Q,U(_,_(1,std::min(U.numCols(),rank))),0.0,tmp1);
			int newnumel = tmp1.numCols();
			int newlcnumel = node->getfirstChild()->getContent()->getNumRows();
			int newrcnumel = node->getlastChild()->getContent()->getNumRows();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > newB(newnumel*newrcnumel,newlcnumel);

			for(int i = 1; i <= newnumel; ++i){
				for(int j = 1; j <= newlcnumel; ++j){
					newB(_((i-1)*newrcnumel+1,i*newrcnumel),j) = tmp1(_((j-1)*newrcnumel + 1,j*newrcnumel) ,i);
				}
			}

			node->getContent()->setUorB(newB);
			node->getContent()->setNumRows(newnumel);
			node->getContent()->setLeftChildNumRows(newlcnumel);
			node->getContent()->setRightChildNumRows(newrcnumel);
		}

		
	
	}
    }
}


template <typename T>
void
HTuckerTree<T>::truncate(const HTuckerTree<T>& gram,
                         const flens::DenseVector<flens::Array<T> >& eps)
{
    using flens::_;

    HTuckerTree<T> Rtree;
    Rtree.set_tree(*this);

    typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Matrix;
    typedef flens::DenseVector<flens::Array<T> >                      Vec;
    typedef flens::DenseVector<flens::Array<int> >                    IntVec;


    Vec s;
    Matrix U, Vt;

    std::vector<Vec> svals;
    std::vector<int> ranks;

    std::vector<Matrix> Us;
    std::vector<Matrix> Urefs;
    std::vector<Matrix> Qs;
    std::vector<Matrix> QRmats;


    GeneralTreeNode<HTuckerTreeNode<T> > *node;

    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = this->getGeneralTree().end();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITR = Rtree.getGeneralTree().end();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITgram = gram.getGeneralTree().end();

    GeneralTreeIterator<HTuckerTreeNode<T> > save;
    GeneralTreeIterator<HTuckerTreeNode<T> > saver;
    GeneralTreeIterator<HTuckerTreeNode<T> > savegram;

    assert(eps.length()==TIT.getNode()->level()-1);

    int  depth    = TIT.getNode()->level();
    int  level    = depth;
    bool sorted   = false;
    bool newlevel = true;
    int  counter  = 0;
    int  length   = 0;

    for(;TIT >= this->getGeneralTree().begin();){
        node = TIT.getNode();

        if (node->level()==level) {

            if (newlevel) {
                save     = TIT;
                saver    = TITR;
                savegram = TITgram;

                Us.clear();
                Urefs.clear();
                Qs.clear();
                QRmats.clear();
                svals.clear();
                ranks.clear();

                length   = 0;
                sorted   = false;
                newlevel = false;
            }

            if (node->isLeaf()) {
                //compute QR-Decomposition
                //hier ist die Kopie notwendig, weil qrf die Matrix verändert!
                Matrix Uref = node->getContent()->getUorB();
                int m = Uref.numRows();
                int n = Uref.numCols();
                Vec tau;

                qrf(Uref,tau);
                Matrix Q(m,std::min(m,n));
                Q = Uref(_,_(1,std::min(m,n)));
                orgqr(Q, tau);

                Matrix &Gref = TITgram.getNode()->getContent()->getUorB();
                Matrix tmp1,tmp2,tmp3,tmp4;

                tmp1 = Gref;

                //hier berechnen wir R_t*G_t*R_t^T   => R_t entsprciht U_ref, wenn m > n dann ist die Matrix höher als Breit und enthält nur 0en, welche rausgeschmissen werden können
                if (m < n){
                    //dieser Fall ist komplizierter, weil die Matrix R zerlegt werden muss..
                    tmp2 = Gref(_,_(1,m));
                    tmp3 = Gref(_,_(m+1,n));
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,Uref(_,_(1,m)).upper(),tmp2);
                    tmp1.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0,tmp3,Uref(_,_(m+1,n)),0.0,tmp1);
                    tmp3.resize(0, 0);
                    tmp3 = tmp1 + tmp2;
                    tmp1.resize(0, 0);
                    tmp2.resize(0, 0);
                    tmp1 = tmp3(_(1,m),_);
                    tmp2 = tmp3(_(m+1,n),_);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,Uref(_,_(1,m)).upper(),tmp1);
                    tmp3.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Uref(_,_(m+1,n)),tmp2,0.0,tmp3);
                    tmp2.resize(0, 0);
                    tmp1 = tmp1 + tmp3;
                } else {
                    //d.h. in diesem Fall ist R in n x n und Q in m x n
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,Uref(_(1,n),_).upper(),tmp1);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,Uref(_(1,n),_).upper(),tmp1);
                }

                flens::svd(tmp1,s,U,Vt);

                Us.push_back(U);
                Urefs.push_back(Uref);
                Qs.push_back(Q);
                svals.push_back(s);
                length += s.length();
            } else if (node->isRoot()) {
                sorted = true;
                --level;
                continue;
            } else {
                //compute (Stl^T Rtl x Str^T Rtr) B_t
                int numel        = node->getContent()->getNumRows();
                int rcnumel      = node->getContent()->getRightChildNumRows();

                Matrix &leftMat  = TITR.getNode()->getfirstChild()->getContent()->getUorB();
                Matrix &rightMat = TITR.getNode()->getlastChild()->getContent()->getUorB();
                Matrix &Bref     = node->getContent()->getUorB();
                Matrix QRmat(leftMat.numCols()*rightMat.numCols(),numel);
                Matrix tmp1,tmp2,tmp3,tmp4;

                int lr = leftMat.numCols(); // wegen leftMat^T
                int rr = rightMat.numCols();


                for(int i = 1; i <= numel; ++i){
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Bref(_((i-1)*rcnumel + 1,i*rcnumel),_),leftMat,0.0,tmp1);
                    flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,rightMat,tmp1,0.0,tmp2);
                    //jetzt noch in die neue Matrix einsortieren...
                    for(int j = 1; j <= lr; ++j){
                        QRmat(_((j-1)*rr+1,j*rr),i) = tmp2(_,j);
                    }
                }
                // jetzt die QR Zerlegung anwenden:
                int m = QRmat.numRows();
                int n = QRmat.numCols();
                Vec tau;
                qrf(QRmat,tau);

                Matrix Q(m,std::min(m,n));
                Q = QRmat(_,_(1,std::min(m,n)));
                orgqr(Q,tau);

                // Dann wie oben R_t G_t R_t^T berechnen
                Matrix &Gref = TITgram.getNode()->getContent()->getUorB();

                tmp1.resize(0, 0);
                tmp1 = Gref;
                if (m < n){
                    //dieser Fall ist komplizierter, weil die Matrix R zerlegt werden muss..
                    tmp2.resize(0, 0);
                    tmp3.resize(0, 0);
                    tmp2 = Gref(_,_(1,m));
                    tmp3 = Gref(_,_(m+1,n));
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,QRmat(_,_(1,m)).upper(),tmp2);
                    tmp1.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0,tmp3,QRmat(_,_(m+1,n)),0.0,tmp1);
                    tmp3.resize(0, 0);
                    tmp3 = tmp1 + tmp2;
                    tmp1.resize(0, 0);
                    tmp1 = tmp3(_(1,m),_);
                    tmp2.resize(0, 0);
                    tmp2 = tmp3(_(m+1,n),_);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,QRmat(_,_(1,m)).upper(),tmp1);
                    tmp3.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,QRmat(_,_(m+1,n)),tmp2,0.0,tmp3);
                    tmp2.resize(0, 0);
                    tmp1 = tmp1 + tmp3;
                } else {
                    //d.h. in diesem Fall ist R in n x n und Q in m x n 
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,QRmat(_(1,n),_).upper(),tmp1);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,QRmat(_(1,n),_).upper(),tmp1);
                }

                //Svd wie oben
                flens::svd(tmp1,s,U,Vt);

                Us.push_back(U);
                Qs.push_back(Q);
                QRmats.push_back(QRmat);
                svals.push_back(s);
                length += s.length();
            }

            TIT--;
            TITR--;
            TITgram--;

        } else if (!sorted) {
            TIT     = save;
            TITR    = saver;
            TITgram = savegram;
            sorted  = true;
            counter = 0;

            Vec s_all(length);
            int off = 1;
            for (std::size_t i=0; i<svals.size(); ++i) {
                s_all(_(off,off+svals[i].length()-1)) = svals[i];
                off += svals[i].length();
            }


            IntVec rho(s_all.length());
            flens::sort(s_all, rho);
            double delta = 0.;
            double eps2  = eps(depth-level+1);
            off          = 0;
            for (; off<s_all.length(); ++off) {
                delta += s_all(off+1);
                if (delta>eps2) break;
            }

            IntVec dels(svals.size());
            for (int i=1; i<=off; ++i) {
                int right = 0;
                for (std::size_t j=0; j<svals.size(); ++j) {
                    right += svals[j].length();
                    if (rho(i)<=right) {
                        ++dels(j+1);
                        break;
                    }
                }
            }

            ranks.resize(svals.size());
            for (std::size_t i=0; i<svals.size(); ++i) {
                int rank = svals[i].length() - dels(i+1);
                if (rank==0) rank = 1;
                ranks[i] = rank;
            }

            --level;
        } else {
            if (node->isLeaf()) {

                Matrix tmp1, tmp2, tmp3, tmp4;
                Matrix& Uref = Urefs[counter];
                Matrix& U    = Us[counter];
                Matrix& Q    = Qs[counter];
                int rank     = ranks[counter];
                int m        = node->getContent()->getUorB().numRows();
                int n        = node->getContent()->getUorB().numCols();;

                if ( m < n){
                    tmp2 = U(_,_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank)));
                    tmp3 = tmp2;
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    Uref(_,_(1,m)).upper(),tmp2);
                    flens::blas::mm(cxxblas::Trans,
                                    cxxblas::NoTrans, 1.0,
                                    Uref(_,_(m+1,n)),tmp3,0.0,tmp4);
                    tmp1.resize(0, 0);
                    tmp1 = Matrix(tmp2.numRows()+tmp4.numRows(),tmp2.numCols());
                    tmp1(_(1,tmp2.numRows()),_) = tmp2;
                    tmp1(_(tmp2.numRows()+1,tmp1.numRows()),_) = tmp4;
                } else {
                    tmp1 = U
                           (_(1,n),_(1,std::min(U.numCols(),
                                                (FLENS_DEFAULT_INDEXTYPE)rank)));
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    Uref(_(1,n),_).upper(),tmp1);

                }

                TITR.getNode()->getContent()->setUorB(tmp1);

                tmp1.resize(0, 0);
                flens::blas::mm(cxxblas::NoTrans,
                                cxxblas::NoTrans, 1.0,
                                Q,
                                U
                                (_,_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank))),0.0,tmp1);
                node->getContent()->setUorB(tmp1);
                node->getContent()->setNumRows(tmp1.numCols());
            } else if (node->isRoot()) {

                Matrix &leftMat = TITR.getNode()->getfirstChild()->getContent()->getUorB();
                Matrix &rightMat = TITR.getNode()->getlastChild()->getContent()->getUorB();
                Matrix &Bref = node->getContent()->getUorB();
                Matrix tmp1, tmp2;

                //numel sollte hier immer 1 sein. Ansonsten genau das selbe wie weiter unten...

                flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Bref,leftMat,0.0,tmp1);
                flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,rightMat,tmp1,0.0,tmp2);

                //einsortiert werden muss hier nicht, da wir nicht mit QR weiterverarbeiten ist dies schon die Blockstruktur!

                node->getContent()->setUorB(tmp2);
                node->getContent()->setNumRows(1);
                node->getContent()->setLeftChildNumRows(tmp2.numCols());
                node->getContent()->setRightChildNumRows(tmp2.numRows());
            } else {
                int numel        = node->getContent()->getNumRows();
                Matrix &leftMat  = TITR.getNode()->getfirstChild()
                                   ->getContent()->getUorB();
                Matrix &rightMat = TITR.getNode()->getlastChild()
                                   ->getContent()->getUorB();


                Matrix tmp1, tmp2, tmp3, tmp4;
                Matrix& QRmat = QRmats[counter];
                Matrix& U     = Us[counter];
                Matrix& Q     = Qs[counter];
                int rank      = ranks[counter];
                int m         = leftMat.numCols()*rightMat.numCols();
                int n         = numel;

                if ( m < n){
                    tmp2 = U(_,_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank)));
                    tmp3 = tmp2;
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    QRmat(_,_(1,m)).upper(),tmp2);
                    flens::blas::mm(cxxblas::Trans,
                                    cxxblas::NoTrans, 1.0,
                                    QRmat(_,_(m+1,n)),tmp3,0.0,tmp4);
                    tmp1 = Matrix(tmp2.numRows()+tmp4.numRows(),tmp2.numCols());
                    tmp1(_(1,tmp2.numRows()),_) = tmp2;
                    tmp1(_(tmp2.numRows()+1,tmp1.numRows()),_) = tmp4;
                } else {
                    tmp1 = U(_(1,n),_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank)));
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    QRmat(_(1,n),_).upper(),tmp1);

                }

                TITR.getNode()->getContent()->setUorB(tmp1);

                //Das hier ist Q_t*S_t, die muss nun zurück ins Blockformat!!!
                tmp1.resize(0, 0);
                flens::blas::mm(cxxblas::NoTrans,
                                cxxblas::NoTrans, 1.0,
                                Q,
                                U(_,_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank))),0.0,tmp1);
                int newnumel   = tmp1.numCols();
                int newlcnumel = node->getfirstChild()->getContent()->getNumRows();
                int newrcnumel = node->getlastChild()->getContent()->getNumRows();
                Matrix newB(newnumel*newrcnumel,newlcnumel);

                for(int i = 1; i <= newnumel; ++i){
                    for(int j = 1; j <= newlcnumel; ++j){
                        newB(_((i-1)*newrcnumel+1,i*newrcnumel),j) = tmp1(_((j-1)*newrcnumel + 1,j*newrcnumel) ,i);
                    }
                }

                node->getContent()->setUorB(newB);
                node->getContent()->setNumRows(newnumel);
                node->getContent()->setLeftChildNumRows(newlcnumel);
                node->getContent()->setRightChildNumRows(newrcnumel);
            }

            TIT--;
            TITR--;
            TITgram--;
            ++counter;
            newlevel = true;
        }
    }
}


template <typename T>
void
HTuckerTree<T>::truncate_soft(const HTuckerTree<T>& gram,
                              const double eps)
{
    using flens::_;

    HTuckerTree<T> Rtree;
    Rtree.set_tree(*this);

    typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Matrix;
    typedef flens::DenseVector<flens::Array<T> >                      Vec;

    Vec s;
    Matrix U, Vt;

    std::vector<Vec> svals;
    std::vector<Vec> scals;
    std::vector<int> ranks;

    std::vector<Matrix> Us;
    std::vector<Matrix> Urefs;
    std::vector<Matrix> Qs;
    std::vector<Matrix> QRmats;


    GeneralTreeNode<HTuckerTreeNode<T> > *node;

    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = this->getGeneralTree().end();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITR = Rtree.getGeneralTree().end();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITgram = gram.getGeneralTree().end();

    GeneralTreeIterator<HTuckerTreeNode<T> > save;
    GeneralTreeIterator<HTuckerTreeNode<T> > saver;
    GeneralTreeIterator<HTuckerTreeNode<T> > savegram;

    double alpha2 = eps*eps;

    int  depth    = TIT.getNode()->level();
    int  level    = depth;
    bool sorted   = false;
    bool newlevel = true;
    int  counter  = 0;
    int  length   = 0;

    for(;TIT >= this->getGeneralTree().begin();){
        node = TIT.getNode();

        if (node->level()==level) {

            if (newlevel) {
                save     = TIT;
                saver    = TITR;
                savegram = TITgram;

                Us.clear();
                Urefs.clear();
                Qs.clear();
                QRmats.clear();
                svals.clear();
                scals.clear();
                ranks.clear();

                length   = 0;
                sorted   = false;
                newlevel = false;
            }

            if (node->isLeaf()) {
                //compute QR-Decomposition
                //hier ist die Kopie notwendig, weil qrf die Matrix verändert!
                Matrix Uref = node->getContent()->getUorB();
                int m = Uref.numRows();
                int n = Uref.numCols();
                Vec tau;

                qrf(Uref,tau);
                Matrix Q(m,std::min(m,n));
                Q = Uref(_,_(1,std::min(m,n)));
                orgqr(Q, tau);

                Matrix &Gref = TITgram.getNode()->getContent()->getUorB();
                Matrix tmp1,tmp2,tmp3,tmp4;

                tmp1 = Gref;

                //hier berechnen wir R_t*G_t*R_t^T   => R_t entsprciht U_ref, wenn m > n dann ist die Matrix höher als Breit und enthält nur 0en, welche rausgeschmissen werden können
                if (m < n){
                    //dieser Fall ist komplizierter, weil die Matrix R zerlegt werden muss..
                    tmp2 = Gref(_,_(1,m));
                    tmp3 = Gref(_,_(m+1,n));
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,Uref(_,_(1,m)).upper(),tmp2);
                    tmp1.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0,tmp3,Uref(_,_(m+1,n)),0.0,tmp1);
                    tmp3.resize(0, 0);
                    tmp3 = tmp1 + tmp2;
                    tmp1.resize(0, 0);
                    tmp2.resize(0, 0);
                    tmp1 = tmp3(_(1,m),_);
                    tmp2 = tmp3(_(m+1,n),_);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,Uref(_,_(1,m)).upper(),tmp1);
                    tmp3.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Uref(_,_(m+1,n)),tmp2,0.0,tmp3);
                    tmp2.resize(0, 0);
                    tmp1 = tmp1 + tmp3;
                } else {
                    //d.h. in diesem Fall ist R in n x n und Q in m x n
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,Uref(_(1,n),_).upper(),tmp1);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,Uref(_(1,n),_).upper(),tmp1);
                }

                flens::svd(tmp1,s,U,Vt);

                Us.push_back(U);
                Urefs.push_back(Uref);
                Qs.push_back(Q);
                svals.push_back(s);
                scals.push_back(s);
                length += s.length();
            } else if (node->isRoot()) {
                sorted = true;
                --level;
                continue;
            } else {
                //compute (Stl^T Rtl x Str^T Rtr) B_t
                int numel        = node->getContent()->getNumRows();
                int rcnumel      = node->getContent()->getRightChildNumRows();

                Matrix &leftMat  = TITR.getNode()->getfirstChild()->getContent()->getUorB();
                Matrix &rightMat = TITR.getNode()->getlastChild()->getContent()->getUorB();
                Matrix &Bref     = node->getContent()->getUorB();
                Matrix QRmat(leftMat.numCols()*rightMat.numCols(),numel);
                Matrix tmp1,tmp2,tmp3,tmp4;

                int lr = leftMat.numCols(); // wegen leftMat^T
                int rr = rightMat.numCols();


                for(int i = 1; i <= numel; ++i){
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Bref(_((i-1)*rcnumel + 1,i*rcnumel),_),leftMat,0.0,tmp1);
                    flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,rightMat,tmp1,0.0,tmp2);
                    //jetzt noch in die neue Matrix einsortieren...
                    for(int j = 1; j <= lr; ++j){
                        QRmat(_((j-1)*rr+1,j*rr),i) = tmp2(_,j);
                    }
                }
                // jetzt die QR Zerlegung anwenden:
                int m = QRmat.numRows();
                int n = QRmat.numCols();
                Vec tau;
                qrf(QRmat,tau);

                Matrix Q(m,std::min(m,n));
                Q = QRmat(_,_(1,std::min(m,n)));
                orgqr(Q,tau);

                // Dann wie oben R_t G_t R_t^T berechnen
                Matrix &Gref = TITgram.getNode()->getContent()->getUorB();

                tmp1.resize(0, 0);
                tmp1 = Gref;
                if (m < n){
                    //dieser Fall ist komplizierter, weil die Matrix R zerlegt werden muss..
                    tmp2.resize(0, 0);
                    tmp3.resize(0, 0);
                    tmp2 = Gref(_,_(1,m));
                    tmp3 = Gref(_,_(m+1,n));
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,QRmat(_,_(1,m)).upper(),tmp2);
                    tmp1.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0,tmp3,QRmat(_,_(m+1,n)),0.0,tmp1);
                    tmp3.resize(0, 0);
                    tmp3 = tmp1 + tmp2;
                    tmp1.resize(0, 0);
                    tmp1 = tmp3(_(1,m),_);
                    tmp2.resize(0, 0);
                    tmp2 = tmp3(_(m+1,n),_);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,QRmat(_,_(1,m)).upper(),tmp1);
                    tmp3.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,QRmat(_,_(m+1,n)),tmp2,0.0,tmp3);
                    tmp2.resize(0, 0);
                    tmp1 = tmp1 + tmp3;
                } else {
                    //d.h. in diesem Fall ist R in n x n und Q in m x n 
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,QRmat(_(1,n),_).upper(),tmp1);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,QRmat(_(1,n),_).upper(),tmp1);
                }

                //Svd wie oben
                flens::svd(tmp1,s,U,Vt);

                Us.push_back(U);
                Qs.push_back(Q);
                QRmats.push_back(QRmat);
                svals.push_back(s);
                scals.push_back(s);
                length += s.length();
            }

            TIT--;
            TITR--;
            TITgram--;

        } else if (!sorted) {
            TIT     = save;
            TITR    = saver;
            TITgram = savegram;
            sorted  = true;
            counter = 0;

            ranks.resize(svals.size());
            for (std::size_t i=0; i<svals.size(); ++i) {
                ranks[i] = 0;
                for (int j=svals[i].length(); j>=1; --j) {
                    if (svals[i](j)>alpha2) break;
                    else ++ranks[i];
                }
                ranks[i] = svals[i].length()-ranks[i];
                ranks[i] = std::max(ranks[i], 1);
            }

            for (std::size_t i=0; i<svals.size(); ++i) {
                for (int j=1; j<=ranks[i]; ++j) {
                    auto sigma = std::sqrt(svals[i](j));
                    scals[i](j) = (sigma-eps)/sigma;
                }
            }

            --level;
        } else {
            if (node->isLeaf()) {

                Matrix tmp1, tmp2, tmp3, tmp4;
                Matrix& Uref = Urefs[counter];
                Matrix& U    = Us[counter];
                Matrix& Q    = Qs[counter];
                int rank     = ranks[counter];
                int m        = node->getContent()->getUorB().numRows();
                int n        = node->getContent()->getUorB().numCols();
                Vec ks       = scals[counter];

                for (int i=1; i<=rank; ++i) {
                    flens::blas::scal(ks(i), U(_, i));
                }

                if ( m < n){
                    tmp2 = U(_,_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank)));
                    tmp3 = tmp2;
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    Uref(_,_(1,m)).upper(),tmp2);
                    flens::blas::mm(cxxblas::Trans,
                                    cxxblas::NoTrans, 1.0,
                                    Uref(_,_(m+1,n)),tmp3,0.0,tmp4);
                    tmp1.resize(0, 0);
                    tmp1 = Matrix(tmp2.numRows()+tmp4.numRows(),tmp2.numCols());
                    tmp1(_(1,tmp2.numRows()),_) = tmp2;
                    tmp1(_(tmp2.numRows()+1,tmp1.numRows()),_) = tmp4;
                } else {
                    tmp1 = U
                           (_(1,n),_(1,std::min(U.numCols(),
                                                (FLENS_DEFAULT_INDEXTYPE)rank)));
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    Uref(_(1,n),_).upper(),tmp1);

                }

                TITR.getNode()->getContent()->setUorB(tmp1);

                tmp1.resize(0, 0);
                flens::blas::mm(cxxblas::NoTrans,
                                cxxblas::NoTrans, 1.0,
                                Q,
                                U
                                (_,_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank))),0.0,tmp1);
                node->getContent()->setUorB(tmp1);
                node->getContent()->setNumRows(tmp1.numCols());
            } else if (node->isRoot()) {

                Matrix &leftMat = TITR.getNode()->getfirstChild()->getContent()->getUorB();
                Matrix &rightMat = TITR.getNode()->getlastChild()->getContent()->getUorB();
                Matrix &Bref = node->getContent()->getUorB();
                Matrix tmp1, tmp2;

                //numel sollte hier immer 1 sein. Ansonsten genau das selbe wie weiter unten...

                flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Bref,leftMat,0.0,tmp1);
                flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,rightMat,tmp1,0.0,tmp2);

                //einsortiert werden muss hier nicht, da wir nicht mit QR weiterverarbeiten ist dies schon die Blockstruktur!

                node->getContent()->setUorB(tmp2);
                node->getContent()->setNumRows(1);
                node->getContent()->setLeftChildNumRows(tmp2.numCols());
                node->getContent()->setRightChildNumRows(tmp2.numRows());
            } else {
                int numel        = node->getContent()->getNumRows();
                Matrix &leftMat  = TITR.getNode()->getfirstChild()
                                   ->getContent()->getUorB();
                Matrix &rightMat = TITR.getNode()->getlastChild()
                                   ->getContent()->getUorB();


                Matrix tmp1, tmp2, tmp3, tmp4;
                Matrix& QRmat = QRmats[counter];
                Matrix& U     = Us[counter];
                Matrix& Q     = Qs[counter];
                int rank      = ranks[counter];
                int m         = leftMat.numCols()*rightMat.numCols();
                int n         = numel;
                Vec ks        = scals[counter];

                for (int i=1; i<=rank; ++i) {
                    flens::blas::scal(ks(i), U(_, i));
                }

                if ( m < n){
                    tmp2 = U(_,_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank)));
                    tmp3 = tmp2;
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    QRmat(_,_(1,m)).upper(),tmp2);
                    flens::blas::mm(cxxblas::Trans,
                                    cxxblas::NoTrans, 1.0,
                                    QRmat(_,_(m+1,n)),tmp3,0.0,tmp4);
                    tmp1 = Matrix(tmp2.numRows()+tmp4.numRows(),tmp2.numCols());
                    tmp1(_(1,tmp2.numRows()),_) = tmp2;
                    tmp1(_(tmp2.numRows()+1,tmp1.numRows()),_) = tmp4;
                } else {
                    tmp1 = U(_(1,n),_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank)));
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    QRmat(_(1,n),_).upper(),tmp1);

                }

                TITR.getNode()->getContent()->setUorB(tmp1);

                //Das hier ist Q_t*S_t, die muss nun zurück ins Blockformat!!!
                tmp1.resize(0, 0);
                flens::blas::mm(cxxblas::NoTrans,
                                cxxblas::NoTrans, 1.0,
                                Q,
                                U(_,_(1,std::min(U.numCols(),(FLENS_DEFAULT_INDEXTYPE)rank))),0.0,tmp1);
                int newnumel   = tmp1.numCols();
                int newlcnumel = node->getfirstChild()->getContent()->getNumRows();
                int newrcnumel = node->getlastChild()->getContent()->getNumRows();
                Matrix newB(newnumel*newrcnumel,newlcnumel);

                for(int i = 1; i <= newnumel; ++i){
                    for(int j = 1; j <= newlcnumel; ++j){
                        newB(_((i-1)*newrcnumel+1,i*newrcnumel),j) = tmp1(_((j-1)*newrcnumel + 1,j*newrcnumel) ,i);
                    }
                }

                node->getContent()->setUorB(newB);
                node->getContent()->setNumRows(newnumel);
                node->getContent()->setLeftChildNumRows(newlcnumel);
                node->getContent()->setRightChildNumRows(newrcnumel);
            }

            TIT--;
            TITR--;
            TITgram--;
            ++counter;
            newlevel = true;
        }
    }
}


template <typename T>
void
HTuckerTree<T>::truncate_elementary(const HTuckerTree<T>& gram,
                                    const flens::DenseVector
                                          <flens::Array<T> >& eps)
{
    using flens::_;

    HTuckerTree<T> Rtree;
    Rtree.set_tree(*this);

    typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Matrix;
    typedef flens::DenseVector<flens::Array<T> >                      Vec;
    typedef flens::DenseVector<flens::Array<int> >                    IntVec;


    Vec s;
    Matrix U, Vt;

    std::vector<Vec> svals;
    std::vector<int> ranks;

    std::vector<Matrix> Us;
    std::vector<Matrix> Urefs;
    std::vector<Matrix> Qs;
    std::vector<Matrix> QRmats;


    GeneralTreeNode<HTuckerTreeNode<T> > *node;

    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = this->getGeneralTree().end();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITR = Rtree.getGeneralTree().end();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITgram = gram.getGeneralTree().end();

    GeneralTreeIterator<HTuckerTreeNode<T> > save;
    GeneralTreeIterator<HTuckerTreeNode<T> > saver;
    GeneralTreeIterator<HTuckerTreeNode<T> > savegram;

    assert(eps.length()==TIT.getNode()->level()-1);

    int  depth    = TIT.getNode()->level();
    int  level    = depth;
    bool sorted   = false;
    bool newlevel = true;
    int  counter  = 0;
    int  length   = 0;
    for(;TIT >= this->getGeneralTree().begin();){
        node = TIT.getNode();

        if (node->level()==level) {

            if (newlevel) {
                save     = TIT;
                saver    = TITR;
                savegram = TITgram;

                Us.clear();
                Urefs.clear();
                Qs.clear();
                QRmats.clear();
                svals.clear();
                ranks.clear();

                length   = 0;
                sorted   = false;
                newlevel = false;
            }

            if (node->isLeaf()) {
                //compute QR-Decomposition
                //hier ist die Kopie notwendig, weil qrf die Matrix verändert!
                Matrix Uref = node->getContent()->getUorB();
                int m = Uref.numRows();
                int n = Uref.numCols();
                Vec tau;
                qrf(Uref,tau);

                Matrix Q(m,std::min(m,n));
                Q = Uref(_,_(1,std::min(m,n)));
                orgqr(Q, tau);


                Matrix &Gref = TITgram.getNode()->getContent()->getUorB();
                Matrix tmp1,tmp2,tmp3,tmp4;

                tmp1 = Gref;

                //hier berechnen wir R_t*G_t*R_t^T   => R_t entsprciht U_ref, wenn m > n dann ist die Matrix höher als Breit und enthält nur 0en, welche rausgeschmissen werden können
                if (m < n){
                    //dieser Fall ist komplizierter, weil die Matrix R zerlegt werden muss..
                    tmp2 = Gref(_,_(1,m));
                    tmp3 = Gref(_,_(m+1,n));
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,Uref(_,_(1,m)).upper(),tmp2);
                    tmp1.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0,tmp3,Uref(_,_(m+1,n)),0.0,tmp1);
                    tmp3.resize(0, 0);
                    tmp3 = tmp1 + tmp2;
                    tmp1.resize(0, 0);
                    tmp2.resize(0, 0);
                    tmp1 = tmp3(_(1,m),_);
                    tmp2 = tmp3(_(m+1,n),_);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,Uref(_,_(1,m)).upper(),tmp1);
                    tmp3.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Uref(_,_(m+1,n)),tmp2,0.0,tmp3);
                    tmp2.resize(0, 0);
                    tmp1 = tmp1 + tmp3;
                } else {
                    //d.h. in diesem Fall ist R in n x n und Q in m x n
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,Uref(_(1,n),_).upper(),tmp1);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,Uref(_(1,n),_).upper(),tmp1);
                }

                flens::svd(tmp1,s,U,Vt);

                Us.push_back(U);
                Urefs.push_back(Uref);
                Qs.push_back(Q);
                svals.push_back(s);
                length += s.length();
            } else if (node->isRoot()) {
                sorted = true;
                --level;
                continue;
            } else {
                //compute (Stl^T Rtl x Str^T Rtr) B_t
                int numel        = node->getContent()->getNumRows();

                Matrix &leftMat  = TITR.getNode()->getfirstChild()->getContent()->getUorB();
                Matrix &rightMat = TITR.getNode()->getlastChild()->getContent()->getUorB();
                Matrix QRmat(leftMat.numCols()*rightMat.numCols(),numel);
                Matrix tmp1,tmp2,tmp3,tmp4;

                int lr = leftMat.numCols(); // wegen leftMat^T
                int rr = rightMat.numCols();


                for(int i = 1; i <= numel; ++i){
                    tmp2.resize(rightMat.numCols(), leftMat.numCols());
                    for (int i1=1; i1<=rightMat.numCols(); ++i1) {
                        for (int i2=1; i2<=leftMat.numCols(); ++i2) {
                            tmp2(i1, i2) = rightMat(i, i1)*leftMat(i, i2);
                        }
                    }
                    //jetzt noch in die neue Matrix einsortieren...
                    for(int j = 1; j <= lr; ++j){
                        QRmat(_((j-1)*rr+1,j*rr),i) = tmp2(_,j);
                    }
                }
                // jetzt die QR Zerlegung anwenden:
                int m = QRmat.numRows();
                int n = QRmat.numCols();
                Vec tau;
                qrf(QRmat,tau);

                Matrix Q(m,std::min(m,n));
                Q = QRmat(_,_(1,std::min(m,n)));
                orgqr(Q,tau);

                // Dann wie oben R_t G_t R_t^T berechnen
                Matrix &Gref = TITgram.getNode()->getContent()->getUorB();

                tmp1.resize(0, 0);
                tmp1 = Gref;
                if (m < n){
                    //dieser Fall ist komplizierter, weil die Matrix R zerlegt werden muss..
                    tmp2.resize(0, 0);
                    tmp3.resize(0, 0);
                    tmp2 = Gref(_,_(1,m));
                    tmp3 = Gref(_,_(m+1,n));
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,QRmat(_,_(1,m)).upper(),tmp2);
                    tmp1.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0,tmp3,QRmat(_,_(m+1,n)),0.0,tmp1);
                    tmp3.resize(0, 0);
                    tmp3 = tmp1 + tmp2;
                    tmp1.resize(0, 0);
                    tmp1 = tmp3(_(1,m),_);
                    tmp2.resize(0, 0);
                    tmp2 = tmp3(_(m+1,n),_);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,QRmat(_,_(1,m)).upper(),tmp1);
                    tmp3.resize(0, 0);
                    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,QRmat(_,_(m+1,n)),tmp2,0.0,tmp3);
                    tmp2.resize(0, 0);
                    tmp1 = tmp1 + tmp3;
                } else {
                    //d.h. in diesem Fall ist R in n x n und Q in m x n 
                    flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,QRmat(_(1,n),_).upper(),tmp1);
                    flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,QRmat(_(1,n),_).upper(),tmp1);
                }

                //Svd wie oben
                flens::svd(tmp1,s,U,Vt);

                Us.push_back(U);
                Qs.push_back(Q);
                QRmats.push_back(QRmat);
                svals.push_back(s);
                length += s.length();
            }

            TIT--;
            TITR--;
            TITgram--;

        } else if (!sorted) {
            TIT     = save;
            TITR    = saver;
            TITgram = savegram;
            sorted  = true;
            counter = 0;

            Vec s_all(length);
            int off = 1;
            for (std::size_t i=0; i<svals.size(); ++i) {
                s_all(_(off,off+svals[i].length()-1)) = svals[i];
                off += svals[i].length();
            }


            IntVec rho(s_all.length());
            flens::sort(s_all, rho);
            double delta = 0.;
            double eps2  = eps(depth-level+1);
            off          = 0;
            for (; off<s_all.length(); ++off) {
                delta += s_all(off+1);
                if (delta>eps2) break;
            }

            IntVec dels(svals.size());
            for (int i=1; i<=off; ++i) {
                int right = 0;
                for (std::size_t j=0; j<svals.size(); ++j) {
                    right += svals[j].length();
                    if (rho(i)<=right) {
                        ++dels(j+1);
                        break;
                    }
                }
            }

            ranks.resize(svals.size());
            for (std::size_t i=0; i<svals.size(); ++i) {
                int rank = svals[i].length() - dels(i+1);
                if (rank==0) rank = 1;
                ranks[i] = rank;
            }

            --level;
        } else {
            if (node->isLeaf()) {

                Matrix tmp1, tmp2, tmp3, tmp4;
                Matrix& Uref = Urefs[counter];
                Matrix& U    = Us[counter];
                Matrix& Q    = Qs[counter];
                int rank     = ranks[counter];
                int m        = node->getContent()->getUorB().numRows();
                int n        = node->getContent()->getUorB().numCols();;

                if ( m < n){
                    tmp2 = U(_,_(1,std::min(U.numCols(),rank)));
                    tmp3 = tmp2;
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    Uref(_,_(1,m)).upper(),tmp2);
                    flens::blas::mm(cxxblas::Trans,
                                    cxxblas::NoTrans, 1.0,
                                    Uref(_,_(m+1,n)),tmp3,0.0,tmp4);
                    tmp1.resize(0, 0);
                    tmp1 = Matrix(tmp2.numRows()+tmp4.numRows(),tmp2.numCols());
                    tmp1(_(1,tmp2.numRows()),_) = tmp2;
                    tmp1(_(tmp2.numRows()+1,tmp1.numRows()),_) = tmp4;
                } else {
                    tmp1 = U
                           (_(1,n),_(1,std::min(U.numCols(),
                                                rank)));
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    Uref(_(1,n),_).upper(),tmp1);

                }

                TITR.getNode()->getContent()->setUorB(tmp1);

                tmp1.resize(0, 0);
                flens::blas::mm(cxxblas::NoTrans,
                                cxxblas::NoTrans, 1.0,
                                Q,
                                U
                                (_,_(1,std::min(U.numCols(),rank))),0.0,tmp1);
                node->getContent()->setUorB(tmp1);
                node->getContent()->setNumRows(tmp1.numCols());
            } else if (node->isRoot()) {
                Matrix &leftMat  = TITR.getNode()->
                                   getfirstChild()->getContent()->getUorB();
                Matrix &rightMat = TITR.getNode()->
                                   getlastChild()->getContent()->getUorB();
                Matrix tmp1;

                //numel sollte hier immer 1 sein. Ansonsten genau das selbe wie weiter unten...

                flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,rightMat,leftMat,0.0,tmp1);

                //einsortiert werden muss hier nicht, da wir nicht mit QR weiterverarbeiten ist dies schon die Blockstruktur!

                node->getContent()->setUorB(tmp1);
                node->getContent()->setNumRows(1);
                node->getContent()->setLeftChildNumRows(tmp1.numCols());
                node->getContent()->setRightChildNumRows(tmp1.numRows());

            } else {
                int numel        = node->getContent()->getNumRows();
                Matrix &leftMat  = TITR.getNode()->getfirstChild()
                                   ->getContent()->getUorB();
                Matrix &rightMat = TITR.getNode()->getlastChild()
                                   ->getContent()->getUorB();


                Matrix tmp1, tmp2, tmp3, tmp4;
                Matrix& QRmat = QRmats[counter];
                Matrix& U     = Us[counter];
                Matrix& Q     = Qs[counter];
                int rank      = ranks[counter];
                int m         = leftMat.numCols()*rightMat.numCols();
                int n         = numel;

                if ( m < n){
                    tmp2 = U(_,_(1,std::min(U.numCols(),rank)));
                    tmp3 = tmp2;
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    QRmat(_,_(1,m)).upper(),tmp2);
                    flens::blas::mm(cxxblas::Trans,
                                    cxxblas::NoTrans, 1.0,
                                    QRmat(_,_(m+1,n)),tmp3,0.0,tmp4);
                    tmp1 = Matrix(tmp2.numRows()+tmp4.numRows(),tmp2.numCols());
                    tmp1(_(1,tmp2.numRows()),_) = tmp2;
                    tmp1(_(tmp2.numRows()+1,tmp1.numRows()),_) = tmp4;
                } else {
                    tmp1 = U(_(1,n),_(1,std::min(U.numCols(),rank)));
                    flens::blas::mm(cxxblas::Left,
                                    cxxblas::Trans, 1.0,
                                    QRmat(_(1,n),_).upper(),tmp1);

                }

                TITR.getNode()->getContent()->setUorB(tmp1);

                //Das hier ist Q_t*S_t, die muss nun zurück ins Blockformat!!!
                tmp1.resize(0, 0);
                flens::blas::mm(cxxblas::NoTrans,
                                cxxblas::NoTrans, 1.0,
                                Q,
                                U(_,_(1,std::min(U.numCols(),rank))),0.0,tmp1);
                int newnumel   = tmp1.numCols();
                int newlcnumel = node->getfirstChild()->getContent()->getNumRows();
                int newrcnumel = node->getlastChild()->getContent()->getNumRows();
                Matrix newB(newnumel*newrcnumel,newlcnumel);

                for(int i = 1; i <= newnumel; ++i){
                    for(int j = 1; j <= newlcnumel; ++j){
                        newB(_((i-1)*newrcnumel+1,i*newrcnumel),j) = tmp1(_((j-1)*newrcnumel + 1,j*newrcnumel) ,i);
                    }
                }

                node->getContent()->setUorB(newB);
                node->getContent()->setNumRows(newnumel);
                node->getContent()->setLeftChildNumRows(newlcnumel);
                node->getContent()->setRightChildNumRows(newrcnumel);
            }

            TIT--;
            TITR--;
            TITgram--;
            ++counter;
            newlevel = true;
        }
    }
}


template <typename T>
HTuckerTree<T> truncate(const HTuckerTree<T> & tree, const int rank){
    using flens::_;

	HTuckerTree<T> treecopy = tree;
	HTuckerTree<T> Stree(tree);
	treecopy.orthogonalize();
	HTuckerTree<T> gram = gramians_orthogonal(treecopy);
	flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > U(1,1);
	flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > Vt(1,1);
	flens::DenseVector<flens::Array<T> > s;
	GeneralTreeNode<HTuckerTreeNode<T> > *node;
	
    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TITgram = gram.getGeneralTree().begin();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITs= Stree.getGeneralTree().begin();
    for(; TITgram <= gram.getGeneralTree().end(); TITgram ++,TITs++){
		flens::svd(TITgram.getNode()->getContent()->getUorB(),s,U,Vt);
		TITs.getNode()->getContent()->setUorB(U(_,_(1,std::min(U.numRows(),
                                                      (FLENS_DEFAULT_INDEXTYPE)rank))));
	}
    }

    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = treecopy.getGeneralTree().begin();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITs = Stree.getGeneralTree().begin();
	for(; TIT <= treecopy.getGeneralTree().end(); TIT++,TITs++){
		node = TIT.getNode();
		
		if(node->isLeaf()){
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp(1,1);
			flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,node->getContent()->getUorB(),TITs.getNode()->getContent()->getUorB(),0.0,tmp);
			node->getContent()->setUorB(tmp);
			node->getContent()->setNumRows(tmp.numCols());
		} else {
			int rcnumel = node->getContent()->getRightChildNumRows();
			int lcnumel = node->getContent()->getLeftChildNumRows();
			
			int numel = node->getContent()->getNumRows();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Sref = TITs.getNode()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Slref = TITs.getNode()->getfirstChild()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Srref = TITs.getNode()->getlastChild()->getContent()->getUorB();
			int newnumel = Sref.numCols();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Bref = node->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp(rcnumel*newnumel,lcnumel);
			
			for(int l = 1; l <= newnumel; ++l){
				for(int j = 1; j <= rcnumel; ++j){
					for(int k = 1; k <= lcnumel; ++k){
						T sum = 0.0;
						for(int i = 1; i <= numel; ++i){
							sum += Bref((i-1)*rcnumel + j,k)*Sref(i,l);
						}
						tmp((l-1)*rcnumel + j,k) = sum;
					}
				}
			}
			//hier haben wir nun als B_t*S_t berechnet. Fehlt noch (S_tl^T x S_tr^T)*B_t*S_t
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > newB(newnumel*Srref.numCols(),Slref.numCols());
			for(int i = 1; i <= newnumel; ++i){
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp1(1,1);
				flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,tmp(_((i-1)*rcnumel+1,i*rcnumel),_),Slref,0.0,tmp1);
				flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp2(Srref.numCols(),tmp1.numRows());
				flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,Srref,tmp1,0.0,tmp2);
				newB(_((i-1)*Srref.numCols()+1,i*Srref.numCols()),_) = tmp2;
			}
			node->getContent()->setUorB(newB);
			node->getContent()->setNumRows(newnumel);
			node->getContent()->setLeftChildNumRows(Slref.numCols());
			node->getContent()->setRightChildNumRows(Srref.numCols());
		
		}
		
		
	}
    }

	return treecopy;
}


template <typename T>
HTuckerTree<T> truncate_nonorthogonal(const HTuckerTree<T> & tree, const int rank){
    using flens::_;

	HTuckerTree<T> treecopy = tree;
	
	
	HTuckerTree<T> gram = gramians_nonorthogonal(treecopy);
	HTuckerTree<T> Rtree(tree);

	flens::DenseVector<flens::Array<T> > s;
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > U, Vt;
	
	GeneralTreeNode<HTuckerTreeNode<T> > *node;
	
    {
    GeneralTreeIterator<HTuckerTreeNode<T> > TIT = treecopy.getGeneralTree().end();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITR = Rtree.getGeneralTree().end();
    GeneralTreeIterator<HTuckerTreeNode<T> > TITgram = gram.getGeneralTree().end();
    for(;TIT >= treecopy.getGeneralTree().begin(); TIT--,TITR--,TITgram--){
		node = TIT.getNode();
		if(node->isLeaf()){
			//compute QR-Decomposition
			//hier ist die Kopie notwendig, weil qrf die Matrix verändert!
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Uref = node->getContent()->getUorB();
			int m = Uref.numRows();
			int n = Uref.numCols();
			flens::DenseVector<flens::Array<T> > tau;
			qrf(Uref,tau);
			
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Q(m,std::min(m,n));
			Q(_(1,m),_(1,std::min(n,m))) = Uref(_,_(1,std::min(m,n)));
			orgqr(Q,tau);


			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Gref = TITgram.getNode()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp1,tmp2,tmp3,tmp4;

            tmp1 = Gref;
			//hier berechnen wir R_t*G_t*R_t^T   => R_t entsprciht U_ref, wenn m > n dann ist die Matrix höher als Breit und enthält nur 0en, welche rausgeschmissen werden können
			if (m < n){
				//dieser Fall ist komplizierter, weil die Matrix R zerlegt werden muss..
				tmp2 = Gref(_,_(1,m));
				tmp3 = Gref(_,_(m+1,n));
				flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,Uref(_,_(1,m)).upper(),tmp2);
                tmp1.resize(0, 0);
				flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0,tmp3,Uref(_,_(m+1,n)),0.0,tmp1);
                tmp3.resize(0, 0);
				tmp3 = tmp1 + tmp2;
                tmp1.resize(0, 0);
                tmp2.resize(0, 0);
				tmp1 = tmp3(_(1,m),_);
				tmp2 = tmp3(_(m+1,n),_);
				flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,Uref(_,_(1,m)).upper(),tmp1);
                tmp3.resize(0, 0);
				flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Uref(_,_(m+1,n)),tmp2,0.0,tmp3);
                tmp2.resize(0, 0);
				tmp2 = tmp1 + tmp3;
                tmp1.resize(0, 0);
				tmp1 = tmp2;
			} else {
				//d.h. in diesem Fall ist R in n x n und Q in m x n 
				flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,Uref(_(1,n),_).upper(),tmp1);
				flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,Uref(_(1,n),_).upper(),tmp1);
			}

			//flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,Uref(_(1,std::min(m,n)),_).upper(),tmp1);
			//flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,Uref(_(1,std::min(m,n)),_).upper(),tmp1);

			//flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0, Gref,Uref(_(1,std::min(m,n)),_).upper(),0.0,tmp1);
			//flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Uref(_(1,std::min(m,n)),_).upper(),tmp1,0.0,tmp2);

			flens::svd(tmp1,s,U,Vt);
			

			//Die S_t^T*  R_t werden hier transponiert gespeicher! => also R_t^T*S_t
			if ( m < n){
                tmp2.resize(0, 0);
                tmp3.resize(0, 0);
				tmp2 = U(_,_(1,std::min(U.numCols(),rank)));
				tmp3 = tmp2;
				flens::blas::mm(cxxblas::Left,cxxblas::Trans,1.0,Uref(_,_(1,m)).upper(),tmp2);
				flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,Uref(_,_(m+1,n)),tmp3,0.0,tmp4);
                tmp1.resize(0, 0);
				tmp1 = flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> >(tmp2.numRows()+tmp4.numRows(),tmp2.numCols());
				tmp1(_(1,tmp2.numRows()),_) = tmp2;
				tmp1(_(tmp2.numRows()+1,tmp1.numRows()),_) = tmp4;
			} else {
                tmp1.resize(0, 0);
				tmp1 = U(_(1,n),_(1,std::min(U.numCols(),rank)));
				flens::blas::mm(cxxblas::Left,cxxblas::Trans,1.0,Uref(_(1,n),_).upper(),tmp1);

			}
			
			
			//flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,U(_,_(1,std::min(U.numCols(),rank))),Uref(_(1,std::min(m,n)),_).upper(),0.0,tmp1);
			TITR.getNode()->getContent()->setUorB(tmp1);

            tmp1.resize(0, 0);
			flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Q,U(_,_(1,std::min(U.numCols(),rank))),0.0,tmp1);
			node->getContent()->setUorB(tmp1);
			node->getContent()->setNumRows(tmp1.numCols());
		} else if(node->isRoot()){
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &leftMat = TITR.getNode()->getfirstChild()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &rightMat = TITR.getNode()->getlastChild()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Bref = node->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp1, tmp2;

			//numel sollte hier immer 1 sein. Ansonsten genau das selbe wie weiter unten...
			

			flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Bref,leftMat,0.0,tmp1);
			flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,rightMat,tmp1,0.0,tmp2);
			
			//einsortiert werden muss hier nicht, da wir nicht mit QR weiterverarbeiten ist dies schon die Blockstruktur!

			node->getContent()->setUorB(tmp2);
			node->getContent()->setNumRows(1);
			node->getContent()->setLeftChildNumRows(tmp2.numCols());
			node->getContent()->setRightChildNumRows(tmp2.numRows());
			
		} else {
			//compute (Stl^T Rtl x Str^T Rtr) B_t
			int numel = node->getContent()->getNumRows();
			int rcnumel = node->getContent()->getRightChildNumRows();
			
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &leftMat = TITR.getNode()->getfirstChild()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &rightMat = TITR.getNode()->getlastChild()->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Bref = node->getContent()->getUorB();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > QRmat(leftMat.numCols()*rightMat.numCols(),numel);
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tmp1,tmp2,tmp3,tmp4;

			int lr = leftMat.numCols(); // wegen leftMat^T
			int rr = rightMat.numCols();
		

			for(int i = 1; i <= numel; ++i){
				flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Bref(_((i-1)*rcnumel + 1,i*rcnumel),_),leftMat,0.0,tmp1);
				flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,rightMat,tmp1,0.0,tmp2);
				//jetzt noch in die neue Matrix einsortieren...
				for(int j = 1; j <= lr; ++j){
					QRmat(_((j-1)*rr+1,j*rr),i) = tmp2(_,j);
				}
			}
			// jetzt die QR Zerlegung anwenden:
			int m = QRmat.numRows();
			int n = QRmat.numCols();
			flens::DenseVector<flens::Array<T> > tau;
			qrf(QRmat,tau);
			
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Q(m,std::min(m,n));
			Q(_(1,m),_) = QRmat(_,_(1,std::min(m,n)));
			orgqr(Q,tau);

			// Dann wie oben R_t G_t R_t^T berechnen
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &Gref = TITgram.getNode()->getContent()->getUorB();

            tmp1.resize(0, 0);
            tmp1 = Gref;
			if (m < n){
				//dieser Fall ist komplizierter, weil die Matrix R zerlegt werden muss..
                tmp2.resize(0, 0);
                tmp3.resize(0, 0);
				tmp2 = Gref(_,_(1,m));
				tmp3 = Gref(_,_(m+1,n));
				flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,QRmat(_,_(1,m)).upper(),tmp2);
                tmp1.resize(0, 0);
				flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0,tmp3,QRmat(_,_(m+1,n)),0.0,tmp1);
                tmp3.resize(0, 0);
				tmp3 = tmp1 + tmp2;
                tmp1.resize(0, 0);
				tmp1 = tmp3(_(1,m),_);
                tmp2.resize(0, 0);
				tmp2 = tmp3(_(m+1,n),_);
				flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,QRmat(_,_(1,m)).upper(),tmp1);
                tmp3.resize(0, 0);
				flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,QRmat(_,_(m+1,n)),tmp2,0.0,tmp3);
                tmp2.resize(0, 0);
				tmp2 = tmp1 + tmp3;
                tmp1.resize(0, 0);
				tmp1 = tmp2;
			} else {
				//d.h. in diesem Fall ist R in n x n und Q in m x n 
				flens::blas::mm(cxxblas::Right,cxxblas::Trans,1.0,QRmat(_(1,n),_).upper(),tmp1);
				flens::blas::mm(cxxblas::Left,cxxblas::NoTrans,1.0,QRmat(_(1,n),_).upper(),tmp1);
			}
			//flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.0, Gref,QRmat(_(1,std::min(m,n)),_).upper(),0.0,tmp1);
			//flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,QRmat(_(1,std::min(m,n)),_).upper(),tmp1,0.0,tmp2);

			//Svd wie oben

			flens::svd(tmp1,s,U,Vt);
            std::cout << " Normal trunc, s sizes in inner = " << s.length()
                      << std::endl;
			//tmp1 = U(_,_(1,std::min(U.numCols(),rank)));
			if ( m < n){
                tmp2.resize(0, 0);
				tmp2 = U(_,_(1,std::min(U.numCols(),rank)));
                tmp3.resize(0, 0);
				tmp3 = tmp2;
				flens::blas::mm(cxxblas::Left,cxxblas::Trans,1.0,QRmat(_,_(1,m)).upper(),tmp2);
				flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,QRmat(_,_(m+1,n)),tmp3,0.0,tmp4);
                tmp1.resize(0, 0);
				tmp1 = flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> >(tmp2.numRows()+tmp4.numRows(),tmp2.numCols());
				tmp1(_(1,tmp2.numRows()),_) = tmp2;
				tmp1(_(tmp2.numRows()+1,tmp1.numRows()),_) = tmp4;
			} else {
                tmp1.resize(0, 0);
				tmp1 = U(_(1,n),_(1,std::min(U.numCols(),rank)));
				flens::blas::mm(cxxblas::Left,cxxblas::Trans,1.0,QRmat(_(1,n),_).upper(),tmp1);

			}
			//flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.0,U(_,_(1,std::min(U.numCols(),rank))),Uref(_(1,std::min(m,n)),_).upper(),0.0,tmp1);
			TITR.getNode()->getContent()->setUorB(tmp1);

			//Das hier ist Q_t*S_t, die muss nun zurück ins Blockformat!!!
            tmp1.resize(0, 0);
			flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.0,Q,U(_,_(1,std::min(U.numCols(),rank))),0.0,tmp1);
			int newnumel = tmp1.numCols();
			int newlcnumel = node->getfirstChild()->getContent()->getNumRows();
			int newrcnumel = node->getlastChild()->getContent()->getNumRows();
			flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > newB(newnumel*newrcnumel,newlcnumel);

			for(int i = 1; i <= newnumel; ++i){
				for(int j = 1; j <= newlcnumel; ++j){
					newB(_((i-1)*newrcnumel+1,i*newrcnumel),j) = tmp1(_((j-1)*newrcnumel + 1,j*newrcnumel) ,i);
				}
			}

			node->getContent()->setUorB(newB);
			node->getContent()->setNumRows(newnumel);
			node->getContent()->setLeftChildNumRows(newlcnumel);
			node->getContent()->setRightChildNumRows(newrcnumel);
		}

		
	
	}
    }

    return treecopy;
}


template <typename T>
HTuckerTree2Tensor<T>::HTuckerTree2Tensor():isset(false),tree(NULL){}

template <typename T>
HTuckerTree2Tensor<T>::HTuckerTree2Tensor(const int min, const int max, const int d):isset(false){}

template <typename T>
HTuckerTree2Tensor<T>::HTuckerTree2Tensor(const DimensionIndex & _min, const DimensionIndex &_max):isset(false){}

template <typename T>
void
HTuckerTree2Tensor<T>::setTree(const HTuckerTree<T> &_tree){
	tree = &_tree;
	isset = true;
}

template <typename T>
T
HTuckerTree2Tensor<T>::operator() (const DimensionIndex &vals) const{
	if(isset){
		return tree->evaluate(vals);
	}
	return 0;
}


template <typename T>
bool
HTuckerTree2Tensor<T>::vecEval()const{
	return true;
}


template <typename T>
void  
HTuckerTree2Tensor<T>::vec(const DimensionIndex & vals, const int dim, flens::DenseVector<flens::Array<type> > & vec) const{
	if(isset){
		vec = tree->vec_evaluate(vals,dim);
	} else {
		vec = flens::DenseVector<flens::Array<T> >(1);
	}
}







	
} // namespace htucker
