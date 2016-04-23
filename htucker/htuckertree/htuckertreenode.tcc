namespace htucker{

template <typename T>
HTuckerTreeNode<T>::HTuckerTreeNode():UorB(1,1){}

template <typename T>
HTuckerTreeNode<T>::HTuckerTreeNode(const DimensionIndex &_index):index(_index),UorB(1,1){}


template <typename T>
HTuckerTreeNode<T>::HTuckerTreeNode(const HTuckerTreeNode<T> &_copy):index(_copy.getIndex()),UorB(_copy.getUorB()),UorB_rcnumel(_copy.UorB_rcnumel),UorB_lcnumel(_copy.UorB_lcnumel),UorB_numel(_copy.UorB_numel){}

template <typename T>
HTuckerTreeNode<T> &
HTuckerTreeNode<T>::operator= (const HTuckerTreeNode<T> &_copy){
	this->index = _copy.getIndex();
	this->UorB = _copy.getUorB();
	this->UorB_numel = _copy.getNumRows();
	this->UorB_lcnumel = _copy.getLeftChildNumRows();
	this->UorB_rcnumel = _copy.getRightChildNumRows();
}

template <typename T>
const DimensionIndex &
HTuckerTreeNode<T>::getIndex() const{
	return index;
}


template <typename T>
DimensionIndex &
HTuckerTreeNode<T>::getIndex(){
	return index;
}


template <typename T>
const flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &
HTuckerTreeNode<T>::getUorB() const{
	return UorB;
}

template <typename T>
flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &
HTuckerTreeNode<T>::getUorB(){
	return UorB;
}


template <typename T>
void 
HTuckerTreeNode<T>::setUorB(const flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &_UorB){
	if (UorB.numRows()!=_UorB.numRows() || UorB.numCols()!=_UorB.numCols()) {
        UorB.resize(_UorB.numRows(), _UorB.numCols());
    }
    UorB = _UorB;
}


//template <typename T>
//flens::DenseVector<flens::Array<T> > & 
//HTuckerTreeNode<T>::getEvaluate() const{
//	return evaluate;
//};
//
//template <typename T>
//void 
//HTuckerTreeNode<T>::setEvaluate(const flens::DenseVector<flens::Array<T> > & _evaluate){
//	evaluate = _evaluate;
//};


template <typename T>
void
HTuckerTreeNode<T>::setIndex(const DimensionIndex & _ind){
	index = _ind;
}

template <typename T>
int 
HTuckerTreeNode<T>::getNumRows() const{
	return UorB_numel;
}
		
template <typename T>
int 
HTuckerTreeNode<T>::getLeftChildNumRows() const{
	return UorB_lcnumel;
}

template <typename T>
int 
HTuckerTreeNode<T>::getRightChildNumRows() const{
	return UorB_rcnumel;
}

template <typename T>
void 
HTuckerTreeNode<T>::setNumRows(const int numel){
	UorB_numel = numel;
}

template <typename T>
void 
HTuckerTreeNode<T>::setLeftChildNumRows(const int lcnumel){
	UorB_lcnumel = lcnumel;
}

template <typename T>
void 
HTuckerTreeNode<T>::setRightChildNumRows(const int rcnumel){
	UorB_rcnumel = rcnumel;
}

template <typename T>
std::ostream& operator<<(std::ostream& Stream, const HTuckerTreeNode<T> & B)
{
	Stream << B.getIndex();
	return Stream;
}

} // namespace htucker


