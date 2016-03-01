Matrix2Tensor::Matrix2Tensor(const int _min, const int _max, const int _d):d(_d),minval(DimensionIndex(_min,_d)),maxval(DimensionIndex(_max,_d)),matrixset(false){}

Matrix2Tensor::Matrix2Tensor(const DimensionIndex &_min, const DimensionIndex &_max):d(_min.length()),minval(_min),maxval(_max),matrixset(false){}

Matrix2Tensor::Matrix2Tensor(const Matrix2Tensor &copy):d(copy.dim()),minval(copy.getmin()),maxval(copy.getmax()),matrixset(copy.getmatrixset()),M(&copy.getMatrix()),conv(copy.getconverter()){
}

void 
Matrix2Tensor::setMatrix(flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > &_M){
	assert(_M.numRows() == _M.numCols());
	//M = flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> >(_M.numRows(),_M.numCols());
	M = &_M;
	DimensionIndex sqrtdim(minval.length());
	for(int i = sqrtdim.length()-1; i>= 0; --i){
		sqrtdim[i] = sqrt(maxval[i] - minval[i] + 1.0) + minval[i] - 1;
	}
	
	conv = m2vindexconverter(minval,sqrtdim);
	matrixset = true;
}

void 
Matrix2Tensor::setMatrix(flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > &_M, const DimensionIndex &_minidxrow, const DimensionIndex &_maxidxrow){
	//M = flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> >(_M.numRows(),_M.numCols());
	M = &_M;
	DimensionIndex coldim(minval.length());
	for(int i = coldim.length()-1; i>= 0; --i){
		coldim[i] = (maxval[i]-minval[i] + 1)/(_maxidxrow[i] - _minidxrow[i] + 1) + minval[i] - 1;
	}
	
	conv = m2vindexconverter(minval,coldim);
	conv.setRowIndices(_minidxrow,_maxidxrow);
	matrixset = true;
}

void
Matrix2Tensor::vec() const{
}



int 
Matrix2Tensor::dim() const{
	return d;
}

const DimensionIndex &
Matrix2Tensor::getmin() const{
	return minval;
}

const DimensionIndex &
Matrix2Tensor::getmax() const{
	return maxval;
}

flens::GeMatrix<flens::FullStorage<double,cxxblas::ColMajor> > &
Matrix2Tensor::getMatrix() const{
	return *M;
}

bool 
Matrix2Tensor::getmatrixset() const{
	return matrixset;
}

const m2vindexconverter &
Matrix2Tensor::getconverter() const{
	return conv;
}


double 
Matrix2Tensor::operator()(const DimensionIndex &vals) const{
	if(matrixset){
		int r,c;
		conv.getMatrixIndices(vals,r,c);
		return M->operator()(r,c);
	}
	return 0.0;
}

bool
Matrix2Tensor::vecEval() const{
	return false;
}

Matrix2Tensor &
Matrix2Tensor::operator=(const Matrix2Tensor &copy){
	this->d = copy.dim();
	minval = copy.getmin();
	maxval = copy.getmax();
	matrixset = copy.getmatrixset();
	M = &(copy.getMatrix());
	conv = copy.getconverter();
    return *this;
}
