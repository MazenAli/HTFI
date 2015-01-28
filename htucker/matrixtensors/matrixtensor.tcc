namespace htucker{

template <typename mattype>
MatrixTensor<mattype>::MatrixTensor(const int _d):d(_d){
	this->matrixlist = (mattype *) malloc (d*sizeof(mattype));
	
	for(int i = 0; i < d;++i){
		mattype * bas = new (this->matrixlist + i) mattype();
	}
};

template <typename mattype>
MatrixTensor<mattype>::~MatrixTensor(){
	for(int i = d-1; i >= 0 ; i--){
		matrixlist[i].~mattype();
	}
	free(matrixlist);
};

template <typename mattype>
void 
MatrixTensor<mattype>::setMatrix(const int _d,const mattype & mat){
	assert(_d>= 1 && _d <= d);
	matrixlist[_d-1] = mat;
};


template <typename mattype>
mattype & 
MatrixTensor<mattype>::getMatrix(const int _d) const{
	return matrixlist[_d-1];
};

template <typename mattype>
int
MatrixTensor<mattype>::dim() const{
	return d;
};






}; //namespace htucker