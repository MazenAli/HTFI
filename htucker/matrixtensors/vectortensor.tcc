namespace htucker{

template <typename T>
VectorTensor<T>::VectorTensor(const int _d):d(_d){
	this->vectorlist= (flens::DenseVector<flens::Array<T> > *) malloc (d*sizeof(flens::DenseVector<flens::Array<T> >));
	
	for(int i = 0; i < d;++i){
		flens::DenseVector<flens::Array<T> > * bas = new (this->vectorlist + i) flens::DenseVector<flens::Array<T> >();
	}
}

template <typename T>
VectorTensor<T>::~VectorTensor(){
	for(int i = d-1; i >= 0 ; i--){
		vectorlist[i].~DenseVector();
	}
	free(vectorlist);
}

template <typename T>
void 
VectorTensor<T>::setVector(const int _d,const flens::DenseVector<flens::Array<T> > & vec){
	assert(_d>= 1 && _d <= d);
	vectorlist[_d-1] = vec;
}


template <typename T>
flens::DenseVector<flens::Array<T> > & 
VectorTensor<T>::getVector(const int _d) const{
	return vectorlist[_d-1];
}

template <typename T>
int
VectorTensor<T>::dim() const{
	return d;
}


} //namespace htucker
