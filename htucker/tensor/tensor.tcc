//constructors

template <typename TensorFunction>
tensor<TensorFunction>::tensor():dimension(1),tf(0,0,1){}

template <typename TensorFunction>
tensor<TensorFunction>::tensor(const int _min, const int _max, const int d):dimension(d),tf(_min,_max,d){
	minval = DimensionIndex(d);
	maxval = DimensionIndex(d);
	for(int i = 0; i < minval.length(); ++i){
		minval[i] = _min;
		maxval[i] = _max;
	}
}


template <typename TensorFunction>
tensor<TensorFunction>::tensor(const DimensionIndex &_min, const DimensionIndex &_max):dimension(_min.length()),minval(_min),maxval(_max),tf(_min,_max){
	assert(minval.length() == dimension && maxval.length() == dimension);
}

template <typename TensorFunction>
tensor<TensorFunction>::tensor(const tensor<TensorFunction> &copy):dimension(copy.dim()),minval(copy.getmin()),maxval(copy.getmax()),tf(copy.getTensorFunction()){
}



template <typename TensorFunction>
typename TensorFunction::type
tensor<TensorFunction>::LinfNorm(const int n) const{
	typename TensorFunction::type maxi = -1;
	DimensionIndex eval(dimension);
	typename TensorFunction::type evaluate;
	for(int i = 1; i<= n ; ++i){
		eval.setRandom(minval,maxval);
		evaluate = abs(tf(eval));
		if(evaluate > maxi){
			maxi = evaluate;
		}
	}
	return maxi;
}

template <typename TensorFunction>
typename TensorFunction::type
tensor<TensorFunction>::operator()(const DimensionIndex &vals) const{
	return tf(vals);
}

template <typename TensorFunction>
int
tensor<TensorFunction>::dim() const{
	return dimension;
}

template <typename TensorFunction>
DimensionIndex & 
tensor<TensorFunction>::getmax() const{
	return maxval;
}

template <typename TensorFunction>
DimensionIndex &
tensor<TensorFunction>::getmin() const{
	return minval;
}

template <typename TensorFunction>
bool
tensor<TensorFunction>::vecEval() const{
	return tf.vecEval();
}

template <typename TensorFunction>
TensorFunction& 
tensor<TensorFunction>::getTensorFunction() const{
	return tf;
}


template <typename TensorFunction>
flens::DenseVector<flens::Array<typename TensorFunction::type> >
tensor<TensorFunction>::vec(const DimensionIndex & vals, const int dim) const{
	flens::DenseVector<flens::Array<typename TensorFunction::type> > v(maxval[dim-1]+1);
	tf.vec(vals,dim,v);
	return v;
}

template <typename TensorFunction>
tensor<TensorFunction> &
tensor<TensorFunction>::operator=(const tensor<TensorFunction> &copy){
	dimension = copy.dim();
	minval = copy.getmin();
	maxval = copy.getmax();
	tf = copy.getTensorFunction();
}
