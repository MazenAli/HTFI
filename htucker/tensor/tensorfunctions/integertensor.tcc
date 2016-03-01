IntegerTensor::IntegerTensor(const int _min, const int _max, const int _d):d(_d),minval(DimensionIndex(_min,_d)),maxval(DimensionIndex(_max,_d)),offset(0){
	assert(_min == 0);
}

IntegerTensor::IntegerTensor(const DimensionIndex &_min, const DimensionIndex &_max):d(_min.length()),minval(_min),maxval(_max),offset(0){
	assert(_min.length() == _max.length());
	for(int i = 0; i < _min.length(); ++i){
		assert(_max[i] == _max[0] && _min[i] == 0);
	}
	
}

IntegerTensor::IntegerTensor(const IntegerTensor &copy):d(copy.dim()),minval(copy.getmin()),maxval(copy.getmax()),offset(copy.getoffset()){}



double 
IntegerTensor::operator()(const DimensionIndex &vals) const{
	double res = 0;
	int j = 1;
	for(int i = vals.length() - 1; i >= 0; --i){
		res+=  j*(vals[i]);
		j*=(maxval[0]+1);
	}
	return res + offset;

}

bool IntegerTensor::vecEval() const {
	return true;
}

void 
IntegerTensor::setoffset(const int _offset){
	offset = _offset;
}


int 
IntegerTensor::getoffset() const{
	return offset;
}

const DimensionIndex &
IntegerTensor::getmax() const{
	return maxval;
}

const DimensionIndex &
IntegerTensor::getmin() const{
	return minval;
}

int 
IntegerTensor::dim() const{
	return d;
}

void
IntegerTensor::vec(const DimensionIndex & vals, const int dim, flens::DenseVector<flens::Array<double> > & vec) const{
	double res = 0;
	int j = 1;
	int jdim;
	for(int i = vals.length() - 1; i >= 0; --i){
		if(i == dim-1){
			jdim = j;
		} else {
			res+=  j*(vals[i]);
		}
		j*=(maxval[0]+1);
	}

	int fi = vec.firstIndex();
	int la = vec.lastIndex() - vec.firstIndex();
	for(int i = 0; i<= la; ++i){
		vec(i+fi) = res + jdim * i + offset;
	}
}

IntegerTensor&
IntegerTensor::operator=(const IntegerTensor &copy){
	this->d = copy.dim();
	this->minval = copy.getmin();
	this->maxval = copy.getmax();
	this->offset = copy.getoffset();
    return *this;
}
