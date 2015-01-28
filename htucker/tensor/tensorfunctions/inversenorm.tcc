inverseNorm::inverseNorm(const int _min, const int _max,const int _d):d(_d),minval(DimensionIndex(_min,_d)),maxval(DimensionIndex(_max,_d)){};

inverseNorm::inverseNorm(const DimensionIndex &_min, const DimensionIndex &_max):d(_min.length()),minval(_min),maxval(_max){
	assert(d == _min.length() && d == _max.length());
};

inverseNorm::inverseNorm(const inverseNorm & copy):d(copy.dim()),minval(copy.getmin()),maxval(copy.getmax()){}

double
inverseNorm::operator()(const DimensionIndex &vals) const{
	assert(vals.length() == d);
	double res = 0.0;
	for(int i = 0; i < vals.length(); ++i){
		res+= vals[i]*vals[i];
	}
	return 1.0/sqrt(res);
};

int 
inverseNorm::dim() const{
	return d;
}

const DimensionIndex &
inverseNorm::getmin() const{
	return minval;
}

const DimensionIndex &
inverseNorm::getmax() const{
	return maxval;
}

bool 
inverseNorm::vecEval() const{
	return false;
}

void  
inverseNorm::vec(const DimensionIndex & vals, const int dim, flens::DenseVector<flens::Array<type> > & vec) const{
	
}


inverseNorm&
inverseNorm::operator=(const inverseNorm &copy){
	this->d = copy.dim();
	this->minval = copy.getmin();
	this->maxval = copy.getmax();
    return *this;
}
