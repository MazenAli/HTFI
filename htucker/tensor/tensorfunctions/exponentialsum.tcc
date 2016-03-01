exponentialSum::exponentialSum(const int _min, const int _max, const int _d):d(_d),minval(DimensionIndex(_min,_d)),maxval(DimensionIndex(_max,_d)){}

exponentialSum::exponentialSum(const DimensionIndex &_min, const DimensionIndex &_max):d(_min.length()),minval(_min),maxval(_max){
	assert(_min.length()==_max.length());
}

exponentialSum::exponentialSum(const exponentialSum & copy):d(copy.dim()),minval(copy.getmin()),maxval(copy.getmax()){}

int 
exponentialSum::dim() const{
	return d;
}

const DimensionIndex & 
exponentialSum::getmin() const{
	return minval;
}

const DimensionIndex &
exponentialSum::getmax() const{
	return maxval;
}

double 
exponentialSum::operator() (const DimensionIndex &vals) const{
	assert(vals.length() == d);
	double sum = 0.0;
	double val;
	for(int i = 0; i< d; ++i){
		val = (vals[i]/(maxval[i] - minval[i] + 1.0));
		sum += val*val;
	}
	return exp(-sqrt(sum));

}

bool 
exponentialSum::vecEval() const {
	return false;
}

void
exponentialSum::vec() const{
}


exponentialSum & 
exponentialSum::operator=(const exponentialSum &copy){
	this->d = copy.dim();
	this->minval = copy.getmin();
	this->maxval = copy.getmax();
    return *this;
}
