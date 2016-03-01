Vector2Tensor::Vector2Tensor(const int _min, const int _max, const int _d):d(_d),minval(DimensionIndex(_min,_d)),maxval(DimensionIndex(_max,_d)),vectorset(false){}

Vector2Tensor::Vector2Tensor(const DimensionIndex &_min, const DimensionIndex &_max):d(_min.length()),minval(_min),maxval(_max),vectorset(false){}

Vector2Tensor::Vector2Tensor(const Vector2Tensor & copy):d(copy.dim()),minval(copy.getmin()),maxval(copy.getmax()),vectorset(copy.VectorSet())
{
    // Correction: actually create a copy // Mazen
    flens::DenseVector<flens::Array<double> > _v(copy.getVector());
    v = &_v;
}

void 
Vector2Tensor::setVector(flens::DenseVector<flens::Array<double > > &_v){
	vectorset = true;
	v = &_v;
}

double 
Vector2Tensor::operator()(const DimensionIndex &vals) const{
	int pos = vals.computeBEvalue(minval,maxval);
	return v->operator()(pos);
}

bool
Vector2Tensor::vecEval() const{
	return false;	
}

void  
Vector2Tensor::vec() const{
}

bool 
Vector2Tensor::VectorSet() const{
	return vectorset;
}

int 
Vector2Tensor::dim() const{
	return d;
}

const DimensionIndex & 
Vector2Tensor::getmin() const{
	return minval;
}

const DimensionIndex &
Vector2Tensor::getmax() const{
	return maxval;
}

const flens::DenseVector<flens::Array<double > > &
Vector2Tensor::getVector() const{
	return *v;
}
	
Vector2Tensor&
Vector2Tensor::operator=(const Vector2Tensor &copy){
	d = copy.dim();
	minval = copy.getmin();
	maxval = copy.getmax();

    //Correction: actually create a copy // Mazen
    flens::DenseVector<flens::Array<double> > _v(copy.getVector());
    v = &_v;
	vectorset = copy.VectorSet();
	return *this;
}
