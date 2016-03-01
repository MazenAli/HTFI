namespace htucker{

	IdentityTensor::IdentityTensor(const DimensionIndex &_maxvals):minvals(1,_maxvals.length()),maxvals(_maxvals){}

	IdentityTensor::IdentityTensor(const DimensionIndex & _minvals, const DimensionIndex & _maxvals):minvals(_minvals),maxvals(_maxvals){}

	int 
	IdentityTensor::dim() const{
		return maxvals.length();
	}
} //namespace htucker
