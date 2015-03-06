namespace htucker{

	IdentityTensor::IdentityTensor(const DimensionIndex &_maxvals):maxvals(_maxvals),minvals(1,_maxvals.length()){};

	IdentityTensor::IdentityTensor(const DimensionIndex & _minvals, const DimensionIndex & _maxvals):maxvals(_maxvals),minvals(_minvals){};

	int 
	IdentityTensor::dim() const{
		return maxvals.length();
	}
} //namespace htucker
