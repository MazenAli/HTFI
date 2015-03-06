#ifndef HTUCKER_MATRIXTENSORS_IDENTITYTENSOR_H
#define HTUCKER_MATRIXTENSORS_IDENTITYTENSOR_H 1


#include <htucker/dimensionindex/dimensionindex.h>


namespace htucker{

	class IdentityTensor{
		public:
		
			
		DimensionIndex minvals,maxvals;

		IdentityTensor(const DimensionIndex & _maxvals);

		IdentityTensor(const DimensionIndex & _minvals, const DimensionIndex & _maxvals);

		int dim() const;
		
	};


} //namespace htucker

#include <htucker/matrixtensors/identitytensor.tcc>


#endif //HTUCKER_MATRIXTENSORS_IDENTITYTENSOR_H
