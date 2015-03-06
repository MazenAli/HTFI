#ifndef HTUCKER_MATRIXTENSORS_VECTRORTENSOR_H
#define HTUCKER_MATRIXTENSORS_VECTORTENSOR_H 1

#include <htucker/htuckertree/htuckerelement.h>
#include <htucker/matrixtensors/identitytensor.h>

namespace htucker{

	template <typename T>
	class VectorTensor{
			
		flens::DenseVector<flens::Array<T> > * vectorlist;
		int d;

		public:

		VectorTensor(const int _d);

		~VectorTensor();

		void 
		setVector(const int _d,const flens::DenseVector<flens::Array<T> > & vec);

		flens::DenseVector<flens::Array<T> > & 
		getVector(const int _d) const;

		int dim() const;

		
	};

} //namespace htucker

#include <htucker/matrixtensors/vectortensor.tcc>


#endif // HTUCKER_MATRIXTENSORS_VECTORTENSOR_H
