#ifndef HTUCKER_MATRIXTENSORS_MATRIXTENSOR_H
#define HTUCKER_MATRIXTENSORS_MATRIXTENSOR_H 1

#include <htucker/htuckertree/htuckerelement.h>
#include <htucker/matrixtensors/identitytensor.h>

namespace htucker{

	template <typename mattype>
	class MatrixTensor{
			
		mattype * matrixlist;
		int d;

		public:

		MatrixTensor(const int _d);

		~MatrixTensor();

		void 
		setMatrix(const int _d,const mattype & mat);

		mattype & 
		getMatrix(const int _d) const;

		int dim() const;

		
	};

} //namespace htucker

#include <htucker/matrixtensors/matrixtensor.tcc>
#include <htucker/matrixtensors/htuckerclosure.h>
#include <htucker/matrixtensors/htuckertreepart.h>
#include <htucker/matrixtensors/vectortensor.h>
#include <htucker/matrixtensors/operators.h>
//#include <htucker/matrixtensors/matrixtensor_tests.h>
//#include <htucker/matrixtensors/htuckerclosure_tests.h>

#endif // HTUCKER_MATRIXTENSORS_MATRIXTENSOR_H