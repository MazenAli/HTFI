#ifndef HTUCKER_MATRIXTENSORS_HTUCKERCLOSURE_H
#define HTUCKER_MATRIXTENSORS_HTUCKERCLOSURE_H 1

#include <htucker/htuckertree/htuckerelement.h>
#include <htucker/matrixtensors/matrixtensor.h>
namespace htucker{




template <typename op,class A, class B>
class HTuckerClosure{
	
	A   left;
	B   right;

	public:

	int d;

	HTuckerClosure(const A & _left, const B & _right, const int _d);

	A & getLeft() const;

	B & getRight() const;

	int dim() const;



};




}

#include <htucker/matrixtensors/htuckerclosure.tcc>


#endif //HTUCKER_MATRIXTENSORS_HTUCKERCLOSURE_H
