#ifndef HTUCKER_MATRIXTENSORS_HTUCKERTREEPART_H
#define HTUCKER_MATRIXTENSORS_HTUCKERTREEPART_H 1


namespace htucker{

	template <typename T>
	class HTuckerTreePart{
		public:
		
		HTuckerTree<T> & httree;
		int mindim, maxdim;


		HTuckerTreePart(const HTuckerTree<T> & tree, const int _min, const int _max);
		
	};


}; //namespace htucker

#include <htucker/matrixtensors/htuckertreepart.tcc>


#endif //HTUCKER_MATRIXTENSORS_HTUCKERTREEPART_H
