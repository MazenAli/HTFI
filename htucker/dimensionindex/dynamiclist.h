#ifndef HTUCKER_DIMENSIONINDEX_DYNAMICLIST_H
#define HTUCKER_DIMENSIONINDEX_DYNAMICLIST_H 1

#include <iostream>

namespace htucker{

template <typename T>
class dynamiclistelement{

	public:

	T value;
	dynamiclistelement * next;
	dynamiclistelement * previous;
	
	dynamiclistelement(const T &val): value(val),next(NULL),previous(NULL){};

};

template <typename T>
class dynamiclist{

	int count;
	

	public:

	dynamiclistelement<T> * root;
	dynamiclistelement<T> * last;

	dynamiclist();

	void
	append(const T &val);
	
	T
	operator[](const int i) const;

	int
	length() const;
};

#include <htucker/dimensionindex/dynamiclist.tcc>


} // namespace htucker

#endif // HTUCKER_DIMENSIONINDEX_DYNAMICLIST_H
