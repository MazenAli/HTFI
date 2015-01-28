#ifndef HTUCKER_DIMENSIONINDEX_DENSEVECTORLIST_H
#define  HTUCKER_DIMENSIONINDEX_DENSEVECTORLIST_H 1

#include <flens/flens.h>

namespace htucker{

template <typename T>
class DenseVectorListElement{
public: 
	flens::DenseVector<flens::Array<T> > elem;
	DenseVectorListElement * next;
    
	DenseVectorListElement(const flens::DenseVector<flens::Array<T> > &vector):elem(vector), next(NULL){};
};

template <typename T>
class DenseVectorList{
private:
	DenseVectorListElement<T> * node;
	DenseVectorListElement<T> * iterator;
	int len;
	int pos;
	DenseVectorListElement<T> *lastnode;
public:
	
	//constructors
	DenseVectorList();
	
	DenseVectorList(const flens::DenseVector<flens::Array<T> > &vector);

	~DenseVectorList();

	//appends a new DenseVector at the End of the List
	void 
		add(const flens::DenseVector<flens::Array<T> > &vector);

	//removes Element i from the list
	void 
	remove(const int i);

	//empty list
	void 
	empty();

	// get len
	int 
	length() const;
	
	// () Operator return the DenseVector at position i of the List
	flens::DenseVector<flens::Array<T> > * 
	operator()(const int i);
	
};
} // namespace htucker
#include <htucker/dimensionindex/densevectorlist.tcc>

#endif // HTUCKER_DIMENSIONINDEX_DENSEVECTORLIST_H
