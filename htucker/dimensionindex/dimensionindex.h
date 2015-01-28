#ifndef HTUCKER_DIMENSIONINDEX_DIMENSIONINDEX_H
#define HTUCKER_DIMENSIONINDEX_DIMENSIONINDEX_H 1

#include <htucker/dimensionindex/dynamiclist.h>
#include <htucker/dimensionindex/densevectorlist.h>
#include <htucker/dimensionindex/dimensionindexiteratortypes.h>

namespace htucker{

class DimensionIndexList;

template <typename T>
class DimensionIndexIterator{}; 

class diinitializer;

class DimensionIndex{
	int * index;
	int len;
public:
	
	DimensionIndex();

	DimensionIndex(const int _d);

	DimensionIndex(const int * _index, const int d_);

	DimensionIndex(const int constval, const int _d);

	DimensionIndex(const DimensionIndex &copy);

	~DimensionIndex();

	//before calling setRandom, please set   srand ( time(NULL) );  !

	void setRandom(const int _min, const int _max);
	
	void setRandom(const DimensionIndex &_min, const DimensionIndex &_max);
	
	void setRandom(const int _pos, const int _min, const int _max);

	void setRandom(const DimensionIndex &_dimidx, const DimensionIndex &_min, const DimensionIndex &_max);

	void setRandom(const DimensionIndexList &_dimlist);

	void setRandom(const DimensionIndexList &_dimlist,const DimensionIndex &_activedims);

	void setValue(const int _value);

	void setValue(const int _pos, const int _value);

	void setValue(const DimensionIndex &_values, const DimensionIndex &_activedims);

	void setValueAsc();

	void setValueAsc(const int _min);

	bool equals(const DimensionIndex &_anotherindex, const DimensionIndex &_compareDims) const;

	bool equals(const DimensionIndexList &_dimlist, const DimensionIndex &_compareDims) const;

	int computeBEvalue(const int _min, const int _max) const;

	int computeBEvalue(const DimensionIndex & _min, const DimensionIndex & _max) const;

	int computeLEvalue(const int _min, const int _max) const;

	int computeLEvalue(const DimensionIndex & _min, const DimensionIndex &_max) const;
	
	DimensionIndex
	getComplement(const int dim) const;

	DimensionIndex
	join(const DimensionIndex & rhs) const;

	DimensionIndex
	intersect(const DimensionIndex &rhs) const;

	DimensionIndexIterator<Iterator1D> 
	getIterator(const int _dim, const int _min, const int _max) const;

	DimensionIndexIterator<Iterator1D> 
	getReverseIterator(const int _dim, const int _min, const int _max) const;

	DimensionIndexIterator<IteratorXD> 
	getIterator(const DimensionIndex & _activeDims, const DimensionIndex & _min, const DimensionIndex & _max) const;

	DimensionIndexIterator<IteratorXD> 
	getReverseIterator(const DimensionIndex & _activeDims, const DimensionIndex & _min, const DimensionIndex & _max) const;

	DimensionIndexIterator<IteratorList> 
	getIterator(const DimensionIndex &_activeset,  const DimensionIndexList & _idxlist) const;

	DimensionIndexIterator<IteratorList> 
	getReverseIterator(const DimensionIndex &_activeset, const DimensionIndexList & _idxlist) const;

	void Reverse();

	int length() const;

	int& operator[] (const int x) const;

	//DimensionIndex& 
	diinitializer operator=(const DimensionIndex &rhs);

	//DimensionIndex& 
	diinitializer operator=(const int arr);

	bool operator==(const DimensionIndex &other);
	 
};

}

#include <htucker/dimensionindex/diinitializer.h>


#include <htucker/dimensionindex/dimensionindexlist.h>

#include <htucker/dimensionindex/dimensionindexiterator.h>

#include <htucker/dimensionindex/dimensionindex.tcc>

#include <htucker/dimensionindex/dimensionindexlist.tcc>

#include <htucker/dimensionindex/dimensionindexiterator.tcc>

#endif // HTUCKER_DIMENSIONINDEX__DIMENSIONINDEX_H