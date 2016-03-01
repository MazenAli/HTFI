#ifndef HTUCKER_DIMENSIONINDEX_DIMENSIONINDEXITERATOR_H
#define HTUCKER_DIMENSIONINDEX_DIMENSIONINDEXITERATOR_H 1


namespace htucker{

template <>
class DimensionIndexIterator<Iterator1D>{
	DimensionIndex idx;
	int dim, minit, maxit;
	
public:
	DimensionIndexIterator(const DimensionIndex & _idx, const int _dim, const int _min, const int _max);

	DimensionIndexIterator(const DimensionIndexIterator<Iterator1D> &iter);

	void operator++(int);

	void operator--(int);

	void setFirst();

	void setLast();

	bool inRange() const;

	const DimensionIndex & getIndex() const;

	int getMin() const;

	int getMax() const;

	int getDim() const;

	DimensionIndexIterator<Iterator1D>& operator=(const DimensionIndexIterator<Iterator1D> &rhs);

};

template <>
class DimensionIndexIterator<IteratorXD>{
	DimensionIndex idx,minit,maxit,active;	
public:
	DimensionIndexIterator(const DimensionIndex & _idx, const DimensionIndex & _activeDims, const DimensionIndex & _min, const DimensionIndex & _max);

	DimensionIndexIterator(const DimensionIndexIterator<IteratorXD> &iter);

	void operator++(int);

	void operator--(int);

	void setFirst();

	void setLast();

	bool inRange() const;

	const DimensionIndex & getIndex() const;

	const DimensionIndex & getMin() const;

	const DimensionIndex & getMax() const;

	const DimensionIndex & getActiveDims() const;

	DimensionIndexIterator<IteratorXD>& operator=(const DimensionIndexIterator<IteratorXD> &rhs);

};

template <>
class DimensionIndexIterator<IteratorXDALL>{
	DimensionIndex idx,minit,maxit;	
public:
	DimensionIndexIterator(const DimensionIndex & _idx,  const DimensionIndex & _min, const DimensionIndex & _max);

	DimensionIndexIterator(const DimensionIndexIterator<IteratorXDALL> &iter);

	void operator++(int);

	void operator--(int);

	void setFirst();

	void setLast();

	bool inRange() const;

	const DimensionIndex & getIndex() const;

	const DimensionIndex & getMin() const;

	const DimensionIndex & getMax() const;

	DimensionIndexIterator<IteratorXDALL>& operator=(const DimensionIndexIterator<IteratorXDALL> &rhs);

};


template <>
class DimensionIndexIterator<IteratorListALL>{
	DimensionIndex idx;
	DimensionIndexList idxlist;
	int pos;
public:
	DimensionIndexIterator(const DimensionIndex & _idx, const DimensionIndexList & _idxlist);

	DimensionIndexIterator(const DimensionIndexIterator<IteratorListALL> &iter);

	void operator++(int);

	void operator--(int);

	void setFirst();

	void setLast();

	bool inRange() const;

	const DimensionIndex & getIndex() const;

	const DimensionIndexList & getList() const;

	int getPos() const;

	DimensionIndexIterator<IteratorListALL>& operator=(const DimensionIndexIterator<IteratorListALL> &rhs);

};



template <>
class DimensionIndexIterator<IteratorList>{
	DimensionIndex idx,active;	
	DimensionIndexList idxlist;
	int pos;
public:
	DimensionIndexIterator(const DimensionIndex & _idx, const DimensionIndex &_activeset,  const DimensionIndexList & _idxlist);

	DimensionIndexIterator(const DimensionIndexIterator<IteratorList> &iter);

	void operator++(int);

	void operator--(int);

	void setFirst();

	void setLast();

	bool inRange() const;

	const DimensionIndex & getIndex() const;

	const DimensionIndex & getActiveDims() const;

	const DimensionIndexList & getList() const;

	int getPos() const;

	DimensionIndexIterator<IteratorList>& operator=(const DimensionIndexIterator<IteratorList> &rhs);

};

}

#endif //HTUCKER_DIMENSIONINDEX_DIMENSIONINDEXITERATOR_H
