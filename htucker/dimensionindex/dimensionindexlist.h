#ifndef HTUCKER_DIMENSIONINDEXLIST_H
#define HTUCKER_DIMENSIONINDEXLIST_H 1


namespace htucker{

class DimensionIndexList{
	DimensionIndex * list;
	int len;
	int size;
public:
	DimensionIndexList();

	DimensionIndexList(const DimensionIndex &_firstElement);

	DimensionIndexList(const DimensionIndexList &copylist);

	~DimensionIndexList();

	void add(const DimensionIndex  &_newElement);

	void remove(const int _pos);

	void empty();

	int length() const;

	DimensionIndex& operator[] (const int x) const;

	DimensionIndexList &
	operator= (const DimensionIndexList & copylist);

};

} //namespace htucker

#endif //HTUCKER_DIMENSIONINDEXLIST_H
