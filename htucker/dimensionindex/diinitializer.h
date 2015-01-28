#ifndef HTUCKER_DIMENSIONINDEX_DIINITIALIZER_H
#define HTUCKER_DIMENSIONINDEX_DIINITIALIZER_H 1



namespace htucker{


class diinitializer{
	
	const DimensionIndex & di;
	int count;

	public:

	diinitializer(const DimensionIndex & _di);

	diinitializer(const DimensionIndex & _di, const int  val);

	diinitializer
	operator, (const int val);
	
};




} // namespace htucker

#include <htucker/dimensionindex/diinitializer.tcc>

#endif // HTUCKER_DIMENSIONINDEX_DIINITIALIZER_H
