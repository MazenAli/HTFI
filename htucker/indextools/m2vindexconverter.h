#ifndef HTUCKER_INDEXTOOLS_M2VINDEXCONVERTER_H
#define HTUCKER_INDEXTOOLS_M2VINDEXCONVERTER_H 1

using namespace flens;

namespace htucker{

	//class for the conversion of an index of a vector into a matrix index and vice versa

	class m2vindexconverter{
		DimensionIndex minidx,maxidx;
		DimensionIndex minidxrow,maxidxrow;
		bool rowinfo;

		public:
			
			// We assume that the Index Structure of columns and rows are identical, _min and _max 
			// are the bounds of the rows or the columns we can compute the max of the vector as
			// (max[i] - min[i] + 1)^2 - 1 + min[i], if not, set numberofcolumnentrys

			//If they are not identical, we can set the values for the rows separately
			m2vindexconverter(){};
			
			m2vindexconverter(const DimensionIndex &_min, const DimensionIndex &_max);

			m2vindexconverter(const m2vindexconverter &copy);

			//braucht copy-Konstruktor

			void 
			getMatrixIndices(const DimensionIndex &vals, DimensionIndex &row, DimensionIndex &column) const;
			
			void
			getMatrixIndices(const DimensionIndex & vals, int & row, int & column) const;

			void 
			setRowIndices(const DimensionIndex &_minrows, const DimensionIndex &_maxrows);
			
			DimensionIndex
			getMax() const; //Max of the resulting vector!

			const DimensionIndex&
			getMin() const; //Min of the resulting vector!

			const DimensionIndex& 
			getInputMax() const;

			const DimensionIndex&
			getRowInfoMin() const;

			const DimensionIndex&
			getRowInfoMax() const;

			bool
			hasRowInfo() const;
			
			m2vindexconverter&
			operator=(const m2vindexconverter &copy);



	};

#include <htucker/indextools/m2vindexconverter.tcc>
} //end of namespace htucker

#endif // HTUCKER_INDEXTOOLS_QM2VINDEXCONVERTER_H
