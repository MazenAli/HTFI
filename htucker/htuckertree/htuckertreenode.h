#ifndef HTUCKER_HTUCKERTREE_HTUCKERTREENODE_H
#define HTUCKER_HTUCKERTREE_HTUCKERTREENODE_H 1

#include <htucker/dimensionindex/dimensionindex.h>
//#include <htucker/rankoneapproximation/roa.h>



namespace htucker{

template <typename T>
class HTuckerTreeNode{
	private:
		DimensionIndex index;
		flens::GeMatrix<flens::FullStorage<T,flens::ColMajor> > UorB;
		//flens::DenseVector<flens::Array<T> > evaluate;
		int UorB_rcnumel, UorB_lcnumel, UorB_numel;

	public:
		HTuckerTreeNode();

		HTuckerTreeNode(const DimensionIndex &_index);

		HTuckerTreeNode(const HTuckerTreeNode<T> &_copy);

		HTuckerTreeNode &
		operator= (const HTuckerTreeNode<T> &_copy);

		DimensionIndex &
		getIndex() const;

		flens::GeMatrix<flens::FullStorage<T,flens::ColMajor> > &
		getUorB() const;

		void 
		setUorB(const flens::GeMatrix<flens::FullStorage<T,flens::ColMajor> > &_UorB);

		//flens::DenseVector<flens::Array<T> > & 
		//getEvaluate() const;

		//void 
		//setEvaluate(const flens::DenseVector<flens::Array<T> > & _evaluate);

		void setIndex(const DimensionIndex &_ind);

		int getNumRows() const;
		
		int getLeftChildNumRows() const;

		int getRightChildNumRows() const;

		void setNumRows(const int numel);

		void setLeftChildNumRows(const int lcnumel);

		void setRightChildNumRows(const int rcnumel);
};

} //namespace htucker

#include <htucker/htuckertree/htuckertreenode.tcc>


#endif // HTUCKER_HTUCKERTREE_HTUCKERTREENODE_H