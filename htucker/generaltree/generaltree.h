#ifndef HTUCKER_GENERALTREE_GENERALTREE_H
#define HTUCKER_GENERALTREE_GENERALTREE_H 1

#include <htucker/generaltree/generaltreeiterator.h>
#include <htucker/generaltree/generaltreenode.h>


namespace htucker{

//This is baisically a wrapper for GeneralTreeNode, containig some function that concern the whole tree
template <typename NType>
class GeneralTree{

public:
	GeneralTreeNode<NType> * root;

	GeneralTree();

	GeneralTree(const NType &_rootcontent);

	GeneralTree(const GeneralTree<NType> &_copy);

	template <typename NType2>
	GeneralTree(const GeneralTree<NType2> &_copy);

	~GeneralTree();

	void print() const;

	void empty();

	GeneralTreeIterator<NType>
	begin() const;

	GeneralTreeIterator<NType>
	end() const;

	GeneralTree<NType> & 
	operator= (const GeneralTree<NType> &_copy);

};

} //namespace htucker

#include <htucker/generaltree/generaltree.tcc>

#endif // HTUCKER_GENERALTREE_GENERALTREE_H
