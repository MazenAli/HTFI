#ifndef HTUCKER_GENERALTREE_GENERALTREEITERATOR_H
#define HTUCKER_GENERALTREE_GENERALTREEITERATOR_H 1


#include <htucker/generaltree/generaltreenode.h>

namespace htucker{

//An Interator Class, to walk through the tree 
template <typename NType>
class GeneralTreeIterator{
	GeneralTreeNode<NType> * begin;
	GeneralTreeNode<NType> *levelstep,*lastnode;
public:

	GeneralTreeIterator();
	
	GeneralTreeIterator(GeneralTreeNode<NType> * _begin);

	void operator++(int);

	void operator--(int);

	bool operator<= (const GeneralTreeIterator<NType> &Right) const;

	bool operator>= (const GeneralTreeIterator<NType> &Right) const;

	bool operator== (const GeneralTreeIterator<NType> &Right) const;

	bool operator!= (const GeneralTreeIterator<NType> &Right) const;


	GeneralTreeNode<NType> *
	getNode() const;

};

} // namespace htucker

#include <htucker/generaltree/generaltreeiterator.tcc>

#endif // HTUCKER_GENERALTREE_GENERALTREEITERATOR_H
