#ifndef HTUCKER_GENERALTREE_GENERALTREENODE_H
#define HTUCKER_GENERALTREE_GENERALTREENODE_H 1

namespace htucker{

template <typename NType>
class GeneralTreeNode{

private: 		
	GeneralTreeNode<NType> * parent;
	GeneralTreeNode<NType> * firstChild;
	GeneralTreeNode<NType> * lastChild;
	GeneralTreeNode<NType> * nextSibling;
	GeneralTreeNode<NType> * lastSibling;
	GeneralTreeNode<NType> * levelleft;
	GeneralTreeNode<NType> * levelright;
	NType content;
public:

	//Constructor 
	GeneralTreeNode(const NType &_content);

	GeneralTreeNode(const GeneralTreeNode<NType> &_copy);

	//Set Pointers
	void setParent(const GeneralTreeNode<NType> * _parent);

	void setfirstChild(const GeneralTreeNode<NType> *  const _firstChild);

	void setlastChild(const GeneralTreeNode<NType> *  const _lastChild);

	void setnextSibling(const GeneralTreeNode<NType> * const _nextSibling);

	void setpreviousSibling(const GeneralTreeNode<NType> * const _lastSibling);

	void setlevelleft(const GeneralTreeNode<NType> * const _levelleft);

	void setlevelright(const GeneralTreeNode<NType> *  const _levelright);

	void setContent(const NType &_content);

	//Get Pointer

	GeneralTreeNode<NType> *
	getParent() const;

	GeneralTreeNode<NType> *
	getfirstChild() const;

	GeneralTreeNode<NType> *
	getlastChild() const;

	GeneralTreeNode<NType> *
	getnextSibling() const;

	GeneralTreeNode<NType> *
	getpreviousSibling() const;

	GeneralTreeNode<NType> *
	getlevelleft() const;

	GeneralTreeNode<NType> *
	getlevelright() const;

	NType * getContent();

	//Append

	void appendChild(const NType &_content);

	//Remove

	bool removeChild(const int pos);


	bool isLeaf() const;

	bool isRoot() const;

	bool isInner() const;

	int level() const;

	//Output:

	void
	printnode() const;
	
};

} // namespace htucker

#include <htucker/generaltree/generaltreenode.tcc>

#endif // HTUCKER_GENERALTREE_GENERALTREENODE_H
