namespace htucker{


template <typename NType>
GeneralTree<NType>::GeneralTree():root(NULL){};

template <typename NType>
GeneralTree<NType>::GeneralTree(const NType &_rootcontent){
	root = new GeneralTreeNode<NType>(_rootcontent);
};

template <typename NType>
GeneralTree<NType>::GeneralTree(const GeneralTree<NType> &_copy){
	root = new GeneralTreeNode<NType>(*(_copy.root->getContent()));
	GeneralTreeNode<NType> * newbase = this->root;
	GeneralTreeNode<NType> * base = _copy.root;
	GeneralTreeNode<NType> * save;
	GeneralTreeNode<NType> * nextlev = NULL;
	GeneralTreeNode<NType> * newnextlev = NULL;
	while(1){
		save = base->getfirstChild();
		if(save != NULL){
			do{
				newbase->appendChild(*(save->getContent()));
				if(nextlev == NULL && newnextlev == NULL){
					nextlev = save;
					newnextlev = newbase->getfirstChild();
				}
				save = save->getnextSibling();
			} while(save != base->getfirstChild());
		}
		if(base->getlevelright() != NULL){
			base = base->getlevelright();
			newbase = newbase->getlevelright();
		} else if(nextlev != NULL && newnextlev != NULL){
			base = nextlev;
			newbase = newnextlev;
			nextlev = NULL;
			newnextlev = NULL;
		} else {
			break;
		}
	}
};


template <typename NType>
template <typename NType2>
GeneralTree<NType>::GeneralTree(const GeneralTree<NType2> &_copy){
	NType emptynode;
	root = new GeneralTreeNode<NType>(emptynode);
	GeneralTreeNode<NType> * newbase = this->root;
	GeneralTreeNode<NType2> * base = _copy.root;
	GeneralTreeNode<NType2> * save;
	GeneralTreeNode<NType2> * nextlev = NULL;
	GeneralTreeNode<NType> * newnextlev = NULL;
	while(1){
		save = base->getfirstChild();
		if(save != NULL){
			do{
				newbase->appendChild(emptynode);
				if(nextlev == NULL && newnextlev == NULL){
					nextlev = save;
					newnextlev = newbase->getfirstChild();
				}
				save = save->getnextSibling();
			} while(save != base->getfirstChild());
		}
		if(base->getlevelright() != NULL){
			base = base->getlevelright();
			newbase = newbase->getlevelright();
		} else if(nextlev != NULL && newnextlev != NULL){
			base = nextlev;
			newbase = newnextlev;
			nextlev = NULL;
			newnextlev = NULL;
		} else {
			break;
		}
	}
};


template <typename NType>
GeneralTree<NType>::~GeneralTree(){
	empty();
};

template <typename NType>
void 
GeneralTree<NType>::print() const{
	if(root != NULL){
		root->printnode();
	} else {
		cout << "empty GeneralTree" << endl;
	}
};

template <typename NType>
void 
GeneralTree<NType>::empty(){
	if(this->root != NULL){
		while(this->root->getfirstChild() != NULL){
			this->root->removeChild(1);
		}
		delete root;
		root = NULL;
	}
};

template <typename NType>
GeneralTreeIterator<NType>
GeneralTree<NType>::begin() const{
	GeneralTreeIterator<NType> ret(root);
	return ret;
};

template <typename NType>
GeneralTreeIterator<NType>
GeneralTree<NType>::end() const{
	
	GeneralTreeNode<NType> * save = root;
	
	while(1){
		if(save->getfirstChild() != NULL){
			save = save->getfirstChild();
			while(save->getlevelleft() != NULL){
				save = save->getlevelleft();
			}
		} else if(save->getlevelright() != NULL) {
			save = save->getlevelright();
		} else {
			break;
		}
	}

	GeneralTreeIterator<NType> ret(save);
	return ret;
};

template <typename NType>
GeneralTree<NType> & 
GeneralTree<NType>::operator= (const GeneralTree<NType> &_copy){
	if(this == &_copy){
		return *this;
	}
	
	if(this->root != NULL){
		this->empty();
	}
	root = new GeneralTreeNode<NType>(*(_copy.root->getContent()));
	GeneralTreeNode<NType> * newbase = this->root;
	GeneralTreeNode<NType> * base = _copy.root;
	GeneralTreeNode<NType> * save;
	GeneralTreeNode<NType> * nextlev = NULL;
	GeneralTreeNode<NType> * newnextlev = NULL;
	while(1){
		save = base->getfirstChild();
		if(save != NULL){
			do{
				newbase->appendChild(*(save->getContent()));
				if(nextlev == NULL && newnextlev == NULL){
					nextlev = save;
					newnextlev = newbase->getfirstChild();
				}
				save = save->getnextSibling();
			} while(save != base->getfirstChild());
		}
		if(base->getlevelright() != NULL){
			base = base->getlevelright();
			newbase = newbase->getlevelright();
		} else if(nextlev != NULL && newnextlev != NULL){
			base = nextlev;
			newbase = newnextlev;
			nextlev = NULL;
			newnextlev = NULL;
		} else {
			break;
		}
	}
	return *this;
	
};

} //  namespace lawa
