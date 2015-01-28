namespace htucker{

template <typename NType>
GeneralTreeIterator<NType>::GeneralTreeIterator():begin(NULL),levelstep(NULL),lastnode(NULL){};

template <typename NType>
GeneralTreeIterator<NType>::GeneralTreeIterator(GeneralTreeNode<NType> * _begin):begin(_begin),levelstep(NULL),lastnode(NULL){};

template <typename NType>
void 
GeneralTreeIterator<NType>::operator++(int){
	lastnode = begin;
	if(begin->getfirstChild() != NULL && levelstep == NULL){
		levelstep = begin->getfirstChild();
	}
	if(begin->getlevelright() != NULL){
		begin = begin->getlevelright();
	} else {
		begin = levelstep; //this can be NULL
		levelstep = NULL;
	}
};


template <typename NType>
void 
GeneralTreeIterator<NType>::operator--(int){
	lastnode = begin;
	if(begin->getParent() != NULL && levelstep == NULL){
		levelstep = begin->getParent();
	}
	if(begin->getlevelleft() != NULL){
		begin = begin->getlevelleft();
	} else {
		begin = levelstep; //this can be NULL
		if(begin != NULL){
			while(begin->getlevelright() != NULL){
				begin = begin->getlevelright();
			}
		}
		levelstep = NULL;
	}
};


template <typename NType>
bool
GeneralTreeIterator<NType>::operator<= (const GeneralTreeIterator<NType> &Right) const{
	if(lastnode != Right.getNode()){	
		return true;
	} else {
		return false;
	}
};

template <typename NType>
bool
GeneralTreeIterator<NType>::operator>= (const GeneralTreeIterator<NType> &Right) const{
	if(lastnode != Right.getNode()){	
		return true;
	} else {
		return false;
	}
};

template <typename NType>
bool
GeneralTreeIterator<NType>::operator== (const GeneralTreeIterator<NType> &Right) const{
	if(Right.getNode() == begin){
		return true;
	} else {
		return false;
	}
};

template <typename NType>
bool
GeneralTreeIterator<NType>::operator!= (const GeneralTreeIterator<NType> &Right) const{
	if(Right.getNode() != begin){
		return true;
	} else {
		return false;
	}
};




template <typename NType>
GeneralTreeNode<NType> *
GeneralTreeIterator<NType>::getNode() const{
	return begin;
};


} // namespace htucker
