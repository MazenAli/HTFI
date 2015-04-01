namespace htucker{

	//Constructor 
	template <typename NType>
	GeneralTreeNode<NType>::GeneralTreeNode(const NType &_content):content(_content),parent(NULL),firstChild(NULL),lastChild(NULL),nextSibling(NULL),lastSibling(NULL),levelleft(NULL),levelright(NULL){};

	template <typename NType>
	GeneralTreeNode<NType>::GeneralTreeNode(const GeneralTreeNode<NType> &_copy):content(_copy.getContent()),parent(_copy.getParent()),firstChild(_copy.getfirstChild()),lastChild(_copy.getlastChild()),nextSibling(_copy.getnextSibling()),lastSibling(_copy.getpreviousSibling()), levelleft(_copy.getlevelleft()), levelright(_copy.getlevelright()){};
	
	//Set Pointers
	template <typename NType>
	void 
	GeneralTreeNode<NType>::setParent(const GeneralTreeNode<NType> * const _parent){
		parent = (GeneralTreeNode<NType> *) _parent;
	};

	template <typename NType>
	void 
	GeneralTreeNode<NType>::setfirstChild(const GeneralTreeNode<NType> * const _firstChild){
		firstChild = (GeneralTreeNode<NType> *) _firstChild;
	};

	template <typename NType>
	void 
	GeneralTreeNode<NType>::setlastChild(const GeneralTreeNode<NType> * const _lastChild){
		lastChild = (GeneralTreeNode<NType> *) _lastChild;
	};

	template <typename NType>
	void 
	GeneralTreeNode<NType>::setnextSibling(const GeneralTreeNode<NType> * const _nextSibling){
		nextSibling = (GeneralTreeNode<NType> *) _nextSibling;
	};

	template <typename NType>
	void 
	GeneralTreeNode<NType>::setpreviousSibling(const GeneralTreeNode<NType> * const _lastSibling){
		lastSibling = (GeneralTreeNode<NType> *) _lastSibling;
	};

	template <typename NType>
	void 
	GeneralTreeNode<NType>::setlevelleft(const GeneralTreeNode<NType> * const _levelleft){
		levelleft = (GeneralTreeNode<NType> *) _levelleft;
	};

	template <typename NType>
	void 
	GeneralTreeNode<NType>::setlevelright(const GeneralTreeNode<NType> * const _levelright){
		levelright = (GeneralTreeNode<NType> *) _levelright;
	};

	template <typename NType>
	void 
	GeneralTreeNode<NType>::setContent(const NType &_content){
		content = _content;
	};


	//Get Pointer

	template <typename NType>
	GeneralTreeNode<NType> * const
	GeneralTreeNode<NType>::getParent() const{
		return parent;
	};

	template <typename NType>
	GeneralTreeNode<NType> * const
	GeneralTreeNode<NType>::getfirstChild() const{
		return firstChild;
	};

	template <typename NType>
	GeneralTreeNode<NType> * const
	GeneralTreeNode<NType>::getlastChild() const{
		return lastChild;
	};

	template <typename NType>
	GeneralTreeNode<NType> * const
	GeneralTreeNode<NType>::getnextSibling() const{
		return nextSibling;
	};

	template <typename NType>
	GeneralTreeNode<NType> * const
	GeneralTreeNode<NType>::getpreviousSibling() const{
		return lastSibling;
	};

	template <typename NType>
	GeneralTreeNode<NType> * const
	GeneralTreeNode<NType>::getlevelleft() const{
		return levelleft;
	};

	template <typename NType>
	GeneralTreeNode<NType> * const
	GeneralTreeNode<NType>::getlevelright() const{
		return levelright;
	};
	
	template <typename NType>
	NType * 
	GeneralTreeNode<NType>::getContent(){
		return  &content;
	};

	
	//Append
	template <typename NType>
	void 
	GeneralTreeNode<NType>::appendChild(const NType &_content){
		GeneralTreeNode<NType> * newnode = new GeneralTreeNode<NType>(_content);
		newnode->parent = this;
		if(this->firstChild == NULL || this->lastChild == NULL){
			newnode->nextSibling = newnode;
			newnode->lastSibling = newnode;
			
			GeneralTreeNode<NType> * lsave;
			lsave = this->levelright;
			while(lsave != NULL){
				if(lsave->firstChild != NULL){
					lsave->firstChild->levelleft = newnode;
					newnode->levelright = lsave->firstChild;
					break;
				}
				lsave = lsave->levelright;
			}
			lsave = this->levelleft;
			while(lsave != NULL){
				if(lsave->lastChild != NULL){
					lsave->lastChild->levelright = newnode;
					newnode->levelleft = lsave->lastChild;
					break;
				}
				lsave = lsave->levelleft;
			}

			this->firstChild = newnode;
			this->lastChild = newnode;
			
		} else {
			newnode->lastSibling = this->lastChild;
			newnode->nextSibling = this->firstChild;
			newnode->levelright = this->lastChild->levelright;
			newnode->levelleft = this->lastChild;
			if(this->lastChild->levelright != NULL){
				this->lastChild->levelright->levelleft = newnode;
			}
			this->lastChild->levelright = newnode;
			this->lastChild->nextSibling = newnode;
			this->firstChild->lastSibling = newnode;
			this->lastChild = newnode;
		}
	};


	//Remove
	template <typename NType>
	bool 
	GeneralTreeNode<NType>::removeChild(const int pos){	
		//std::cout << "removeChild(" << pos << ")   level = " << this->level() << std::endl;
		if(this == NULL){ //this->firstChild
			//std::cout <<  "this->firstChild == NULL" << std::endl;
			return false;
		} else {
			GeneralTreeNode<NType> * Ssave;
			Ssave = this->firstChild;
			int count = 1;
			do {
				if(count == pos){
					if(Ssave->getfirstChild() != NULL){
						//Da steht noch was drin.... also alles leeren...
						
						while(Ssave->firstChild != NULL){
							if(! Ssave->removeChild(1)){
								return false;
							}
						}
					} 
					if(this->firstChild != this->lastChild){
						Ssave->lastSibling->nextSibling = Ssave->nextSibling;
						Ssave->nextSibling->lastSibling = Ssave->lastSibling;
					}
					if(Ssave->levelright != NULL){
						Ssave->levelright->levelleft = Ssave->levelleft;	
					}
					if(Ssave->levelleft != NULL){
						Ssave->levelleft->levelright = Ssave->levelright;
					}
					if(this->firstChild != this-> lastChild){
						if(Ssave == this->firstChild){
							this->firstChild = Ssave->nextSibling;
						} 
						if(Ssave == this->lastChild){
							this->lastChild = Ssave->lastSibling;
						}
					} else {
						this->firstChild = NULL;
						this->lastChild = NULL;
					}
					delete Ssave;
					

					return true;
				}
				count = count + 1;
				Ssave = Ssave->nextSibling;
			} while(Ssave != this->firstChild);
			return false;
		}
		
	};




	template <typename NType>
	bool 
	GeneralTreeNode<NType>::isLeaf() const{
		if(this->firstChild == NULL){
			return true;
		} else {
			return false;
		}
	};

	template <typename NType>
	bool 
	GeneralTreeNode<NType>::isRoot() const{
		if(this->parent == NULL){
			return true;
		} else {
			return false;
		}
	};

	template <typename NType>
	bool 
	GeneralTreeNode<NType>::isInner() const{
		if((!this->isRoot()) && (! this->isLeaf())){
			return true;
		} else {
			return false;
		}
	};

	template <typename NType>
	int
	GeneralTreeNode<NType>::level() const{
		GeneralTreeNode<NType> * save;
		int count = 1;
		if(this->parent != NULL){
			save = this->parent;
			count = 2;
		} else {
			return count;
		};
		while(save->parent != NULL){
			count += 1;
			save = save->parent;
		}
		return count;
	}


	template <typename NType>
	void
	GeneralTreeNode<NType>::printnode() const{
		int depth = this->level();
		for(int i = 1; i< depth; ++i){
			std::cout << "   ";
		}
		std::cout << "|-- " << this->content << std::endl;
		
		if(this->firstChild != NULL){
			GeneralTreeNode<NType> * save;
			save = this->firstChild;
			do {
				save->printnode();
				save = save-> nextSibling;
			} while(save != this->firstChild);
		}
			
	};


} // namespace lawa

