namespace htucker{

DimensionIndexList::DimensionIndexList(){
	list = new DimensionIndex[10];
	size = 10;
	len = 0;
}

DimensionIndexList::DimensionIndexList(const DimensionIndex  &_firstElement){
	list = new DimensionIndex[10];
	size = 10;
	len = 0;
	list[0] = _firstElement;
	len = 1;
}

DimensionIndexList::DimensionIndexList(const DimensionIndexList &copylist){
	list = new DimensionIndex[10];
	size = 10;
	len = 0;
	for(int i = 0; i <copylist.length(); ++i){
		this->add(copylist[i]);
	}
}

DimensionIndexList::~DimensionIndexList(){
	delete [] list;
}



void
DimensionIndexList::add(const DimensionIndex &_newElement){
	if(len < size){
		list[len] = _newElement;
		len ++;
	} else {
		int newsize = size + 10;
		DimensionIndex* cpy = new DimensionIndex[newsize];
		for(int i = 0; i < len ; ++i){
			cpy[i] = list[i];
		}
		//memcpy(cpy, list, size*sizeof(DimensionIndex) );
		size = newsize;
		delete [] list;
		list = cpy;
		list[len] = _newElement;
		len ++;
	}
}

void
DimensionIndexList::remove(const int  _pos){
	assert(_pos >= 0 && _pos < len);
	for(int i = _pos; i < len - 1; ++i){
		list[i] = list[i+1];
	}
	len--;
}

void
DimensionIndexList::empty(){
	delete [] list;
	list = new DimensionIndex[10];
	size = 10;
	len = 0;
}

int
DimensionIndexList::length() const{
	return len;
}

// This function not really const // Mazen
DimensionIndex& 
DimensionIndexList::operator[] (const int x) const{
	assert(x>= 0 && x < len);
	return list[x];
}

DimensionIndexList &
DimensionIndexList::operator= (const DimensionIndexList & copylist){
	if(this->length() > 0){
		this->empty();
	}
	for(int i = 0; i < copylist.length(); ++i){
		this->add(copylist[i]);
	}
	return *this;
}

std::ostream& operator<<(std::ostream& Stream, const DimensionIndexList& B)
{
		Stream << "{" ;
		for(int i = 0; i< B.length() - 1; ++i){
			Stream << i << ": " << B[i] << ", " << std::endl;
		}
		if(B.length() > 0){
			Stream << (B.length() - 1) << ": " << B[B.length() - 1] ;
		}
		Stream << "} ";
	return Stream;
}

} // namespace htucker

