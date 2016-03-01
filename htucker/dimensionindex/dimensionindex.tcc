#include <cassert>

namespace htucker{

DimensionIndex::DimensionIndex():index(NULL), len(0){}

DimensionIndex::DimensionIndex(const int  _d): len(_d){
	assert(_d >= 0);
	index = new int[_d];	
	for(int i = 0; i< _d; ++i){
		index[i] = 0;
	}
}

DimensionIndex::DimensionIndex(const int * _index,const int  _d){
	len = _d;
	index = new int[len];
	for(int i = 0; i< len; ++i){
		index[i] = _index[i];
	}
}

DimensionIndex::DimensionIndex(const int  constval,const int  _d){
	len = _d;
	index = new int[len];
	for(int i = 0; i< len; ++i){
		index[i] = constval;
	}
}

DimensionIndex::DimensionIndex(const DimensionIndex &copy){
	len = copy.length();
	index = new int[len];
	if(len > 0){
		for(int i = 0; i< len; ++i){
			index[i] = copy[i];
		}
	}
}


void
DimensionIndex::setRandom(const int  _min, const int  _max){
	for(int i = 0; i< len; ++i){
		index[i] = (rand() % (_max - _min + 1)) + _min;
	}
}

void
DimensionIndex::setRandom(const DimensionIndex &_min, const DimensionIndex &_max){
	assert(_min.length() == this->length() && _max.length() == this->length());

	for(int i = 0; i< this->length(); ++i){
		index[i] = (rand() % (_max[i] - _min[i] + 1)) + _min[i];
	}
}

void DimensionIndex::setRandom(const int  _pos, const int  _min, const int  _max){
	assert(_pos >= 0 && _pos < this->length());
	index[_pos - 1] = (rand() % (_max - _min + 1)) + _min;
}

void DimensionIndex::setRandom(const DimensionIndex & _dimidx, const DimensionIndex &_min, const DimensionIndex &_max){
	assert(this->length() == _min.length() && this->length() == _max.length());
	for(int i = 0; i< _dimidx.length(); ++i){
		index[_dimidx[i]-1] = (rand() % (_max[_dimidx[i]-1] - _min[_dimidx[i]-1] + 1)) + _min[_dimidx[i]-1];
	}
}

void DimensionIndex::setRandom(const DimensionIndexList &_dimlist){
	int pos = rand() % _dimlist.length();
	assert(_dimlist[pos].length() == len);
	for(int i = 0; i < len; ++i){
		index[i] = _dimlist[pos][i];
	}
}

void DimensionIndex::setRandom(const DimensionIndexList &_dimlist,const DimensionIndex &_activedims){
	DimensionIndex save = _dimlist[rand() % _dimlist.length()];
	for(int i = 0; i< _activedims.length(); ++i){
		index[_activedims[i] - 1] = save[_activedims[i] - 1];
	}
}

void 
DimensionIndex::setValue(const int _value){
	for(int i = 0; i< len; ++i){
		index[i] = _value;
	}
}

void 
DimensionIndex::setValue(const int _pos, const int _value){
	index[_pos - 1] = _value;
}

void
DimensionIndex::setValue(const flens::DenseVector<flens::Array<int> >& vals)
{
    typedef flens::DenseVector<flens::Array<int> >::IndexType  IndexType;

    assert(vals.length() <= len);
    IndexType first = vals.firstIndex();
    for (IndexType i=vals.firstIndex(); i!=vals.endIndex(); i+=vals.inc()) {
        index[(int)(i-first)] = vals(i);
    }
}

void 
DimensionIndex::setValue(const DimensionIndex &_values, const DimensionIndex &_activedims){
	for(int i = 0; i < _activedims.length(); ++i){
		index[_activedims[i] - 1] = _values[_activedims[i]-1];	
	}
}

void 
DimensionIndex::setValueAsc(){
	for(int i = 0; i< len; ++i){
		index[i] = i+1;
	}
}

void 
DimensionIndex::setValueAsc(const int _min){
	for(int i = 0; i < len; i++){
		index[i] = i + _min;
	}
}

bool 
DimensionIndex::equals(const DimensionIndex &_anotherindex, const DimensionIndex &_compareDims) const{
	for(int i = 0; i < _compareDims.length(); ++i){
		if(_anotherindex[_compareDims[i]-1] != this->operator[](_compareDims[i]-1)){
			return false;
		}
	}
	return true;
}

bool 
DimensionIndex::equals(const DimensionIndexList &_dimlist, const DimensionIndex &_compareDims) const{
	for(int i = 0; i < _dimlist.length(); ++i){
		if(this->equals(_dimlist[i],_compareDims)){
			return true;
		}
	}
	return false;
}

int 
DimensionIndex::computeBEvalue(const int _min, const int _max) const{
	int value = 0;
	for(int i = 0; i <len; ++i){
		if(i > 0){
			value = value * (_max - _min + 1);
		}
		value = value + (this->operator [](i) - _min);
	}
	return value + 1;
}

int
DimensionIndex::computeBEvalue(const DimensionIndex &_min, const DimensionIndex &_max) const {
	int value = 0;
	for(int i = 0; i <len; ++i){
		if(i > 0){
			value = value * (_max[i] - _min[i] + 1);
		}
		value = value + (this->operator [](i) - _min[i]);
	}
	return value + 1;
}

int 
DimensionIndex::computeLEvalue(const int _min, const int _max) const {
	int value = 0;
	for(int i = len - 1; i >= 0; --i){
		if(i < len-1){
			value = value * (_max - _min + 1);
		}
		value = value + (this->operator [](i) - _min);
	}
	return value + 1;
}

int
DimensionIndex::computeLEvalue(const DimensionIndex &_min, const DimensionIndex &_max) const{
	int value = 0;
	for(int i = len - 1; i >= 0; --i){
		if(i < len - 1){
			value = value * (_max[i] - _min[i] + 1);
		}
		value = value + (this->operator [](i) - _min[i]);
	}
	return value + 1;
}

DimensionIndex
DimensionIndex::getComplement(const int dim) const{
	assert(dim >= len);

	DimensionIndex ret(dim - len);
	DimensionIndex lon(dim);
	for(int i = 0; i < dim ; ++i){
		lon[i] = i+1;
	}
	for(int i = 0; i < len; ++i){
		lon[index[i]-1] = 0;
	}
	int pos = 0;
	for(int i = 0 ; i < dim; ++i){
		if(lon[i] != 0){
			ret[pos] = lon[i];
			pos++;
		}
	}
	return ret;
}

DimensionIndex
DimensionIndex::join(const DimensionIndex & rhs) const{
	DimensionIndex ret(len+rhs.length());
	for(int i = 0; i < len; ++i){
		ret[i] = index[i];
	}
	for(int i = 0; i < rhs.length(); ++i){
		ret[i + len] = rhs[i];
	}
	return ret;
}

DimensionIndex
DimensionIndex::intersect(const DimensionIndex &rhs) const{
	DimensionIndex neu(std::max(rhs.length(),len));
	int count = 0;
	for(int i = 0; i < len; ++i){
		for(int j = 0; j < rhs.length(); ++j){
			if(index[i] == rhs[j]){
				neu[count] = index[i];
				count++;
				break;
			}
		}
	}
	DimensionIndex ret(count);
	for(int i = 0; i < count; ++i){
		ret[i] = neu[i];
	}
	return ret;
}

DimensionIndexIterator<Iterator1D> 
DimensionIndex::getIterator(const int _dim, const int _min, const int _max) const{
	DimensionIndexIterator<Iterator1D> ret(*this,_dim,_min,_max);
	ret.setFirst();
	return ret;
}

DimensionIndexIterator<Iterator1D> 
DimensionIndex::getReverseIterator(const int _dim, const int _min, const int _max) const{
	DimensionIndexIterator<Iterator1D> ret(*this,_dim,_min,_max);
	ret.setLast();
	return ret;
}

DimensionIndexIterator<IteratorXD> 
DimensionIndex::getIterator(const DimensionIndex & _activeDims, const DimensionIndex & _min, const DimensionIndex & _max) const{
	DimensionIndexIterator<IteratorXD> ret(*this,_activeDims, _min, _max);
	ret.setFirst();
	return ret;
}

DimensionIndexIterator<IteratorXD> 
DimensionIndex::getReverseIterator(const DimensionIndex & _activeDims, const DimensionIndex & _min, const DimensionIndex & _max) const{
	DimensionIndexIterator<IteratorXD> ret(*this,_activeDims, _min, _max);
	ret.setLast();
	return ret;
}

DimensionIndexIterator<IteratorList> 
DimensionIndex::getIterator(const DimensionIndex &_activeset,  const DimensionIndexList & _idxlist) const{
	DimensionIndexIterator<IteratorList> ret(*this, _activeset, _idxlist);
	ret.setFirst();
	return ret;
}

DimensionIndexIterator<IteratorList> 
DimensionIndex::getReverseIterator(const DimensionIndex &_activeset,  const DimensionIndexList & _idxlist) const{
	DimensionIndexIterator<IteratorList> ret(*this, _activeset, _idxlist);
	ret.setLast();
	return ret;
}

void 
DimensionIndex::Reverse(){
	DimensionIndex save = *this;
	for(int i = len-1;i >= 0; i--){
		index[len-1 - i] = save[i];
	}
}





int 
DimensionIndex::length() const{
	return len;
}

int& 
DimensionIndex::operator[] (const int x) const{
	assert(x >= 0 && x < this->length());
	return index[x];
}

//DimensionIndex&
diinitializer
DimensionIndex::operator=(const DimensionIndex &rhs){

	if(&rhs == this){
		//return *this;
	} else {
		if(rhs.length() == this->len){
			for(int i = 0; i < len; ++i){
				index[i] = rhs[i];
			}
		} else {
			//std::cout << index << std::endl;
			//std::cout << NULL << std::endl;
			if(index != NULL){
				delete [] index;
			}
			this->len = rhs.length();
			index = new int[len];
			for(int i = 0; i < len; ++i){
				index[i] = rhs[i];
			}
		}
		//return *this;
	}
	return diinitializer(*this);
}

//DimensionIndex&
diinitializer 
DimensionIndex::operator=(const int  arr){
	diinitializer init(*this,arr);
	return init;
	
	
}

bool 
DimensionIndex::operator==(const DimensionIndex &other){
	if(other.length() != len){
		return false;
	}
	for(int i = 0; i < len; ++i){
		if(other[i]  != index[i]){
			return false;
		}
	}
	return true;
}




std::ostream& operator<<(std::ostream& Stream, const DimensionIndex& B)
{
	Stream << "{";
	for(int i = 0; i< B.length() - 1; ++i){
		Stream << B[i] << ", ";
	}
	if(B.length() > 0){
		Stream << B[B.length() - 1] ;
	}
	Stream <<  "} ";
	return Stream;
}


DimensionIndex::~DimensionIndex(){
	//std::cout << (*this) <<  "    Destruct Index" << std::endl;
	delete [] index;
	
	//index = NULL;
}


} // namespace htucker
