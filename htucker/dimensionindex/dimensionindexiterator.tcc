namespace htucker{

DimensionIndexIterator<Iterator1D>::DimensionIndexIterator(const DimensionIndex &_idx, const int _dim, const int _min, const int _max):idx(_idx), dim(_dim), maxit(_max), minit(_min){
		assert(idx.length() > _dim - 1 && _dim > 0);
};

DimensionIndexIterator<Iterator1D>::DimensionIndexIterator(const DimensionIndexIterator<Iterator1D> &iter):idx(iter.getIndex()),dim(iter.getDim()),maxit(iter.getMax()),minit(iter.getMin()){};

void
DimensionIndexIterator<Iterator1D>::operator ++(int){
	idx[dim - 1]++;
};

void
DimensionIndexIterator<Iterator1D>::operator --(int){
	idx[dim - 1]--;
};

void 
DimensionIndexIterator<Iterator1D>::setFirst(){
	idx[dim-1] = minit;
};

void 
DimensionIndexIterator<Iterator1D>::setLast(){
	idx[dim-1] = maxit;
};

bool 
DimensionIndexIterator<Iterator1D>::inRange() const{
	return (idx[dim-1] >= minit && idx[dim-1] <= maxit);
};

const DimensionIndex &
DimensionIndexIterator<Iterator1D>::getIndex() const{
	return idx;
};

int 
DimensionIndexIterator<Iterator1D>::getMin() const{
	return minit;
};

int 
DimensionIndexIterator<Iterator1D>::getMax() const{
	return maxit;
};

int 
DimensionIndexIterator<Iterator1D>::getDim() const{
	return dim;
};

DimensionIndexIterator<Iterator1D>& 
DimensionIndexIterator<Iterator1D>::operator=(const DimensionIndexIterator<Iterator1D> &rhs){
	if(&rhs != this){
		idx = rhs.getIndex();
		minit = rhs.getMin();
		maxit = rhs.getMax();
		dim = rhs.getDim();
	}
	return *this;
};



//___________________________________________________________________________________________




DimensionIndexIterator<IteratorXD>::DimensionIndexIterator(const DimensionIndex & _idx, const DimensionIndex & _activeDims, const DimensionIndex & _min, const DimensionIndex & _max):idx(_idx), active(_activeDims), maxit(_max), minit(_min){
	assert(maxit.length() == idx.length() && minit.length() == idx.length());
};

DimensionIndexIterator<IteratorXD>::DimensionIndexIterator(const DimensionIndexIterator<IteratorXD> &iter):idx(iter.getIndex()),active(iter.getActiveDims()),maxit(iter.getMax()),minit(iter.getMin()){
	assert(maxit.length() == idx.length() && minit.length() == idx.length());
};

void
DimensionIndexIterator<IteratorXD>::operator ++(int){
	/*idx[active[0]-1]++;
	for(int i = 0; i< active.length() - 1; ++i){
		if(idx[active[i]-1] > maxit[active[i]-1]){
			idx[active[i]-1] = minit[active[i]-1];
			idx[active[i+1]-1] ++;
		}
	}*/

	idx[active[active.length()-1]-1]++;
	for(int i = active.length()-1; i >= 0; --i){
		if(idx[active[i]-1] > maxit[active[i]-1]){
			if(i > 0){
			   idx[active[i]-1] = minit[active[i]-1];
			   idx[active[i-1]-1] ++;
			} 
		}
	}
};

void
DimensionIndexIterator<IteratorXD>::operator --(int){
	/*idx[active[0]-1] --;
	for(int i = 0; i<active.length() -1; ++i){
		if(idx[active[i]-1] < minit[active[i]-1]){
			idx[active[i]-1] = maxit[active[i]-1];
			idx[active[i+1]-1] --;
		}
	}*/
	idx[active[active.length()-1]-1]--;
	for(int i = active.length()-1; i >= 0; --i){
		if(idx[active[i]-1] < minit[active[i]-1]){
			if(i > 0){
			   idx[active[i]-1] = maxit[active[i]-1];
			   idx[active[i-1]-1] --;
			} 
		}
	}
};

void 
DimensionIndexIterator<IteratorXD>::setFirst(){
	for(int i = 0; i < active.length(); ++i){
		idx[active[i]-1] = minit[active[i]-1];
	}
};

void 
DimensionIndexIterator<IteratorXD>::setLast(){
	for(int i = 0; i < active.length(); ++i){
		idx[active[i]-1] = maxit[active[i]-1];
	}
};

bool 
DimensionIndexIterator<IteratorXD>::inRange() const{
	for(int i = 0; i< active.length(); ++i){
		if(!(idx[active[i]-1] >= minit[active[i]-1] && idx[active[i]-1] <= maxit[active[i]-1])){
			return false;
		}
	}
	return true;
};

const DimensionIndex & 
DimensionIndexIterator<IteratorXD>::getIndex() const{
	return idx;
};

const DimensionIndex &
DimensionIndexIterator<IteratorXD>::getMin() const{
	return minit;
};

const DimensionIndex &
DimensionIndexIterator<IteratorXD>::getMax() const{
	return maxit;
};

const DimensionIndex &
DimensionIndexIterator<IteratorXD>::getActiveDims() const{
	return active;
};

DimensionIndexIterator<IteratorXD>& 
DimensionIndexIterator<IteratorXD>::operator=(const DimensionIndexIterator<IteratorXD> &rhs){
	if(&rhs != this){
		idx = rhs.getIndex();
		minit = rhs.getMin();
		maxit = rhs.getMax();
		active = rhs.getActiveDims();
	}
	return *this;
};




//___________________________________________________________________________________________




DimensionIndexIterator<IteratorXDALL>::DimensionIndexIterator(const DimensionIndex & _idx, const DimensionIndex & _min, const DimensionIndex & _max):idx(_idx), maxit(_max), minit(_min){
	assert(maxit.length() == idx.length() && minit.length() == idx.length());
};

DimensionIndexIterator<IteratorXDALL>::DimensionIndexIterator(const DimensionIndexIterator<IteratorXDALL> &iter):idx(iter.getIndex()),maxit(iter.getMax()),minit(iter.getMin()){
	assert(maxit.length() == idx.length() && minit.length() == idx.length());
};

void
DimensionIndexIterator<IteratorXDALL>::operator ++(int){
	idx[0]++;
	for(int i = 0; i< idx.length() - 1; ++i){
		if(idx[i] > maxit[i]){
			idx[i] = minit[i];
			idx[i+1] ++;
		}
	}
};

void
DimensionIndexIterator<IteratorXDALL>::operator --(int){
	idx[0] --;
	for(int i = 0; i<idx.length() -1; ++i){
		if(idx[i] < minit[i]){
			idx[i] = maxit[i];
			idx[i+1] --;
		}
	}
};

void 
DimensionIndexIterator<IteratorXDALL>::setFirst(){
	for(int i = 0; i < idx.length(); ++i){
		idx[i] = minit[i];
	}
};

void 
DimensionIndexIterator<IteratorXDALL>::setLast(){
	for(int i = 0; i < idx.length(); ++i){
		idx[i] = maxit[i];
	}
};

bool 
DimensionIndexIterator<IteratorXDALL>::inRange() const{
	for(int i = 0; i< idx.length(); ++i){
		if(!(idx[i] >= minit[i] && idx[i] <= maxit[i])){
			return false;
		}
	}
	return true;
};

const DimensionIndex & 
DimensionIndexIterator<IteratorXDALL>::getIndex() const{
	return idx;
};

const DimensionIndex &
DimensionIndexIterator<IteratorXDALL>::getMin() const{
	return minit;
};

const DimensionIndex &
DimensionIndexIterator<IteratorXDALL>::getMax() const{
	return maxit;
};


DimensionIndexIterator<IteratorXDALL>& 
DimensionIndexIterator<IteratorXDALL>::operator=(const DimensionIndexIterator<IteratorXDALL> &rhs){
	if(&rhs != this){
		idx = rhs.getIndex();
		minit = rhs.getMin();
		maxit = rhs.getMax();
	}
	return *this;
};







//___________________________________________________________________________________________




DimensionIndexIterator<IteratorListALL>::DimensionIndexIterator(const DimensionIndex & _idx,  const DimensionIndexList & _idxlist):idx(_idx), idxlist(_idxlist),pos(0){};

DimensionIndexIterator<IteratorListALL>::DimensionIndexIterator(const DimensionIndexIterator<IteratorListALL> &iter):idx(iter.getIndex()),idxlist(iter.getList()),pos(iter.getPos()){};

void
DimensionIndexIterator<IteratorListALL>::operator ++(int){
	pos ++;
	if(inRange()){
		idx = idxlist[pos];
	}
};

void
DimensionIndexIterator<IteratorListALL>::operator --(int){
	pos --;
	if(inRange()){
		idx = idxlist[pos];
	}
};

void 
DimensionIndexIterator<IteratorListALL>::setFirst(){
	idx = idxlist[0];
};

void 
DimensionIndexIterator<IteratorListALL>::setLast(){
	idx = idxlist[idxlist.length() -1];
};

bool 
DimensionIndexIterator<IteratorListALL>::inRange() const{
	return (pos >= 0 && pos < idxlist.length());
};

const DimensionIndex & 
DimensionIndexIterator<IteratorListALL>::getIndex() const{
	return idx;
};

const DimensionIndexList &
DimensionIndexIterator<IteratorListALL>::getList() const{
	return idxlist;
};

int
DimensionIndexIterator<IteratorListALL>::getPos() const{
	return pos;
};


DimensionIndexIterator<IteratorListALL>& 
DimensionIndexIterator<IteratorListALL>::operator=(const DimensionIndexIterator<IteratorListALL> &rhs){
	if(&rhs != this){
		idx = rhs.getIndex();
		idxlist = rhs.getList();
		pos = rhs.getPos();
	}
	return *this;
};



//___________________________________________________________________________________________




DimensionIndexIterator<IteratorList>::DimensionIndexIterator(const DimensionIndex & _idx, const DimensionIndex & _activedims, const DimensionIndexList & _idxlist):idx(_idx), active(_activedims),idxlist(_idxlist),pos(0){};

DimensionIndexIterator<IteratorList>::DimensionIndexIterator(const DimensionIndexIterator<IteratorList> &iter):idx(iter.getIndex()),active(iter.getActiveDims()),idxlist(iter.getList()),pos(iter.getPos()){};

void
DimensionIndexIterator<IteratorList>::operator ++(int){
	pos ++;
	if(inRange()){
		idx.setValue(idxlist[pos],active);
	}
};

void
DimensionIndexIterator<IteratorList>::operator --(int){
	pos --;
	if(inRange()){
		idx.setValue(idxlist[pos],active);
	}
};

void 
DimensionIndexIterator<IteratorList>::setFirst(){
	idx.setValue(idxlist[0],active);
	pos = 0;
};

void 
DimensionIndexIterator<IteratorList>::setLast(){
	idx.setValue(idxlist[idxlist.length() -1],active);
	pos = idxlist.length() - 1;
};

bool 
DimensionIndexIterator<IteratorList>::inRange() const{
	return (pos >= 0 && pos < idxlist.length());
};

const DimensionIndex & 
DimensionIndexIterator<IteratorList>::getIndex() const{
	return idx;
};

const DimensionIndex & 
DimensionIndexIterator<IteratorList>::getActiveDims() const{
	return active;
};


const DimensionIndexList &
DimensionIndexIterator<IteratorList>::getList() const{
	return idxlist;
};

int
DimensionIndexIterator<IteratorList>::getPos() const{
	return pos;
};


DimensionIndexIterator<IteratorList>& 
DimensionIndexIterator<IteratorList>::operator=(const DimensionIndexIterator<IteratorList> &rhs){
	if(&rhs != this){
		idx = rhs.getIndex();
		idxlist = rhs.getList();
		pos = rhs.getPos();
		active = rhs.getActiveDims();
	}
	return *this;
};


} //namespace lawa
