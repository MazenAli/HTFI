

m2vindexconverter::m2vindexconverter(const DimensionIndex &_min, const DimensionIndex &_max):minidx(_min),maxidx(_max),rowinfo(false){}

m2vindexconverter::m2vindexconverter(const m2vindexconverter & copy):minidx(copy.getMin()),maxidx(copy.getInputMax()),rowinfo(copy.hasRowInfo()),maxidxrow(copy.getRowInfoMax()),minidxrow(copy.getRowInfoMin()){}


void 
m2vindexconverter::getMatrixIndices(const DimensionIndex &vals, DimensionIndex & row, DimensionIndex & column) const{
	DimensionIndex r(vals.length());
	DimensionIndex c(vals.length());
	//cout << minidx << endl;

	//cout << maxidx << endl;
	if(!rowinfo){
		for(int i = 0; i < vals.length(); ++i){
			c[i] = floor((vals[i] - minidx[i]+0.0)/(maxidx[i]-minidx[i]+1.0)) + minidx[i];
			r[i] = vals[i] - (c[i] - minidx[i])*(maxidx[i]-minidx[i]+1);
		}
		row = r;
		column = c;
	} else {
		for(int i = 0; i < vals.length(); ++i){
			c[i] = floor((vals[i] - minidx[i]+0.0)/(maxidxrow[i]-minidxrow[i]+1.0)) + minidx[i];
			r[i] = vals[i] - (c[i] - minidx[i])*(maxidxrow[i]-minidxrow[i]+1) - maxidx[i] + maxidxrow[i];
		}
		row = r;
		column = c;
	}
}




void
m2vindexconverter::getMatrixIndices(const DimensionIndex & vals, int & row, int & column) const{
	DimensionIndex r,c;
	getMatrixIndices(vals,r,c);
	if(rowinfo){
		row = r.computeBEvalue(minidxrow,maxidxrow);
	} else {
		row = r.computeBEvalue(minidx,maxidx);
	}
	column = c.computeBEvalue(minidx,maxidx);
}


void 
m2vindexconverter::setRowIndices(const DimensionIndex &_minidxrow, const DimensionIndex &_maxidxrow){
	assert(_minidxrow.length() == minidx.length() && _maxidxrow.length() == maxidx.length());
	rowinfo = true;
	minidxrow = _minidxrow;
	maxidxrow = _maxidxrow;
	
}

DimensionIndex
m2vindexconverter::getMax() const{
	DimensionIndex vecmax(maxidx.length());
	if(rowinfo){
		for(int i = vecmax.length()-1; i >= 0; --i){
			vecmax[i] = (maxidx[i]-minidx[i]+1)*(maxidxrow[i]-minidxrow[i]+1) + minidx[i]-1;		
		}
	} else {
		for(int i = vecmax.length()-1; i >= 0; --i){
			vecmax[i] = (maxidx[i]-minidx[i]+1)*(maxidx[i]-minidx[i]+1) + minidx[i]-1;		
		}
	}
	return vecmax;
}

const DimensionIndex&
m2vindexconverter::getMin() const{
	return minidx;
}

const DimensionIndex& 
m2vindexconverter::getInputMax() const{
	return maxidx;
}

const DimensionIndex&
m2vindexconverter::getRowInfoMin() const{
	return minidxrow;
}

const DimensionIndex&
m2vindexconverter::getRowInfoMax() const{
	return maxidxrow;
}

bool
m2vindexconverter::hasRowInfo() const{
	return rowinfo;
}

m2vindexconverter&
m2vindexconverter::operator=(const m2vindexconverter &copy){
	minidx = copy.getMin();
	maxidx = copy.getInputMax();
	minidxrow = copy.getRowInfoMin();
	maxidxrow = copy.getRowInfoMax();
	rowinfo = copy.hasRowInfo();
	return *this;
}

