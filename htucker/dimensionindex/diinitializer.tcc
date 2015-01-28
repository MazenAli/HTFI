namespace htucker{


diinitializer::diinitializer(const DimensionIndex & _di):di(_di),count(0){};

diinitializer::diinitializer(const DimensionIndex & _di, const int val):di(_di),count(0){
	assert(di.length() > 0);
	di[0] = val;
	count ++;
}


diinitializer
diinitializer::operator, (const int val){
	assert(count < di.length());
	di[count] = val;
	count ++;
	return *this;
}

}



