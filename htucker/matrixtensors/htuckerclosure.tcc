namespace htucker{


template <typename op,class A, class B>
HTuckerClosure<op,A,B>::HTuckerClosure(const A & _left, const B & _right, const int _d):left(_left),right(_right),d(_d){
}


template <typename op,class A, class B>
A & 
HTuckerClosure<op,A,B>::getLeft() const{
	return left;
}



template <typename op,class A, class B>
B & 
HTuckerClosure<op,A,B>::getRight() const{
	return right;
}

template <typename op,class A, class B>
int 
HTuckerClosure<op,A,B>::dim() const{
	return d;
}



}
