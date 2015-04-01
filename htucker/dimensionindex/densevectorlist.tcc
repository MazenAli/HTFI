#include <flens/flens.h>

//constructors
namespace htucker{

template <typename T>
DenseVectorList<T>::DenseVectorList(){
	node = NULL;
	lastnode = NULL;
	len = 0;
	iterator = NULL;
	pos = 0;
}

template <typename T>
DenseVectorList<T>::DenseVectorList(const flens::DenseVector<flens::Array<T> > &vector){
	node = new DenseVectorListElement<T>(vector);
	iterator = node;
	lastnode = node;
	len = 1;
	pos = 1;
};

template <typename T>
DenseVectorList<T>::~DenseVectorList(){
	//empty();
};


template <typename T>
void 
DenseVectorList<T>::add(const flens::DenseVector<flens::Array<T> > &vector){
	if(len > 0){
		lastnode->next = new DenseVectorListElement<T>(vector);
		lastnode = lastnode->next;
		len++;
	} else {
		node = new DenseVectorListElement<T>(vector);
		lastnode = node;
		iterator = node;
		len = 1;
		pos = 1;
	}
};

template <typename T>
void 
DenseVectorList<T>::remove(const int i){
	if(len == 1 && i == 1){
		delete(node);
		node = NULL;
		iterator = NULL;
		pos = 0;
		len = 0;
	} else {
		DenseVectorListElement<T> * help;
		if(i < pos){
			pos = 1;
			iterator = node;
		}
		if( i == 1){
			//remove first element
			iterator = node;
			pos = 1;
			node = node->next;
			delete(iterator);
			iterator = node;
			len--;
		} else{
			if(i == pos + 1){
				help = iterator->next;
				iterator->next = iterator->next->next;
				delete help;
				help = NULL;
				len --;
			}else {
				while(iterator->next != 0){
					iterator = iterator->next;
					pos++;
					if(i == pos+1) {
						if(i == len){
							delete(lastnode);
							iterator->next = NULL;
							lastnode = iterator;
							len--;
						} else {
							help = iterator->next;
							iterator->next = iterator->next->next;
							delete help;
							help = NULL;
							len --;
						}
					}
				}
			}
		}
	}
};

template <typename T>
void 
DenseVectorList<T>::empty(){
	int end = len;
	for(int i = 1; i <= end; ++i){
		this->remove(1);
	}
};


template <typename T>
int 
DenseVectorList<T>::length() const{
	return len;
};

template <typename T>
flens::DenseVector<flens::Array<T> > * 
DenseVectorList<T>::operator()(const int i){
	if(len > 0){
		if(i < pos){
			pos = 1;
			iterator = node;
		}
		if(i == pos) return &(iterator->elem);
		while(iterator->next != 0){
			iterator = iterator->next;
			pos++;
			if(i == pos) return &(iterator->elem);
		}
	}
	return NULL;
};

template <typename T>
std::ostream &
operator<<(std::ostream &out, DenseVectorList<T> &vl)
{
	out << "VL.length = " << vl.length() << std::endl;
	for(int i = 1; i <= vl.length(); ++i){
		out << "Vector "<<(i+1) << ":  " << *vl(i) ;
	}
	return out;
};

} //namespace htucker
