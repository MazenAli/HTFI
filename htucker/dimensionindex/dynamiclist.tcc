#include <cassert>

template <typename T>
dynamiclist<T>::dynamiclist():count(0),root(NULL),last(NULL){}

template <typename T>
void
dynamiclist<T>::append(const T &val){
	if(root == NULL){
		root = new dynamiclistelement<T>(val);
		last = root;
	} else {
		last->next = new dynamiclistelement<T>(val);
		last->next->previous = last;
		last = last->next;
	}

	count++;
}

template <typename T>
T
dynamiclist<T>::operator[](const int i) const{
	dynamiclistelement<T> * pos = root;
	assert(i<count && i >= 0);
	for(int j = 0; j < i; ++j){
		pos = pos->next;
	}
	return pos->value;
}


template <typename T>
int
dynamiclist<T>::length() const{
	return count;
}
