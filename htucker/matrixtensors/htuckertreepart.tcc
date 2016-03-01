namespace htucker{
	template <typename T>
	HTuckerTreePart<T>::HTuckerTreePart(const HTuckerTree<T> & tree, const int _min, const int _max):httree(tree),mindim(_min),maxdim(_max){}
}
