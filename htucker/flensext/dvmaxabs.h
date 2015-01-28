#ifndef HTUCKER_FLENSEXT_DVMAXABS_H
#define HTUCKER_FLENSEXT_DVMAXABS_H 1

namespace flens{
	template <typename T>
	T dvmaxabs(const flens::DenseVector<flens::Array<T> > &vec, int & pos){
		T max = -1;
		int li = vec.lastIndex();
		for(int i = vec.firstIndex(); i<= li; ++i){
			double absval =  abs(vec(i));
			if(absval > max){
				max = absval;
				pos = i;
			}
		}
		return max;
	}
}

#endif // HTUCKER_FLENSEXT_DVMAXABS_H 