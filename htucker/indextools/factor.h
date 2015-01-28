#ifndef HTUCKER_INDEXTOOLS_FACTOR_H
#define HTUCKER_INDEXTOOLS_FACTOR_H 1

namespace htucker{

	DimensionIndex factor(const int i){
		for(int j = 2; j < sqrt(i+0.0); ++j){
			if(i % j == 0){
				return factor(i/j).join(DimensionIndex(j,1));
			}
		}
		return DimensionIndex(i,1);
	};
}

#endif //HTUCKER_INDEXTOOLS_FACTOR_H 