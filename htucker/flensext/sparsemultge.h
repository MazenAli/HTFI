#ifndef HTUCKER_FLENSEXT_SPARSEMULTGE_H
#define HTUCKER_FLENSEXT_SPARSEMULTGE_H 1


namespace flens{

using namespace cxxblas;

template <typename T> 
GeMatrix<FullStorage<T,ColMajor> >
operator*(SparseGeMatrix<flens::extensions::CRS<T,flens::CRS_General> > & mat, GeMatrix<FullStorage<T,ColMajor> > & gemat){
	assert(mat.numCols() == gemat.numRows());
	GeMatrix<FullStorage<T,ColMajor> > ret(mat.numRows(),gemat.numCols());
	flens::extensions::CRS<T,flens::CRS_General> eng = mat.engine();
	T val = 0;
	for(int i = 1; i <= mat.numRows(); ++i){
		for(int j = 1; j<= gemat.numCols(); ++j){
			val = 0;
			for(int k = eng.rows(i); k < eng.rows(i+1); ++k){
				val += eng.values(k)*gemat(eng.columns(k),j);
			}
			ret(i,j) = val;
		}
	}
	return ret;
};

}

#endif //HTUCKER_FLENSEXT_SPARSEMULTGE
