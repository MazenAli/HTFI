#ifndef HTUCKER_FLENSEXT_OUTPUTSPARSEMATRIX_H
#define HTUCKER_FLENSEXT_OUTPUTSPARSEMATRIX_H 1

#include <iostream>
#include <extensions/flens/flens.h>


namespace flens{
	
    template <typename T>
	std::ostream& operator<<(std::ostream& Stream, flens::SparseGeMatrix<flens::extensions::CRS<T> > &B)
	{
		Stream << B.numRows() << " x " << B.numCols() << ": " <<  std::endl;
		const flens::DenseVector<flens::Array<T> > & ref = (B.engine().values);
		const flens::DenseVector<flens::Array<T> > & col = (B.engine().columns);
		const flens::DenseVector<flens::Array<T> > & refrows = B.engine().rows ;
		int j = 1;
		for(int i = ref.firstIndex(); i <= ref.lastIndex(); ++i){
			
			if(j < refrows.length() && refrows(j+1) == i){
				j = j + 1;
			}
			Stream << "(" << j << ", " << col(i) << ") = " << ref(i) << std::endl;
		}
		
		return Stream;
	}
} //namespace flens

#endif //HTUCKER_FLENSEXT_OUTPUTSPARSEMATRIX_H
