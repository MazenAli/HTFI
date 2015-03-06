#ifndef HTUCKER_FLENSEXT_OUTPUTSPARSEMATRIX_H
#define HTUCKER_FLENSEXT_OUTPUTSPARSEMATRIX_H 1

#include <iostream>
#include <extensions/flens/flens.h>


namespace flens{
    using namespace std;
	
    template <typename T>
	ostream& operator<<(ostream& Stream, SparseGeMatrix<flens::extensions::CRS<T> > &B)
	{
		Stream << B.numRows() << " x " << B.numCols() << ": " <<  endl;
		const DenseVector<Array<T> > & ref = (B.engine().values);
		const DenseVector<Array<T> > & col = (B.engine().columns);
		const DenseVector<Array<T> > & refrows = B.engine().rows ;
		int j = 1;
		for(int i = ref.firstIndex(); i <= ref.lastIndex(); ++i){
			
			if(j < refrows.length() && refrows(j+1) == i){
				j = j + 1;
			}
			Stream << "(" << j << ", " << col(i) << ") = " << ref(i) << endl;
		}
		
		return Stream;
	}
} //namespace flens

#endif //HTUCKER_FLENSEXT_OUTPUTSPARSEMATRIX_H
