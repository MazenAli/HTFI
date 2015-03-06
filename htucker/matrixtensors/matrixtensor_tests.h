#ifndef HTUCKER_MATRIXTENSORS_MATRIXTENSOR_TESTS_H
#define HTUCKER_MATRIXTENSORS_MATRIXTENSOR_TESTS_H 1

#include <htucker/matrixtensors/matrixtensor.h>


namespace htucker{


	void test_matrixtensor(){
		MatrixTensor<SparseGeMatrix<flens::extensions::CRS<double,flens::CRS_General> > > mt(2);

		SparseGeMatrix<flens::extensions::CRS<double, flens::CRS_General> > A(3,3);
		A(1,1) = 1;
		A(2,2) = 2;
		A(3,3) = 3;
		A.finalize();
		mt.setMatrix(1,A);
		A.~SparseGeMatrix<flens::extensions::CRS<double, flens::CRS_General> >();
		SparseGeMatrix<flens::extensions::CRS<double,flens::CRS_General> > B(2,2);
		B(1,2) = 1.5;
		B(2,1) = 3;
		B.finalize();
		/*cout << B << endl;
		mt.setMatrix(2,B);
		cout << mt.getMatrix(1) << endl;*/
		
	}


}

#endif // HTUCKER_MATRIXTENSORS_MATRIXTENSOR_TESTS_H
