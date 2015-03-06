#ifndef LAWA_METHODS_HTUCKER_MATRIXTENSORS_HTUCKERCLOSURE_TESTS_H
#define LAWA_METHODS_HTUCKER_MATRIXTENSORS_HTUCKERCLOSURE_TESTS_H 1

#include<lawa/methods/htucker/matrixtensors/matrixtensor.h>

using namespace std;

namespace lawa{


void test_htuckerclosure(){
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
	cout << B << endl;
	mt.setMatrix(2,B);
	
	emptytensor et;
	tensor<double> empty(et,2,1,2);
	HTuckerTree<double,SVD> testtree(2,empty);

	auto mat = mt + mt;
	mat*testtree;
	
};


};

#endif // LAWA_METHODS_HTUCKER_MATRIXTENSORS_HTUCKERCLOSURE_TESTS_H
