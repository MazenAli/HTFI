#ifndef LAWA_APPLICATIONS_CDO_CDO_H
#define LAWA_APPLICATIONS_CDO_CDO_H 1

#include <applications/cdo/kscdo.h>

using namespace std;
using namespace flens;
using namespace lawa;





int test_lhs(){
	DenseVector<Array<double> > ints(2);
	ints = 1.5, 0.75;

	kscdomatrix lhs(1,6,0.05,0.25,0.05,ints);
	lhs(lhs.maxvals);
};


int test_operator2tensor(){
	
	Basis<double,Orthogonal,Interval,Multi> basis(3);
	BSOperator1D<double,Basis<double,Orthogonal,Interval,Multi> > a(basis,0.05,0.1,0.25);
	WavOperator2Tensor<BSOperator1D<double,Basis<double,Orthogonal,Interval,Multi> > > tens(a,1,7);
	tensor<double> wavtens(tens,tens.minvals.length(),tens.minvals,tens.quadmaxi);
	HTuckerTree<double,SVD> WavTree(tens.minvals.length(),wavtens);
	WavTree.print();
	double eps = 0.000000000001;
	WavTree.approximate(eps);
	setUorB(WavTree);
	WavTree.print();
	cout << "here " << endl;
	double err = WavTree.Linfnorm(wavtens,tens.minvals,tens.maxvals,1000);
	cout << "error = " << err << endl;
	cout << "soweit ist schonmal alles gut gegangen" << endl;


};






double
rhs_f(double x)
{
    return 1.;
}




void test_problem(){
	int di[4] = {1,2};
	DimensionIndex dims(di,2);
	DimensionIndex mi(0,2);
	DimensionIndex ma(3,2);
	DimensionIndex run = mi;

	DenseVector<Array<double> > ints(2);
	ints = 1.5, 0.75;

	transitionintensities ti(ints);

	for(DimensionIndexIterator<IteratorXD> it = run.getIterator(dims,mi,ma); it.inRange(); it++){
		cout << "index = " << it.getIndex() << "  => " << ti(it.getIndex()) <<  endl;
	}

	double eps = 0.000000001;
	tensor<double> lambdatens(ti,2,0,3);
	HTuckerTree<double,SVD> lambdatree(2,lambdatens);
	lambdatree.approximate(eps);
	setUorB(lambdatree);
	cout << endl << endl;
	DimensionIndex onessize(2,2);
	HTuckerTree<double,SVD> lambdadiag = vec2diag(lambdatree*ones<double,SVD>(onessize));
	HTuckerTree<double,SVD> lambdatree2 = lambdatree - lambdadiag;
	HTuckerTree<double,SVD> lambdatree3 = reapproximate(lambdatree2,eps);
	//cout << "reapproximated" << endl;
	for(DimensionIndexIterator<IteratorXD> it = run.getIterator(dims,mi,ma); it.inRange(); it++){
		cout << "index = " << it.getIndex() << "  => " << lambdatree2.evaluate(it.getIndex()) <<  endl;
	}

//-------------------------------------------------------------------------------------------------------------------------------------

	int d = 2;          // (d,d_)-wavelets
    int j0 = 2;         // minimal level
    int J = 3;          // maximal level
    
    /// Basis initialization, using Dirichlet boundary conditions
    Basis<double,Orthogonal,Interval,Multi> basis(d);
    basis.enforceBoundaryCondition<DirichletBC>();
    
    /// Operator initialization
    BSOperator1D <double, Basis<double,Orthogonal,Interval,Multi> > a(basis, 0.05,0.07,0.25);

	Assembler1D <double, Basis<double,Orthogonal,Interval,Multi> > assembler(basis);
	SparseGeMatrix<flens::CRS<double ,flens::CRS_General> >   A = assembler.assembleStiffnessMatrix(a, J);


	MatrixTensor<SparseGeMatrix<flens::CRS<double ,flens::CRS_General> > > mt(1);
	mt.setMatrix(1,A);
	
//---------------------------------------------------------------------------------------------------------------------------------------
	DimensionIndex idleftmax(mt.dim());
	for(int i = 0; i < idleftmax.length(); ++i){
		idleftmax[i] = mt.getMatrix(i+i).numRows();
	}
	IdentityTensor idleft(idleftmax);

	IdentityTensor idright(ma);


	auto op = mat(lambdatree)*idleft + idright*mt;

// jetzt die rechte seite aufstellen...



};


#endif // LAWA_METHODS_HTUCKER_LINEARSOLVER_LINEARSOLVER_TESTS_H