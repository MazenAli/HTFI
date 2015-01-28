#ifndef LAWA_APPLICATIONS_CDO_KSCDO_H
#define LAWA_APPLICATIONS_CDO_KSCDO_H 1

#include <lawa/methods/htucker/dimensionindex/densevectorlist.h>
#include <lawa/methods/htucker/dimensionindex/dimensionindex.h>

using namespace std;
using namespace flens;

namespace lawa{

class kscdomatrix : public TensorFunction{
public:
	double alpha, beta, Tstep;
	int j0,J;
	DimensionIndex cumWavNumber;
	DenseVector<Array<double> > intensities;
	DimensionIndex minvals,maxvals;

	HTuckerTree<double,SVD> lambda;
	

	kscdomatrix(int _j0, int _J, double _alpha, double _beta, double _Tstep, DenseVector<Array<double> > _intensities);
	
	double operator() (DimensionIndex vals);
};


class transitionintensities : public TensorFunction{
public:
	DenseVector<Array<double> > intensities;
	DimensionIndex _min, _max, _inner; //inner berschreibt den umbruchpunkt, d.h. die eindimensionale Dimension der Spalten.
	DimensionIndex zeilendims,spaltendims,zeilendimscum,spaltendimscum; // Anzahl der Elemente pro Dimension...

	transitionintensities(DenseVector<Array<double> > _intensities);
	
	double operator() (DimensionIndex vals);
};

class transitionmatrix: public TensorFunction{
public:
	HTuckerTree<double,SVD> &t1, &t2;
	DimensionIndex _min,_max,_inner;
	DimensionIndex zeilendims, spaltendims, zeilendimscum, spaltendimscum;

	transitionmatrix(HTuckerTree<double,SVD> &_t1, HTuckerTree<double,SVD> &_t2);

	double operator() (DimensionIndex vals);
};


} // end namespace lawa

#include <applications/cdo/kscdo.tcc>

#endif // LAWA_APPLICATIONS_CDO_KSCDO_H