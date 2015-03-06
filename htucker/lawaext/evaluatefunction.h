#ifndef HTUCKER_LAWAEXT_EVALUATEFUNCTION_H
#define HTUCKER_LAWAEXT_EVALUATEFUNCTION_H 1

namespace lawa{
	template <typename _T, FunctionSide _Side, DomainType _Domain, Construction _Cons>
	void
	evaluatefunction(const BasisFunction<_T,_Side,_Domain,_Cons> * BF,const  int j, const long k, const unsigned short deriv, const int steps, flens::DenseVector<flens::Array<_T> > &x){
		assert(steps>=2);
		flens::DenseVector<flens::Array<_T> > y(steps+1);
		flens::Support<_T> supp = BF->support(j,k);
		_T dist = supp.l2 - supp.l1;
		for(int i = 1; i<= steps+1; ++i){
			y(i) = BF->operator()(supp.l1+(i-1)*dist/(steps+0.0),j,k,deriv);
		}
		x = y;
	}

}

#endif // HTUCKER_LAWAEXT_EVALUATEFUNCTION_H
