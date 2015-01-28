#ifndef HTUCKER_FLENSEXT_DVSORT_H
#define HTUCKER_FLENSEXT_DVSORT_H 1

namespace flens{
	template <typename VX, typename VI>
	void
	sort(DenseVector<VX> &x, DenseVector<VI> &rho){
		using std::swap;
		typedef typename DenseVector<VX>::IndexType  IndexType;
		
		const IndexType n = x.length();

		rho = DenseVector<VI>(n);
		for (IndexType i=1; i<=n; ++i) {
			rho(i) = i;
		}
		
		for (IndexType i=1; i<=n; ++i) {
			bool swaped = false;
			for (IndexType j=1; j<=n-i; ++j) {
				if (x(j)>x(j+1)) {
					swap(x(j),x(j+1));
					swap(rho(j),rho(j+1));
					swaped = true;
				}
			}
			if (!swaped) {
				break;
			}
		}
	}
}

#endif // HTUCKER_FLENSEXT_DVSORT_H 